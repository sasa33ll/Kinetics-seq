import argparse
import logging
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from typing import Dict, List

import numpy as np
import pandas as pd
import pysam

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

Base = ('A', 'T', 'C', 'G')
Base_index = {b: i for i, b in enumerate(Base)}


# ---------------------------------------------------------------------------
# GTF helpers
# ---------------------------------------------------------------------------

def read_gtf(file_path: str) -> pd.DataFrame:
    columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    return pd.read_csv(file_path, sep='\t', comment='#', header=None, names=columns)


def parse_attributes(attr_str: str) -> Dict[str, str]:
    attributes = {}
    for item in attr_str.split(';'):
        item = item.strip()
        if item:
            key, value = item.split(' ', 1)
            attributes[key] = value.strip('"')
    return attributes


def build_gene_infos(gtf_file: str) -> dict:
    gtf_data = read_gtf(gtf_file)
    gene_infos = {}
    for _, row in gtf_data[gtf_data["feature"] == "gene"].iterrows():
        attributes = parse_attributes(row["attribute"])
        gene_id = attributes.get("gene_id")
        if gene_id:
            gene_infos[gene_id] = {
                "gene_name": attributes.get("gene_name"),
                "gene_type": attributes.get("gene_type"),
                "chromosome": row["seqname"],
                "start": row["start"],
                "end": row["end"],
                "strand": row["strand"],
                "segment": (row["start"], row["end"]),
            }
    return gene_infos


# ---------------------------------------------------------------------------
# MD tag parser — single module-level definition shared by all callers
# ---------------------------------------------------------------------------

def parse_md_tag(md_tag: str) -> List[tuple]:
    """Return list of (genome_position, ref_base) for each mismatch in MD tag."""
    mismatches = []
    position = 0
    number = ''
    i = 0
    while i < len(md_tag):
        char = md_tag[i]
        if char.isdigit():
            number += char
        else:
            if number:
                position += int(number)
                number = ''
            if char == '^':
                i += 1
                while i < len(md_tag) and md_tag[i] in 'ACGT':
                    i += 1
                continue
            mismatches.append((position, char))
            position += 1
        i += 1
    return mismatches


# ---------------------------------------------------------------------------
# Consensus building
# ---------------------------------------------------------------------------

def process_batch(batch, gene_infos: dict, header_dict: dict,
                  temp_out_dir: str, batch_num: int):
    temp_out = os.path.join(temp_out_dir, f"batch_{batch_num}.bam")
    with pysam.AlignmentFile(
            temp_out, 'wb',
            header=pysam.AlignmentHeader.from_dict(header_dict)) as out:

        for gene, barcode_data in batch:
            gene_info = gene_infos[gene]
            strand = gene_info['strand']

            for barcode, umi_data in barcode_data.items():
                for umi, reads in umi_data.items():

                    # ---- single-read shortcut ----
                    if len(reads) == 1:
                        reads[0].set_tag('ST', strand)
                        out.write(reads[0])
                        continue

                    # ---- multi-read consensus ----
                    tags = {
                        'CB': barcode,
                        'UB': umi,
                        'GX': gene,
                        'GN': gene_info['gene_name'],
                        'ST': strand,
                    }

                    left_pos  = min(r.reference_start for r in reads)
                    right_pos = max(r.reference_end   for r in reads)
                    length    = right_pos - left_pos

                    seq = np.zeros((length, 4), dtype=np.float32)  # base counts
                    qua = np.zeros((length, 4), dtype=np.float32)  # summed Phred scores
                    ref = np.full(length, -1, dtype=np.int8)

                    for read in reads:
                        read_seq       = read.query_sequence.upper()
                        read_qualities = read.query_qualities
                        for read_i, genome_i, genome_base in read.get_aligned_pairs(
                                matches_only=False, with_seq=True):
                            if genome_i is None or genome_base is None:
                                continue
                            genome_base = genome_base.upper()
                            if genome_base == 'N':
                                continue
                            i = genome_i - left_pos
                            if read_i is None:
                                if ref[i] == -1:
                                    ref[i] = Base_index[genome_base]
                                continue
                            read_base = read_seq[read_i]
                            if read_base == 'N':
                                continue
                            if ref[i] == -1:
                                ref[i] = Base_index[genome_base]
                            b_idx = Base_index[read_base]
                            seq[i, b_idx] += 1
                            qua[i, b_idx] += read_qualities[read_i]

                    consensus_length = int((seq > 0).any(axis=1).sum())
                    consensus      = np.zeros(consensus_length, dtype=np.uint8)
                    consensus_qual = np.zeros(consensus_length, dtype=np.uint8)

                    cigar = []
                    last_cigar_op = None
                    cigar_n = 0
                    md = []
                    md_n = 0
                    md_zero = True
                    md_del = False
                    nm = 0
                    consensus_i = 0

                    for i in range(length):
                        ref_base = ref[i]
                        cigar_op = 'N'

                        if ref_base >= 0:
                            seq_base = seq[i]

                            if (seq_base == 0).all():
                                # deletion
                                cigar_op = 'D'
                                if md_n > 0 or md_zero:
                                    md.append(str(md_n))
                                    md_n = 0
                                if not md_del:
                                    md.append('^')
                                    md.append(Base[ref_base])
                                    md_del = True
                            else:
                                md_del = False
                                max_count = seq_base.max()
                                top_bases = (seq_base == max_count).nonzero()[0]
                                # quality-weighted tie-breaking
                                if len(top_bases) > 1:
                                    base = int(top_bases[np.argmax(qua[i, top_bases])])
                                else:
                                    base = int(top_bases[0])

                                cigar_op = 'M'
                                if ref_base == base:
                                    md_n += 1
                                    md_zero = False
                                else:
                                    if md_n > 0 or md_zero:
                                        md.append(str(md_n))
                                        md_n = 0
                                    md.append(Base[ref_base])
                                    md_zero = True
                                    nm += 1

                                consensus[consensus_i] = base
                                # clamp summed quality to uint8 max
                                consensus_qual[consensus_i] = min(int(qua[i, base]), 255)
                                consensus_i += 1

                        if cigar_op == last_cigar_op:
                            cigar_n += 1
                        else:
                            if last_cigar_op:
                                cigar.append(f"{cigar_n}{last_cigar_op}")
                            last_cigar_op = cigar_op
                            cigar_n = 1

                    md.append(str(md_n))
                    cigar.append(f"{cigar_n}{last_cigar_op}")

                    header = pysam.AlignmentHeader.from_dict(header_dict)
                    al = pysam.AlignedSegment(header)
                    al.query_name      = sha256(
                        ''.join(r.query_name for r in reads).encode('utf-8')
                    ).hexdigest()
                    al.query_sequence  = ''.join(Base[b] for b in consensus)
                    al.query_qualities = consensus_qual
                    al.reference_name  = reads[0].reference_name
                    al.reference_id    = reads[0].reference_id
                    al.reference_start = left_pos
                    al.mapping_quality = 255
                    al.cigarstring     = ''.join(cigar)

                    tags.update({'MD': ''.join(md), 'NM': nm})
                    al.set_tags(list(tags.items()))

                    al.is_unmapped      = False
                    al.is_paired        = False
                    al.is_duplicate     = False
                    al.is_qcfail        = False
                    al.is_secondary     = False
                    al.is_supplementary = False

                    out.write(al)


# ---------------------------------------------------------------------------
# BAM merge / sort / index / cleanup helpers
# ---------------------------------------------------------------------------

def merge_bams(output_file: str, temp_out_dir: str, num_batches: int):
    first_bam = os.path.join(temp_out_dir, "batch_0.bam")
    with pysam.AlignmentFile(first_bam, 'rb') as template_bam:
        with pysam.AlignmentFile(output_file, 'wb',
                                 header=template_bam.header) as merged_bam:
            for i in range(num_batches):
                temp_bam = os.path.join(temp_out_dir, f"batch_{i}.bam")
                with pysam.AlignmentFile(temp_bam, 'rb') as tmp:
                    for read in tmp:
                        merged_bam.write(read)


def create_bam_index(bam_file: str):
    try:
        subprocess.check_call(['samtools', 'index', bam_file])
        logging.info("Index for %s created successfully.", bam_file)
    except subprocess.CalledProcessError as e:
        logging.error("Error creating index for %s: %s", bam_file, e)


def delete_temp_bam_files(temp_out_dir: str, num_batches: int):
    for i in range(num_batches):
        temp_bam = os.path.join(temp_out_dir, f"batch_{i}.bam")
        if os.path.exists(temp_bam):
            os.remove(temp_bam)
            logging.info("Deleted temporary BAM file: %s", temp_bam)


def run_samtools_calmd(bam_file: str, ref_fa: str):
    out_path = f"{bam_file}_calmd.bam"
    try:
        with open(out_path, 'wb') as out_bam:
            process = subprocess.Popen(
                ['samtools', 'calmd', '-AEur', bam_file, ref_fa],
                stdout=out_bam,
                stderr=subprocess.PIPE,
            )
            _, stderr = process.communicate()
        if process.returncode != 0:
            logging.error("Error during samtools calmd: %s", stderr.decode())
        else:
            logging.info("Samtools calmd completed for %s.", bam_file)
    except Exception as e:
        logging.error("Error: %s", e)


# ---------------------------------------------------------------------------
# Tag labelling (new vs old RNA)
# ---------------------------------------------------------------------------

def process_tags(input_bam: str, output_bam: str):
    with pysam.AlignmentFile(input_bam, 'rb') as infile, \
         pysam.AlignmentFile(output_bam, 'wb', header=infile.header) as outfile:

        for read in infile.fetch():
            md        = read.get_tag('MD')
            strand    = read.get_tag('ST')
            gene_name = read.get_tag('GN')

            if md.isdigit():
                read.set_tag('GN', gene_name + '--T')
                outfile.write(read)
                continue

            if strand == '+' and 'T' not in md:
                read.set_tag('GN', gene_name + '--T')
                outfile.write(read)
                continue

            if strand == '-' and 'A' not in md:
                read.set_tag('GN', gene_name + '--T')
                outfile.write(read)
                continue

            label = '--T'
            if strand == '+':
                for position, ref_base in parse_md_tag(md):
                    if read.query_sequence[position] == 'C' and ref_base == 'T':
                        label = '--C'
                        break
            elif strand == '-':
                for position, ref_base in parse_md_tag(md):
                    if read.query_sequence[position] == 'G' and ref_base == 'A':
                        label = '--C'
                        break

            read.set_tag('GN', gene_name + label)
            outfile.write(read)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Process BAM files and GTF annotations.")
    parser.add_argument('--bam_path',     required=True,
                        help='Path to the input BAM file.')
    parser.add_argument('--gtf_file',     type=str,
                        default='Homo_sapiens.GRCh38.99.gtf',
                        help='Path to the GTF file.')
    parser.add_argument('--temp_out_dir', required=True,
                        help='Temporary output directory for BAM files.')
    parser.add_argument('--ref_fa',       required=True,
                        help='Path to the reference genome FASTA.')
    args = parser.parse_args()

    bam_path     = args.bam_path
    ref_fa       = args.ref_fa
    gtf_file     = args.gtf_file
    temp_out_dir = args.temp_out_dir

    # Create output directory before any writes
    os.makedirs(temp_out_dir, exist_ok=True)

    run_samtools_calmd(bam_path, ref_fa)
    bam_path = f"{bam_path}_calmd.bam"
    create_bam_index(bam_path)

    temp_out   = os.path.join(temp_out_dir,
                              'Aligned.sortedByCoord.UniqueGene.consensus.bam')
    sorted_bam = os.path.join(temp_out_dir,
                              'Aligned.sortedByCoord.UniqueGene.consensus.sorted.bam')
    tagged_bam = os.path.join(temp_out_dir,
                              'Aligned.sortedByCoord.UniqueGene.consensus.sorted.bam.tagST.bam')

    gene_infos = build_gene_infos(gtf_file)

    # Single BAM open: collect header and reads in one pass
    barcode_umi_groups: dict = {}
    with pysam.AlignmentFile(bam_path, 'rb') as f:
        header_dict = f.header.to_dict()
        for read in f.fetch():
            if read.is_unmapped or read.is_secondary:
                continue
            barcode = read.get_tag('CB')
            umi     = read.get_tag('UB')
            genes   = read.get_tag('GX')
            if barcode == '-' or umi == '-' or genes == '-':
                continue
            (barcode_umi_groups
             .setdefault(genes, {})
             .setdefault(barcode, {})
             .setdefault(umi, [])
             .append(read))

    batches = [
        list(barcode_umi_groups.items())[i:i + 2000]
        for i in range(0, len(barcode_umi_groups), 2000)
    ]

    with ThreadPoolExecutor(max_workers=50) as executor:
        futures = [
            executor.submit(process_batch, batch, gene_infos,
                            header_dict, temp_out_dir, i)
            for i, batch in enumerate(batches)
        ]
        for future in futures:
            future.result()

    merge_bams(temp_out, temp_out_dir, len(batches))
    create_bam_index(temp_out)
    pysam.sort("-o", sorted_bam, temp_out)
    create_bam_index(sorted_bam)
    logging.info("Final sorted BAM file: %s", sorted_bam)

    process_tags(sorted_bam, tagged_bam)
    create_bam_index(tagged_bam)
    logging.info("Tagged BAM file: %s", tagged_bam)

    delete_temp_bam_files(temp_out_dir, len(batches))


if __name__ == "__main__":
    main()
