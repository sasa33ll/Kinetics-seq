import argparse
import os
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

import matplotlib.pyplot as plt
import ngs_tools as ngs
import numpy as np
import pandas as pd
import pysam


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def complement(base: str) -> str:
    return _COMPLEMENT.get(base, base)


# ---------------------------------------------------------------------------
# Per-BAM processing (runs in worker processes)
# ---------------------------------------------------------------------------

def process_bam_file(bam_file_path: str, barcodes: set, n_threads: int):
    """Process one BAM file; return (conversions, alignments, conversion_counts)."""
    conversions = []
    alignments  = []
    conversion_counts: Counter = Counter()

    with pysam.AlignmentFile(bam_file_path, 'rb', threads=n_threads) as bam:
        # Use index statistics to skip empty contigs without re-scanning reads
        valid_contigs = {
            stat.contig
            for stat in bam.get_index_statistics()
            if stat.mapped > 0
        }

        for contig in valid_contigs:
            for record in bam.fetch(contig):
                barcode = record.get_tag('CB').replace('_', '')
                gene_id = record.get_tag('GX')
                umi     = record.get_tag('UB')
                strand  = record.get_tag('ST')

                if barcode == '-' or gene_id == '-' or umi == '-':
                    continue
                if barcode not in barcodes:
                    continue

                ref_base_counts = Counter(record.get_reference_sequence().upper())
                read_sequence   = record.query_sequence.upper()
                mutations: dict = {}

                for read_idx, ref_idx, ref_base in record.get_aligned_pairs(
                        matches_only=True, with_seq=True):
                    read_base = read_sequence[read_idx]
                    ref_base  = ref_base.upper()

                    if 'N' in (ref_base, read_base):
                        continue

                    if strand == '+':
                        conversion = f'{ref_base}{read_base}'
                    else:
                        conversion = f'{complement(ref_base)}{complement(read_base)}'

                    conversion_counts[conversion] += 1

                    if ref_base != read_base:
                        mutations[ref_idx] = {
                            'read_id':    record.query_name,
                            'contig':     contig,
                            'position':   ref_idx,
                            'conversion': conversion,
                        }

                conversions.extend(mutations.values())
                alignments.append({
                    'read_id':  record.query_name,
                    'barcode':  barcode,
                    'umi':      umi,
                    'gene_id':  gene_id,
                    'A': ref_base_counts.get('A', 0),
                    'C': ref_base_counts.get('C', 0),
                    'G': ref_base_counts.get('G', 0),
                    'T': ref_base_counts.get('T', 0),
                    'status':   'unassigned',
                    'assigned': bool(gene_id),
                })

    return conversions, alignments, conversion_counts


# ---------------------------------------------------------------------------
# Temp file cleanup
# ---------------------------------------------------------------------------

def delete_temp_bam_files(temp_dir: str):
    if not os.path.isdir(temp_dir):
        return
    for filename in os.listdir(temp_dir):
        if filename.startswith('temp_') and (
                filename.endswith('.bam') or filename.endswith('.bam.bai')):
            file_path = os.path.join(temp_dir, filename)
            if os.path.exists(file_path):
                os.remove(file_path)


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def main(bam_file: str, barcodes: set, n_threads: int,
         max_workers: int, out_dir: str, gtf_file: str):
    all_conversions: list  = []
    all_alignments:  list  = []
    conversion_counts: Counter = Counter()

    os.makedirs(out_dir, exist_ok=True)
    temp_dir = os.path.join(out_dir, 'temp')

    gene_infos, _ = ngs.gtf.genes_and_transcripts_from_gtf(gtf_file, use_version=False)

    split_bam_paths = [
        value[0]
        for value in ngs.bam.split_bam(
            bam_file, temp_dir, n=64,
            n_threads=n_threads, show_progress=True
        ).values()
    ]

    for bam_path in split_bam_paths:
        pysam.index(bam_path, f'{bam_path}.bai', '-@', str(n_threads))

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(process_bam_file, bam_path, barcodes, n_threads)
            for bam_path in split_bam_paths
        ]
        for future in futures:
            local_conversions, local_alignments, local_counts = future.result()
            all_conversions.extend(local_conversions)
            all_alignments.extend(local_alignments)
            conversion_counts.update(local_counts)

    delete_temp_bam_files(temp_dir)
    return all_conversions, all_alignments, conversion_counts


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Process BAM files and calculate conversion statistics.')
    parser.add_argument('--sample_path', help='Sample directory (contains outs/)')
    parser.add_argument('--bam_file',    help='Path to the BAM file')
    parser.add_argument('--gtf_file',    help='Path to the GTF file')
    parser.add_argument('--out_dir',     help='Output directory')
    parser.add_argument('--max_workers', type=int, default=4,
                        help='Maximum number of worker processes')
    parser.add_argument('--n_threads',   type=int, default=4,
                        help='Number of threads per BAM operation')
    args = parser.parse_args()

    barcodes: set = set(
        pd.read_csv(
            os.path.join(args.sample_path,
                         'outs/filtered_feature_bc_matrix/barcodes.tsv.gz'),
            header=None,
        )[0].tolist()
    )

    all_conversions, all_alignments, conversion_counts = main(
        args.bam_file, barcodes, args.n_threads,
        args.max_workers, args.out_dir, args.gtf_file,
    )

    pd.DataFrame(all_conversions).to_csv(
        os.path.join(args.out_dir, 'all_conversions_list.csv'), index=False)
    pd.DataFrame(all_alignments).to_csv(
        os.path.join(args.out_dir, 'all_alignments_list.csv'), index=False)

    for conversion, count in sorted(conversion_counts.items()):
        print(f'Conversion: {conversion}, Count: {count}')

    # Conversion proportions (exclude identity conversions)
    start_base_totals: Counter = Counter()
    for conv, count in conversion_counts.items():
        start_base_totals[conv[0]] += count

    identity = {'AA', 'CC', 'GG', 'TT'}
    conversion_proportions = {
        conv: count / start_base_totals[conv[0]]
        for conv, count in conversion_counts.items()
        if conv not in identity
    }

    labels      = sorted(conversion_proportions)
    proportions = [conversion_proportions[l] for l in labels]
    colors      = plt.cm.rainbow(np.linspace(0, 1, len(labels)))

    plt.figure(figsize=(8, 4))
    plt.bar(labels, proportions, color=colors)
    plt.xlabel('Conversion Type')
    plt.ylabel('Proportion')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, 'conversion_proportions.pdf'),
                bbox_inches='tight')
    plt.close()

    pd.DataFrame(
        list(conversion_proportions.items()),
        columns=['Conversion_Type', 'Proportion'],
    ).to_csv(os.path.join(args.out_dir, 'conversion_proportions.csv'), index=False)

    pd.DataFrame(
        list(conversion_counts.items()),
        columns=['Conversion_Type', 'Count'],
    ).to_csv(os.path.join(args.out_dir, 'conversion_counts.csv'), index=False)
