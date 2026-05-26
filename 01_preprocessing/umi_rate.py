import argparse
import gzip
import os
from collections import Counter

import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import pysam


def process_gene_data(barcode, umi, gene, gene_dict, gene_umi_dict,
                      bc_gene_dict, bc_gene_umi_dict):
    gene_name = gene[:-3]  # strip '--C' or '--T' suffix (always 3 chars)
    gene_dict.add(gene_name)
    gene_umi_dict.setdefault(gene_name, set()).add(f'{barcode}#{umi}')
    bc_gene_dict.setdefault(barcode, set()).add(gene_name)
    bc_gene_umi_dict.setdefault(barcode, set()).add(f'{gene_name}#{umi}')


def create_expression_matrix(bc_gene_umi: dict) -> pd.DataFrame:
    """Build a gene × barcode UMI count matrix from bc_gene_umi dict.

    bc_gene_umi maps barcode → set of 'gene_name#umi' strings.
    Each unique 'gene_name#umi' pair counts as 1 UMI for that gene.
    """
    records = []
    for barcode, entries in bc_gene_umi.items():
        gene_counts = Counter(entry.split('#')[0] for entry in entries)
        for gene_name, count in gene_counts.items():
            records.append((gene_name, barcode, count))

    if not records:
        return pd.DataFrame()

    df = pd.DataFrame(records, columns=['gene', 'barcode', 'count'])
    matrix = df.pivot_table(index='gene', columns='barcode',
                            values='count', aggfunc='sum', fill_value=0)
    matrix.columns.name = None
    matrix.index.name = None
    return matrix


def calculate_umi_rate(bc_new_gene_umi: dict, bc_old_gene_umi: dict,
                       cell_barcodes: list, mode: str) -> list:
    """Return per-barcode fraction of new (mode='new') or old (mode='old') UMIs."""
    rates = []
    for bc in cell_barcodes:
        new_umi   = len(bc_new_gene_umi.get(bc, set()))
        old_umi   = len(bc_old_gene_umi.get(bc, set()))
        total_umi = new_umi + old_umi
        if total_umi == 0:
            rates.append(0.0)
        elif mode == 'new':
            rates.append(new_umi / total_umi)
        else:
            rates.append(old_umi / total_umi)
    return rates


def main(args):
    result_path = args.result_path
    sample_path = args.sample_path

    new_gene, old_gene = set(), set()
    new_gene_umi, old_gene_umi = {}, {}
    bc_new_gene, bc_old_gene = {}, {}
    bc_new_gene_umi, bc_old_gene_umi = {}, {}

    barcode_path = os.path.join(
        sample_path, 'outs', 'filtered_feature_bc_matrix', 'barcodes.tsv.gz')
    with gzip.open(barcode_path, 'rt') as f:
        cell_barcodes = [line.strip() for line in f]  # ordered list
    cell_barcode_set = set(cell_barcodes)

    with pysam.AlignmentFile(args.bam_path, 'rb') as f:
        for r in f.fetch():
            barcode = r.get_tag('CB').replace('_', '')
            geneid  = r.get_tag('GX')
            umi     = r.get_tag('UB')

            if barcode == '-' or geneid == '-' or umi == '-':
                continue
            if barcode not in cell_barcode_set:
                continue

            try:
                gene = r.get_tag('GN')
            except KeyError:
                continue

            if gene.endswith('--T'):
                process_gene_data(barcode, umi, gene,
                                  old_gene, old_gene_umi,
                                  bc_old_gene, bc_old_gene_umi)
            elif gene.endswith('--C'):
                process_gene_data(barcode, umi, gene,
                                  new_gene, new_gene_umi,
                                  bc_new_gene, bc_new_gene_umi)

    new_gene_matrix = create_expression_matrix(bc_new_gene_umi)
    old_gene_matrix = create_expression_matrix(bc_old_gene_umi)

    os.makedirs(result_path, exist_ok=True)
    new_gene_matrix.to_csv(os.path.join(result_path, 'new_gene_umi.csv'))
    old_gene_matrix.to_csv(os.path.join(result_path, 'old_gene_umi.csv'))

    combined_matrix = new_gene_matrix.add(old_gene_matrix, fill_value=0)
    combined_matrix.to_csv(os.path.join(result_path, 'all_gene_umi.csv'))

    new_gene_umi_rate = calculate_umi_rate(
        bc_new_gene_umi, bc_old_gene_umi, cell_barcodes, mode='new')
    old_gene_umi_rate = calculate_umi_rate(
        bc_new_gene_umi, bc_old_gene_umi, cell_barcodes, mode='old')

    umi_rate_df = pd.DataFrame({
        'New Gene UMI Rate': new_gene_umi_rate,
        'Old Gene UMI Rate': old_gene_umi_rate,
    }, index=cell_barcodes)
    umi_rate_df.to_csv(os.path.join(result_path, 'gene_umi_rates_with_index.csv'))

    plt.figure(figsize=(4, 10))
    box = plt.boxplot(
        [new_gene_umi_rate, old_gene_umi_rate],
        patch_artist=True,
        medianprops={'color': 'black', 'linewidth': 2},
    )
    box['boxes'][0].set(facecolor='red',  alpha=0.5)
    box['boxes'][1].set(facecolor='blue', alpha=0.5)
    plt.xticks([1, 2], ['New Gene', 'Old Gene'], fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim(0, 1)
    plt.ylabel('Rate of Gene UMI Total Count per Cell', fontsize=14)
    plt.savefig(os.path.join(result_path, 'bc_new_old_gene_umi_rate_box.pdf'),
                bbox_inches='tight')
    plt.savefig(os.path.join(result_path, 'bc_new_old_gene_umi_rate_box.png'),
                bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Process BAM files for gene and UMI data analysis.')
    parser.add_argument('--sample_path',
                        help='Base folder containing all samples')
    parser.add_argument('--bam_path',
                        help='Path to the tagged BAM file')
    parser.add_argument('--result_path',
                        help='Output directory for results')
    main(parser.parse_args())
