import argparse
import logging
import os

import matplotlib.pyplot as plt
import ngs_tools as ngs
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import binom

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

CONVERSION_COLUMNS = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG']
BASE_COLUMNS       = ['A', 'C', 'G', 'T']
# Conversions excluded from background rate estimation (T-to-X are the signal)
EXCLUDE_FROM_PE    = {'TC', 'TA', 'TG'}


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_data(mutations_file, reads_file, group_info_file, convs_file):
    mutations_df = pd.read_csv(mutations_file)
    reads_df     = pd.read_csv(reads_file)
    group_meta   = pd.read_csv(group_info_file)
    conversions  = ngs.utils.read_pickle(convs_file)
    return mutations_df, reads_df, group_meta, conversions


# ---------------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------------

def preprocess_data(mutations_df, reads_df):
    count_df = (
        mutations_df
        .groupby(['read_id', 'conversion'])
        .size()
        .reset_index(name='count')
    )
    pivot_df = (
        count_df
        .pivot_table(index='read_id', columns='conversion',
                     values='count', fill_value=0)
        .reset_index()
    )
    merged_df = pd.merge(reads_df, pivot_df, on='read_id', how='left').fillna(0)

    conv_sum = merged_df[CONVERSION_COLUMNS].sum(axis=1)
    merged_df = (
        merged_df
        .assign(conversion_sum=conv_sum, has_conversions=conv_sum > 0)
        .sort_values(['has_conversions', 'conversion_sum'])
        .drop(columns=['conversion_sum', 'has_conversions'])
    )
    return merged_df


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate_conversions(df_sorted, conversions, all_conversions):
    results = []
    for convs in conversions:
        convs       = sorted(convs)
        other_convs = list(set(all_conversions) - set(convs))
        bases       = list({c[0] for c in convs})

        df_filtered = (
            df_sorted[(df_sorted[other_convs] == 0).all(axis=1)]
            if other_convs else df_sorted
        )
        df_combined = df_filtered[['barcode', 'gene_id']].copy()
        df_combined['conversion'] = df_filtered[list(convs)].sum(axis=1).astype(int)
        df_combined['base']       = df_filtered[bases].sum(axis=1)
        df_combined = (
            df_combined
            .groupby(['barcode', 'gene_id', 'conversion', 'base'], sort=False)
            .size()
            .reset_index(name='count')
        )
        results.append(df_combined)
    return results


# ---------------------------------------------------------------------------
# Group info merge
# ---------------------------------------------------------------------------

def merge_group_info(df_sorted, df_combined_list, group_meta):
    # Rename on a copy to avoid mutating the caller's DataFrame
    gm = group_meta.copy()
    gm.columns = ['barcode', 'group']

    df_combined = pd.concat(df_combined_list, ignore_index=True)
    df_sorted   = pd.merge(df_sorted,   gm, on='barcode', how='left')
    df_combined = pd.merge(df_combined, gm, on='barcode', how='left')
    return df_sorted, df_combined


# ---------------------------------------------------------------------------
# Background conversion rate (pe / q)
# ---------------------------------------------------------------------------

def calculate_mean_conversion_rate(df_sorted, group_meta, convs):
    """Compute per-group mean background conversion rate, excluding T-to-X conversions."""
    bg_convs = [c for c in convs if c not in EXCLUDE_FROM_PE]
    group_pe_values = {}

    for group_name, group_data in df_sorted.groupby('group'):
        group_sum = group_data.sum(numeric_only=True).astype(float)

        # Vectorised normalisation — divide each conversion by its reference base total
        # without modifying group_sum in place across iterations
        rates = pd.Series(
            {c: group_sum[c] / group_sum[c[0]] for c in bg_convs if group_sum[c[0]] > 0}
        )
        group_pe_values[group_name] = rates.mean() if len(rates) > 0 else 0.0

    return pd.DataFrame(group_pe_values.items(), columns=['Group', 'Mean Conversion Rate'])


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------

def sample_data(filtered_df, subset_size):
    sampled = []
    for group, group_df in filtered_df.groupby('group'):
        if len(group_df) < subset_size:
            logging.warning(
                "Only %d samples available for group '%s'; using all (requested %d).",
                len(group_df), group, subset_size,
            )
            sampled.append(group_df)
        else:
            sampled.append(group_df.sample(n=subset_size, random_state=1))
    return pd.concat(sampled, ignore_index=True)


# ---------------------------------------------------------------------------
# Likelihood functions
# ---------------------------------------------------------------------------

def LL(params, data):
    """Negative log-likelihood for two-component binomial mixture (p, q, pi)."""
    p, q, pi = params
    t1 = binom.pmf(data[:, 1], data[:, 0], p)
    t2 = binom.pmf(data[:, 1], data[:, 0], q)
    l  = np.clip(pi * t1 + (1 - pi) * t2, 1e-10, 1.0)
    return -np.sum(np.log(l))


def LL_2(pi, data, p, q):
    """Negative log-likelihood for fixed-p-q mixture, optimising pi only."""
    t1 = binom.pmf(data[:, 1], data[:, 0], p)
    t2 = binom.pmf(data[:, 1], data[:, 0], q)
    l  = np.clip(pi * t1 + (1 - pi) * t2, 1e-10, 1.0)
    return -np.sum(np.log(l))


# ---------------------------------------------------------------------------
# p / q estimation
# ---------------------------------------------------------------------------

def estimate_pq_values(df_combined_group_data, group_pe_values):
    """
    For each group, estimate the new-RNA conversion rate p by jointly
    optimising (p, pi) while holding q fixed at the background rate.
    """
    results = []
    for _, row in group_pe_values.iterrows():
        group  = row['Group']
        q      = row['Mean Conversion Rate']
        pqdata = (
            df_combined_group_data[df_combined_group_data['group'] == group]
            [['base', 'conversion']]
            .values
        )
        if len(pqdata) == 0:
            continue

        r_mtx = []
        for _ in range(100):
            init_p  = np.random.uniform(0.01, 0.99)
            init_pi = np.random.uniform(0.01, 0.99)
            # Jointly optimise p (new-RNA rate) and pi (mixing weight)
            result = minimize(
                lambda x: LL([x[0], q, x[1]], pqdata),
                x0=[init_p, init_pi],
                bounds=[(0.01, 0.99), (0.01, 0.99)],
            )
            if result.success and 0 <= result.x[1] <= 1:
                r_mtx.append([result.x[0], result.x[1], result.fun])

        if not r_mtx:
            continue

        r_mtx_df = pd.DataFrame(r_mtx, columns=['p', 'pi', 'nll'])
        best     = r_mtx_df.loc[r_mtx_df['nll'].idxmin()]
        results.append({'group': group, 'q_value': q, 'optimal_p': best['p']})

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Per-gene theta estimation
# ---------------------------------------------------------------------------

def calculate_theta_for_gene(results_df, genedata):
    theta_records = []
    for _, row in results_df.iterrows():
        group      = row['group']
        q_value    = row['q_value']
        optimal_p  = row['optimal_p']

        for gene_id, group_data in genedata[genedata['group'] == group].groupby('gene_id'):
            if len(group_data) < 100:
                continue
            data = group_data[['base', 'conversion']].values
            result = minimize_scalar(
                lambda x: LL_2(x, data, optimal_p, q_value),
                bounds=(0.01, 0.99),
                method='bounded',
            )
            theta = result.x if optimal_p >= q_value else 1.0 - result.x
            theta_records.append({'gene_id': gene_id, 'theta': theta, 'group': group})

    return pd.DataFrame(theta_records).set_index('gene_id')


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

def plot_boxplot_of_theta(all_theta_for_gene, out_dir):
    plt.figure(figsize=(6, 4))
    sns.boxplot(
        x='group', y='theta', data=all_theta_for_gene,
        palette=['#00A087FF', '#4DBBD5FF', '#E64B35FF', '#a65628'],
    )
    plt.title('Boxplot of Theta by Group')
    plt.xlabel('')
    plt.ylabel('Theta')
    plt.grid(False)
    plt.savefig(os.path.join(out_dir, 'Boxplot_of_Theta_by_Group.pdf'),
                bbox_inches='tight')
    plt.close()


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------

def main(mutations_file, reads_file, group_info_file, convs_file, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    mutations_df, reads_df, group_meta, conversions = load_data(
        mutations_file, reads_file, group_info_file, convs_file)

    df_sorted        = preprocess_data(mutations_df, reads_df)
    all_conversions  = sorted(ngs.utils.flatten_iter(conversions))
    df_combined_list = aggregate_conversions(df_sorted, conversions, all_conversions)
    df_sorted, df_combined = merge_group_info(df_sorted, df_combined_list, group_meta)

    # Compute background conversion rate once; reuse for both CSV output and optimisation
    group_pe_values = calculate_mean_conversion_rate(df_sorted, group_meta, CONVERSION_COLUMNS)
    group_pe_values.to_csv(os.path.join(out_dir, 'group_pe_values.csv'), index=False)

    filtered_df          = df_combined[(df_combined['conversion'] != 0) & (df_combined['base'] != 0)]
    df_combined_group_data = sample_data(filtered_df, subset_size=100_000)

    results_df = estimate_pq_values(df_combined_group_data, group_pe_values)

    genedata = df_combined[['base', 'conversion', 'gene_id', 'group']]
    all_theta_for_gene = calculate_theta_for_gene(results_df, genedata)
    all_theta_for_gene.to_csv(os.path.join(out_dir, 'all_theta_for_gene.csv'))

    plot_boxplot_of_theta(all_theta_for_gene, out_dir)
    logging.info("Finished.")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process conversion rate data')
    parser.add_argument('--mutations_file',  help='Path to mutations CSV file')
    parser.add_argument('--reads_file',      help='Path to reads CSV file')
    parser.add_argument('--group_info_file', help='Path to group information CSV file')
    parser.add_argument('--convs_file',      help='Path to conversions pickle file')
    parser.add_argument('--out_dir',         help='Path to output directory')
    args = parser.parse_args()

    main(args.mutations_file, args.reads_file,
         args.group_info_file, args.convs_file, args.out_dir)
