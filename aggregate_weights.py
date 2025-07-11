import argparse
import pandas as pd
import os
from motif_utils.aggregation_utils import (
    create_combined_dataframe, normalize_by_cpgs,
    aggregate_weight_by_mean, aggregate_weights_by_mean_max_frequency,
    aggregate_weights_by_borda_count)


def main():
    """
    Parse CLI args to normalize, group, and aggregate GRNBoost2 outputs,
    then save top-ranked genes with their aggregated weights.
    """
    p = argparse.ArgumentParser(description="Aggregate GRNBoost2 weights")
    p.add_argument('--nruns',    default='all',   help="How many runs to load")
    p.add_argument('--normalize',action='store_true',         help="Normalize by CpG")
    p.add_argument('--group_by', action='store_true',         help="Sum by gene at end")
    p.add_argument('--method',   choices=['mean','borda','freq'], default='mean')
    p.add_argument('--alpha',    type=float, default=1.0,     help="alpha for freq weighting")
    p.add_argument('--beta',     type=float, default=0.5,     help="beta for freq weighting")
    args = p.parse_args()

    grnboost_output_path = 'grnboost_output'

    # Ensure output directory exists
    os.makedirs('aggregation_output', exist_ok=True)

    # Determine number of runs: infer if 'all'
    if args.nruns == 'all':
        files = [f for f in os.listdir(grnboost_output_path) if f.endswith('.tsv')]
        nruns = len(files)
    else:
        nruns = args.nruns

    # Load raw scores matrix
    df, weight_cols = create_combined_dataframe(grnboost_output_path, nruns)
    weight_cols = [f'weight_{i}' for i in range(1, nruns+1)]

    # Optional normalization by CpG sums
    if args.normalize:
        df = normalize_by_cpgs(df, weight_cols)

    # Select aggregation method
    if args.method == 'mean':
        df = aggregate_weight_by_mean(df, weight_cols)
        weight_col = 'mean_weight'
    elif args.method == 'borda':
        df = aggregate_weights_by_borda_count(df, weight_cols)
        weight_col = 'borda_weight'
    else:  # freq
        df = aggregate_weights_by_mean_max_frequency(df, weight_cols,
                                                         alpha=args.alpha,
                                                         beta=args.beta)
        weight_col = 'frequency_weighted'

    # Prepare DataFrame for output
    if isinstance(df, pd.Series):
        df.name = weight_col
        df = df.reset_index().rename(columns={'index': 'gene_cpg'})
    else:
        df = df.reset_index().rename(columns={'index': 'gene_cpg'})

    # Optional grouping by gene keeps highest weight per gene
    if args.group_by:
        df = df.groupby('genes', as_index=False)[weight_col].sum()

    df = df.sort_values(weight_col, ascending=False).drop_duplicates(subset='genes', keep='first')

    # Save final ranked list
    df[['genes', weight_col]].to_csv('aggregation_output/aggregated_weights.tsv', sep='\t', index=False)
    print("Saved results to aggregation_output")

if __name__ == "__main__":
    main()
