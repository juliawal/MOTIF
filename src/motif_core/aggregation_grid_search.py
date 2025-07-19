import pandas as pd
import os
import sys
import json
from src.motif_utils.aggregation_utils import (
    create_combined_dataframe, normalize_by_cpgs,
    create_grid_search_params, 
    aggregate_weight_by_mean,
    aggregate_weights_by_mean_max_frequency,
    aggregate_weights_by_borda_count)


def main(param_grid_path=None):
    """
    Execute grid search over aggregation settings and output per-parameter weight files.
    """
    grnboost_output_path = 'results/grn_inference'
    nruns = 'all'
    nrows = None

    # Load param_grid from file if provided or use default grid
    if param_grid_path and os.path.exists(param_grid_path):
        with open(param_grid_path, 'r') as f:
            param_grid = json.load(f)
        print(f"Loaded parameter grid from: {param_grid_path}")
    else:
        # Default grid
        param_grid = {
            'should_normalize': [True],
            'should_group_by_gene': [False],
            'method': ['mean'],
            'alpha': [10],
            'beta': [2]
        }

    # Ensure output directory exists
    os.makedirs('results/aggregation_grid_search', exist_ok=True)

    # Determine number of runs: infer if 'all'
    if nruns == 'all':
        files = [f for f in os.listdir(grnboost_output_path) if f.endswith('.tsv')]
        nruns = len(files)

    # Load base data matrix
    df, weight_cols = create_combined_dataframe(grnboost_output_path, nruns, nrows)

    # Generate all parameter combinations
    params_list = create_grid_search_params(param_grid)

    for idx, params in enumerate(params_list):
        normalize, group_by, method, alpha, beta = params
        # Save parameter file for this run
        pd.DataFrame([{        
            'should_normalize': normalize,
            'should_group_by_gene': group_by,
            'method': method, 'alpha': alpha, 'beta': beta
        }]).to_csv(f'results/aggregation_grid_search/params_{idx+1}.tsv', sep='\t', index=False)

        # Copy dataframe 
        df_grid = df.copy()
        
        # Reset index to bring 'genes' and 'cpgs' into columns
        df_grid = df_grid.reset_index()

        # Apply normalization if requested
        if normalize: 
            df_grid = normalize_by_cpgs(df_grid, weight_cols)
        # Compute weights based on method
        if method == 'mean':
            df_grid = aggregate_weight_by_mean(df_grid, weight_cols)
        elif method == 'borda':
            df_grid = aggregate_weights_by_borda_count(df_grid, weight_cols)
        else:
            df_grid = aggregate_weights_by_mean_max_frequency(df_grid, weight_cols, alpha, beta)
        # Optionally sum by gene
        if group_by:
            df_grid = df_grid.groupby('genes', as_index=False).sum()
            df_grid = df_grid.drop(columns='cpgs')
        # Save aggregated weights for this grid run
        df_grid.to_csv(f'results/aggregation_grid_search/aggregated_weights_{idx+1}.tsv',
        sep='\t', index=False)


if __name__ == "__main__":
    # For command-line use (optional)
    param_file = sys.argv[1] if len(sys.argv) > 1 else None
    main(param_file)