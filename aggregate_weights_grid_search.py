import pandas as pd
import os
from motif_utils.aggregation_utils import (
    create_combined_dataframe, normalize_by_cpgs,
    create_grid_search_params, 
    aggregate_weight_by_mean,
    aggregate_weights_by_mean_max_frequency,
    aggregate_weights_by_borda_count)


def main():
    """
    Execute grid search over aggregation settings and output per-parameter weight files.
    """
    grnboost_output_path = 'grnboost_output'
    nruns = 3
    nrows = None

    # Ensure output directory exists
    os.makedirs('aggregation_output', exist_ok=True)

    # Define parameter grid for search
    param_grid = {
    'should_normalize': [False],
    'should_group_by_gene': [False],
    'method': ['mean','freq'],
    'alpha': [10],
    'beta': [2]
    }

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
        }]).to_csv(f'aggregation_output/params_{idx+1}.tsv', sep='\t', index=False)
        df_run = df.copy()
        # Apply normalization if requested
        if normalize: 
            df_run = normalize_by_cpgs(df_run, weight_cols)
        # Compute weights based on method
        if method == 'mean':
            df_out = aggregate_weight_by_mean(df_run, weight_cols)
        elif method == 'borda':
            df_out = aggregate_weights_by_borda_count(df_run, weight_cols)
        else:
            df_out = aggregate_weights_by_mean_max_frequency(df_run, weight_cols, alpha, beta)
        # Optionally sum by gene
        if group_by:
            df_out = df_out.groupby('genes', as_index=False).sum()
        # Save aggregated weights for this grid run
        df_out.reset_index()[['genes', 'cpgs', df_out.columns[-1]]].to_csv(
        f'aggregation_output/aggregated_weights_{idx+1}.tsv',
        sep='\t', index=False)


if __name__ == "__main__":
    main()
