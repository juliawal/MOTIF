import pandas as pd
import itertools
import os
import re
import numpy as np


# Common functions for GRNBoost2 output processing and aggregation


def extract_seeds(folder_path):
    """
    Extracts all seed numbers from filenames in the folder of the form: grn_with_seed_<seed>.tsv
    """
    seeds = []
    pattern = re.compile(r'grn_with_seed_(\d+)\.tsv')
    for fname in os.listdir(folder_path):
        match = pattern.match(fname)
        if match:
            seeds.append(int(match.group(1)))
    return sorted(seeds)


def prepare_grnboost_output(seed, folder_path, nrows=None):
    """
    Loads the GRNBoost2 output for a given seed.
    """
    file_path = os.path.join(folder_path, f'grn_with_seed_{seed}.tsv')
    df = pd.read_csv(file_path, sep='\t', dtype={'TF': str}, nrows=nrows)

    if 'importance' not in df.columns:
        raise ValueError(f"File {file_path} does not contain an 'importance' column.")

    df = df.rename(columns={'TF': 'genes', 'target': 'cpgs', 'importance': f'weight_{seed}'})

    return df.set_index(['genes', 'cpgs'])


def create_combined_dataframe(folder_path, nruns=None, nrows=None):
    """
    Loads GRNBoost2 outputs from available seed files and concatenates their weights.
    If nruns is set, randomly selects that many seeds.
    """
    all_seeds = extract_seeds(folder_path)

    if nruns is not None and nruns < len(all_seeds):
        selected_seeds = sorted(np.random.choice(all_seeds, nruns, replace=False))
    else:
        selected_seeds = all_seeds

    print(f"Using seeds: {selected_seeds}")

    dfs = [prepare_grnboost_output(seed, folder_path, nrows) for seed in selected_seeds]
    combined = pd.concat(dfs, axis=1)
    weight_cols = [f'weight_{seed}' for seed in selected_seeds]
    combined[weight_cols] = combined[weight_cols].fillna(0)
    return combined, weight_cols


def normalize_by_cpgs(df, weight_cols):
    """Divide each weight by the sum across its CpG site."""
    cpg_sums = df.groupby('cpgs')[weight_cols].sum()
    merged = df.join(cpg_sums, on='cpgs', rsuffix='_sum')
    for col in weight_cols:
        sum_col = col + '_sum'
        mask = merged[sum_col] != 0
        merged.loc[mask, col] = merged.loc[mask, col] / merged.loc[mask, sum_col]
        merged.loc[~mask, col] = 0
    return merged.drop(columns=[c + '_sum' for c in weight_cols]).fillna(0)


def aggregate_weight_by_mean(df, weight_cols):
    """
        Compute the average importance across all runs for each edge and sort descending.
        """
    df['mean_weight'] = df[weight_cols].mean(axis=1)
    return df.sort_values('mean_weight', ascending=False)


def aggregate_weights_by_mean_max_frequency(df, weight_cols, alpha, beta):
    """
    Combine frequency (number of nonzero runs) with run weights: a weighted mix of max and mean
    (beta balances max vs mean), then scale by freq^alpha to emphasize consistently high edges.
    """
    freq = (df[weight_cols] > 0).sum(axis=1)
    mean = df[weight_cols].mean(axis=1)
    mx = df[weight_cols].max(axis=1)
    df['mean_max_frequency_weight'] = (beta * mx + (1 - beta) * mean) * (freq ** alpha)
    return df.sort_values('mean_max_frequency_weight', ascending=False)


def aggregate_weights_by_borda_count(df, weight_cols):
    """
    Rank each run's weights, then apply Borda count to combine these rankings into a global weight ranking.
    """
    import ranky as rky
    for col in weight_cols:
        df[f'{col}_rank'] = df[col].rank(ascending=False)
    rank_cols = [c for c in df.columns if c.endswith('_rank')]
    aggr_ranking = rky.borda(df[rank_cols], axis=1)
    return aggr_ranking.sort_values(ascending=False)


def create_grid_search_params(param_grid):
    """
    Expand grid into list of parameter tuples.
    """
    combos = []
    for normalize, group_by, method in itertools.product(
            param_grid['should_normalize'],
            param_grid['should_group_by_gene'],
            param_grid['method']):
        if method == 'freq':
            for alpha in param_grid['alpha']:
                for beta in param_grid['beta']:
                    combos.append((normalize, group_by, method, alpha, beta))
        else:
            combos.append((normalize, group_by, method, None, None))
    return combos
