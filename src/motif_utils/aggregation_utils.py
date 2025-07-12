import pandas as pd
import itertools

# Common functions for GRNBoost2 output processing and aggregation

def prepare_grnboost_output(run_index, file_path, nrows=None):
    """
    Reads one GRNBoost2 TSV and returns a DataFrame indexed by (genes, cpgs).
    """
    df = pd.read_csv(
        f"{file_path}/grn_output_{run_index}.tsv",
        sep='\t',
        dtype={'TF': str},
        nrows=nrows
    )
    df = df.rename(columns={'TF': 'genes', 'target': 'cpgs', 'importance': f'weight_{run_index}'})
    return df.set_index(['genes', 'cpgs'])


def create_combined_dataframe(file_path, nruns, nrows=None):
    """
    Loads all runs and concatenates their weights into one DataFrame.
    """
    dfs = [prepare_grnboost_output(i, file_path, nrows) for i in range(1, nruns+1)]
    combined = pd.concat(dfs, axis=1)
    weight_cols = [f'weight_{i}' for i in range(1, nruns+1)]
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
