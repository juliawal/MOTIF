import os
import glob
import pandas as pd
import numpy as np
from src.motif_utils.evaluation_utils import (
    load_hippie_network, get_gene_info, map_symbols_to_entrez,
    calculate_neighbors_and_degrees, get_same_degree_genes,
    calculate_p_values_pr, calculate_p_values_hc)

def evaluate_single_run(params_path, weights_path, method='pagerank'):
    """
    Evaluate one parameter setting by computing the mean p-value for selected TFs
    using the specified proximity method on the HIPPIE network,
    and record the intersecting Entrez IDs.
    """
    # Load parameter setting
    params = pd.read_csv(params_path, sep='\t').iloc[0].to_dict()

    # Load aggregated weights and index by gene
    aggregated_weights = pd.read_csv(weights_path, sep='\t')
    aggregated_weights = aggregated_weights.set_index('genes')
    ens_ids = list(dict.fromkeys(aggregated_weights.index))

    # Load HIPPIE network and gene info
    hippie_graph, hippie_df = load_hippie_network('resources/hippie.txt')
    hippie_genes = set(hippie_df['gene1']).union(hippie_df['gene2'])
    hippie_gene_info = get_gene_info(hippie_genes)
    entrez2symbol = map_symbols_to_entrez(hippie_gene_info)
    symbol2entrez = {sym: eid for eid, sym in entrez2symbol.items()}

    # Map ENSEMBL to symbols
    tfs_info = get_gene_info(ens_ids)
    ens2symbol = {entry['query']: entry.get('symbol') for entry in tfs_info if 'symbol' in entry}

    # Filter and map
    filtered = [eid for eid in ens_ids if eid in ens2symbol and ens2symbol[eid] in symbol2entrez]
    top_tfs = filtered[:30]
    symbols = [ens2symbol[eid] for eid in top_tfs]
    intersec_symbols = [sym for sym in symbols if sym in symbol2entrez]
    intersec_symbol2entrez = {sym: symbol2entrez[sym] for sym in intersec_symbols}

    if not intersec_symbol2entrez:
        print(f"Skipping {params_path} – no intersecting TFs")
        return None

    # Define DNMTs and TETs (Entrez IDs)
    dmts_symbol2entrez = {
        'DNMT1': 1786, 'DNMT3A': 1788, 'DNMT3B': 1789,
        'TET1': 80312, 'TET2': 54790, 'TET3': 200424
    }

    neighbors, tf_degrees = calculate_neighbors_and_degrees(
        hippie_graph, list(intersec_symbol2entrez.values())
    )
    all_degrees = dict(hippie_graph.degree())
    same_degree_map = get_same_degree_genes(tf_degrees, all_degrees)

    # Compute p-values
    if method == 'pagerank':
        p_vals = calculate_p_values_pr(
            hippie_graph, tf_degrees, same_degree_map,
            dmts_symbol2entrez, intersec_symbol2entrez
        )
    elif method == 'harmonic':
        p_vals = calculate_p_values_hc(
            hippie_graph, tf_degrees, same_degree_map,
            dmts_symbol2entrez, intersec_symbol2entrez
        )
    else:
        raise ValueError(f"Unknown method: {method}")

    # Attach metrics to parameters
    params['p_value'] = np.mean(p_vals)
    # Join intersecting Entrez IDs into comma-separated string
    entrez_list = sorted(intersec_symbol2entrez.values())
    params['intersec_entrez'] = ",".join(map(str, entrez_list))
    return params


def run_grid_evaluation(method='pagerank'):
    """
    Iterate over all grid search parameter files, evaluate each, and save combined results with intersecting Entrez IDs.
    """
    param_files = sorted(glob.glob('results/aggregation_grid_search/params_*.tsv'))
    results = []

    for param_path in param_files:
        suffix = os.path.basename(param_path).split('_')[-1].replace('.tsv', '')
        weights_path = f'results/aggregation_grid_search/aggregated_weights_{suffix}.tsv'
        if not os.path.exists(weights_path):
            print(f'Missing weights file for {param_path}')
            continue

        result = evaluate_single_run(param_path, weights_path, method)
        if result:
            results.append(result)

    # Save all results
    df_results = pd.DataFrame(results)
    out_path = 'results/evaluation_grid_search/grid_search_results.tsv'
    df_results.to_csv(out_path, sep='\t', index=False)
    print(f'Saved results to {out_path}')


def main():
    # Ensure output directory exists
    os.makedirs('results/evaluation_grid_search', exist_ok=True)

    # Run grid evaluation
    run_grid_evaluation(method='pagerank')  # Or 'harmonic'


if __name__ == "__main__":
    main()