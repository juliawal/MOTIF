import numpy as np
import pandas as pd
from motif_utils.evaluation_utils import (
    get_gene_info, map_symbols_to_entrez,
    load_hippie_network, calculate_neighbors_and_degrees,
    get_same_degree_genes, calculate_p_values_pr,
    calculate_p_values_hc)


def main():
    """
    Load aggregated TF weights, map to network genes, and compute mean p-values
    using either PageRank or Harmonic Centrality proximity measures.
    """
    # Choose the method for proximity calculation: 'pagerank' or 'harmonic'
    method = 'pagerank'

    # Load aggregated weights and index by gene
    aggregated_weights = pd.read_csv('aggregation_output/aggregated_weights.tsv', sep='\t')
    aggregated_weights = aggregated_weights.set_index('genes')

    # Unique gene IDs in ranked order
    ens_ids = list(dict.fromkeys(aggregated_weights.index))
    
    # Load HIPPIE PPI network and gene info
    hippie_graph, hippie_df = load_hippie_network('motif_utils/hippie.txt')
    hippie_genes = set(hippie_df['gene1']).union(hippie_df['gene2'])
    hippie_gene_info = get_gene_info(hippie_genes)
    entrez2symbol = map_symbols_to_entrez(hippie_gene_info)
    symbol2entrez = {sym: eid for eid, sym in entrez2symbol.items()}
    
    # Map ENSEMBL IDs to gene symbols
    tfs_info = get_gene_info(ens_ids)
    ens2symbol = {entry['query']: entry.get('symbol') for entry in tfs_info if 'symbol' in entry}
    
    # Filter ENS IDs to those with a mapped symbol present in HIPPIE
    filtered = [eid for eid in ens_ids if eid in ens2symbol and ens2symbol[eid] in symbol2entrez]
    ensID_top30TFs = filtered[:30]
    intersec_symbols = [ens2symbol[eid] for eid in ensID_top30TFs if eid in ens2symbol]
    intersec_symbols = [sym for sym in intersec_symbols if sym in symbol2entrez]
    intersec_symbol2entrez = {sym: symbol2entrez[sym] for sym in intersec_symbols}
    
    # Define DNMTs and TETs (Entrez IDs)
    dmts_symbol2entrez = {
        'DNMT1': 1786, 'DNMT3A': 1788, 'DNMT3B': 1789,
        'TET1': 80312, 'TET2': 54790, 'TET3': 200424
    }
    
    if len(intersec_symbol2entrez) == 0:
        raise ValueError("Error: Not enough TFs selected for analysis.")
    
    # Compute neighbors and degrees for selected TFs
    neighbors, tf_degrees = calculate_neighbors_and_degrees(hippie_graph, list(intersec_symbol2entrez.values()))
    all_degrees = dict(hippie_graph.degree())
    same_degree_map = get_same_degree_genes(tf_degrees, all_degrees)
    

    if method == 'pagerank':
        p_values = calculate_p_values_pr(
            hippie_graph, tf_degrees, same_degree_map,
            dmts_symbol2entrez, intersec_symbol2entrez
        )
        print("The mean p-value (PageRank) is", np.mean(p_values))

    elif method == 'harmonic':
        p_values = calculate_p_values_hc(
            hippie_graph, tf_degrees, same_degree_map,
            dmts_symbol2entrez, intersec_symbol2entrez
        )
        print("The mean p-value (Harmonic Centrality) is", np.mean(p_values))


if __name__ == "__main__":
    main()