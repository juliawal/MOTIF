import logging
import random
import numpy as np
import pandas as pd
import networkx as nx
import mygene

def get_gene_info(genes):
    """
    Queries MyGene.info to retrieve gene symbols for a list of gene IDs.
    """
    logging.disable(logging.CRITICAL)
    mg = mygene.MyGeneInfo()
    info = mg.getgenes(list(genes), fields='symbol')
    logging.disable(logging.NOTSET)
    return info


def map_symbols_to_entrez(gene_info):
    """
    Builds a dict mapping Entrez ID to gene symbol from MyGene.info output.
    """
    return {int(entry['_id']): entry['symbol'] for entry in gene_info if 'symbol' in entry}

def load_hippie_network(file_path):
    """
    Loads the HIPPIE PPI network and returns a NetworkX graph and the DataFrame.
    """
    hippie_df = pd.read_csv(file_path, sep='\t', usecols=[1, 3], names=['gene1', 'gene2'])
    G = nx.Graph()
    G.add_edges_from(zip(hippie_df['gene1'], hippie_df['gene2']))
    return G, hippie_df

def calculate_neighbors_and_degrees(G, genes):
    """
    Returns neighbors list and degree dict for each gene in 'genes'.
    """
    neighbors = {g: list(G.neighbors(g)) for g in genes}
    degrees = {g: len(neigh) for g, neigh in neighbors.items()}
    return neighbors, degrees


def get_same_degree_genes(degrees, graph_degrees):
    """
    Maps each gene to other genes with the same or nearest degree in the network.
    """
    mapping = {}
    unique_degs = set(graph_degrees.values())
    for gene, deg in degrees.items():
        step = 0
        same = []
        while not same:
            possible = {deg + step, deg - step} & unique_degs
            same = [g for g, dg in graph_degrees.items() if dg in possible]
            step += 1
        mapping[gene] = same
    return mapping

def calculate_p_values_pr(graph, tf_degrees, same_deg_map, dmts, intersec, n_iter=30, random_iter=1000):
    """
    Computes p-values using personalized PageRank.
    """
    pvals = []
    pers = {v: 1/len(dmts) for v in dmts.values()}
    for _ in range(n_iter):
        rand_sets = [[random.choice(same_deg_map[tf]) for tf in tf_degrees] for _ in range(random_iter)]
        weights = nx.pagerank(graph, personalization=pers)
        true_vals = [weights[k] for k in intersec.values() if k in weights]
        if not true_vals:
            continue
        m_true = np.mean(true_vals)
        rand_means = []
        for rs in rand_sets:
            vals = [weights[k] for k in rs if k in weights]
            if vals:
                rand_means.append(np.mean(vals))
        if not rand_means:
            continue
        p = (np.sum(np.array(rand_means) >= m_true) + 1) / (len(rand_means) + 1)
        pvals.append(p)
    return pvals


def relative_harmonic_centrality(G, sources, targets):
    """
    For each node s in `sources`, computes HC_rel(s) = sum_{t in targets} 1 / dist(s, t),
    ignoring unreachable pairs. Returns dict: {s: HC_rel(s)}.
    """
    hc = {}
    for s in sources:
        lengths = nx.single_source_shortest_path_length(G, s)
        hc[s] = sum(1/lengths[t] for t in targets if t in lengths and lengths[t] > 0)
    return hc


def calculate_p_values_hc(graph, tf_degrees, same_deg_map, dmts, intersec, n_iter=30, random_iter=1000):
    """
    Computes p-values using relative harmonic centrality.
    """
    sources = list(intersec.values())
    targets = list(dmts.values())
    true_hc = relative_harmonic_centrality(graph, sources, targets)
    m_true = np.mean(list(true_hc.values()))
    pvals = []
    for _ in range(n_iter):
        rand_sets = [[random.choice(same_deg_map[src]) for src in sources] for _ in range(random_iter)]
        rand_means = []
        for rs in rand_sets:
            hc_rand = relative_harmonic_centrality(graph, rs, targets)
            rand_means.append(np.mean(list(hc_rand.values())))
        p = (np.sum(np.array(rand_means) >= m_true) + 1) / (random_iter + 1)
        pvals.append(p)
    return pvals