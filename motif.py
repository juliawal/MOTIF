import pandas as pd
import sys
from arboreto_added.utils import load_tf_names
from arboreto_added.algo import grnboost2 as grnb2
import zmq
import argparse
import os

def pepare_data(cpgs_file , genes_file, filter_cpgs = False, aggregate_cpgs = False, gene_list = None):

    #  Load data
    cpgs = pd.read_csv(cpgs_file, sep='\t', index_col=[0]).T
    genes = pd.read_csv(genes_file, sep='\t', index_col=[0]).T
    cpgs = cpgs.dropna(axis=1, how='all')
    genes = genes.dropna(axis=1, how='all')

    # Load info for coding/non-coding genes
    path_gene_types = 'utils/gene_types.txt'
    gene_types = pd.read_csv(path_gene_types, index_col=[0], sep='\t')

    # Get list of coding genes without version info
    gene_types['Transcript type'] = gene_types.index
    _coding = gene_types.loc[gene_types['Transcript type'] == ("protein_coding")]
    __coding = gene_types.loc[gene_types['Transcript type'] == ("protein_coding_CDS_not_defined")]
    coding = pd.concat([_coding, __coding])
    coding_list = coding['Gene stable ID'].tolist()

    # Convert gene names to names without version info
    rounded_genes = [x[:15] for x in list(genes.columns)]
    genes = genes.set_axis(rounded_genes, axis=1)

    # Get only coding genes
    intersecting_columns = list(set(genes.columns) & set(coding_list))
    genes = genes[intersecting_columns]

    # Only use 5000 genes with the highest variance in the gene expression matrix
    var = genes.var(axis=0, skipna=True)
    sort = var.sort_values(ascending=False)
    sort5000 = sort[:5000]
    list5000 = list(sort5000.index)
    genes = genes.loc[:, genes.columns.intersection(list5000)]

    if filter_cpgs == True:
        cpgs = filter_cpgs(cpgs, gene_list)

    if aggregate_cpgs == True:
        cpgs = aggregate_cpgs(cpgs)

    # Set genes and cpgs together to the expression matrix
    cpgs.reset_index(drop=True, inplace=True)
    genes.reset_index(drop=True, inplace=True)
    ex_matrix = pd.concat([genes, cpgs],axis=1)
    print(ex_matrix.head())

    # Get transcription factor (genes) names and target (cpgs) names
    _gene_names = pd.DataFrame(genes.columns)
    _gene_names.index.name = None
    _cpg_names = pd.DataFrame(cpgs.columns)
    _cpg_names.index.name = None
    _gene_names.to_csv('gene_names.tsv', sep='\t', index = False, header = False)
    _cpg_names.to_csv('cpg_names.tsv', sep='\t', index = False, header = False)

    # tf_names and target_names is read using a utility function included in Arboreto
    gene_names = load_tf_names('gene_names.tsv')
    cpg_names = load_tf_names('cpg_names.tsv')

    return ex_matrix, gene_names, cpg_names

def run_grnboost2(ex_matrix, gene_names, cpg_names):

    context = zmq.Context()
    context.term()

    # Compute GRN
    network = grnb2(expression_data=ex_matrix,tf_names = gene_names,target_genes = cpg_names)

    out_file = 'output/grn_output.tsv'

    network.to_csv(out_file, sep='\t', index=False, header=False)

# Filter cpgs
def filter_cpgs(cpgs, gene_list):
    # TODO
    return cpgs

# Aggregate cpgs
def aggregate_cpgs(cpgs):
    # TODO
    return cpgs

def main():
    sys.path.insert(0,'/arboreto_added')

    # Argument parser for command-line execution
    parser = argparse.ArgumentParser(description="Run GRNBoost2 analysis with CpG and gene data to get the Methylation Factor Identifiers.")
    parser.add_argument('--cpgs_file', type=str, required=True, help='Path to the CpGs file (TSV format)')
    parser.add_argument('--genes_file', type=str, required=True, help='Path to the genes file (TSV format)')
    parser.add_argument('--filter_cpgs', action='store_true', help='Apply filtering to CpGs')
    parser.add_argument('--filtering_genes_file', type=str, help='Path to the file containing genes to filter CpGs (TSV format)')
    parser.add_argument('--aggregate_cpgs', action='store_true', help='Apply aggregation to CpGs')
    parser.add_argument('--output_dir', type=str, default='output', help='Directory to save the output')

    args = parser.parse_args()

    # Ensure the output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    try:
        # Prepare data
        print("Preparing data...")
        ex_matrix, gene_names, cpg_names = pepare_data(
            cpgs_file=args.cpgs_file,
            genes_file=args.genes_file,
            filter_cpgs=args.filter_cpgs,
            aggregate_cpgs=args.aggregate_cpgs
        )

        # Run GRNBoost2
        print("Running GRNBoost2...")
        run_grnboost2(ex_matrix, gene_names, cpg_names)
        print(f"GRNBoost2 analysis completed. Results saved in '{args.output_dir}'.")

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()