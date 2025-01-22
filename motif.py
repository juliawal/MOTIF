# Suppress 'Dask dataframe query planning is disabled because dask-expr is not installed'-warning and for this all FutureWarnings globally
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import pandas as pd
import sys
from arboreto_added.utils import load_tf_names
from arboreto_added.algo import grnboost2 as grnb2
import zmq
import argparse
import os
from pybiomart import Dataset

def prepare_data(cpgs_file , genes_file, should_filter_cpgs, filtering_genes_file, filtering_keep_promoters, filtering_keep_bodies, should_aggregate_cpgs, aggregation_genes_file, aggregation_keep_promoters, aggregation_keep_bodies, only_coding_genes, promoter_length, promoter_overlap):

    # Load data
    cpgs = pd.read_csv(cpgs_file, sep='\t', index_col=[0]).T
    genes = pd.read_csv(genes_file, sep='\t', index_col=[0]).T
    cpgs = cpgs.dropna(axis=1, how='all')
    genes = genes.dropna(axis=1, how='all')

    # Get list of coding genes without version info
    path_gene_types = 'utils/gene_types.txt'
    gene_types = pd.read_csv(path_gene_types, index_col=[0], sep='\t')
    gene_types['Transcript type'] = gene_types.index
    _coding = gene_types.loc[gene_types['Transcript type'] == ("protein_coding")]
    __coding = gene_types.loc[gene_types['Transcript type'] == ("protein_coding_CDS_not_defined")]
    coding = pd.concat([_coding, __coding])
    coding_list = coding['Gene stable ID'].tolist()

    # Convert gene names to names without version info
    rounded_genes = [x[:15] for x in list(genes.columns)]
    genes = genes.set_axis(rounded_genes, axis=1)

    # Get only coding genes if coding_genes == True
    if only_coding_genes == True:
        intersecting_columns = list(set(genes.columns) & set(coding_list))
        genes = genes[intersecting_columns]

    # Only use 5000 genes with the highest variance in the gene expression matrix
    var = genes.var(axis=0, skipna=True)
    sort = var.sort_values(ascending=False)
    sort5000 = sort[:5000]
    list5000 = list(sort5000.index)
    genes = genes.loc[:, genes.columns.intersection(list5000)]

    # Only retain CpGs that are in the promoters or bodies of genes listed in filtering_genes_file
    if should_filter_cpgs:
        cpgs, cpgs_in_genes = filter_cpgs(cpgs, filtering_genes_file, filtering_keep_promoters, filtering_keep_bodies, promoter_length, promoter_overlap)

    # Aggregate CpGs that are in the promoters or bodies of genes listed in aggregation_genes_file
    if should_aggregate_cpgs:
        cpgs = aggregate_cpgs(cpgs, filtering_genes_file, aggregation_keep_promoters, aggregation_keep_bodies, promoter_length, promoter_overlap)

    # Set genes and cpgs together to the expression matrix
    cpgs.reset_index(drop=True, inplace=True)
    genes.reset_index(drop=True, inplace=True)
    ex_matrix = pd.concat([genes, cpgs],axis=1)


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

def run_grnboost2(ex_matrix, gene_names, cpg_names, nruns):
    for i in range(1,nruns+1):
        context = zmq.Context()
        context.term()

        # Compute GRN
        output = grnb2(expression_data = ex_matrix, tf_names = gene_names, target_genes = cpg_names)
        out_file = 'output/grn_output_' + str(i) + '.tsv'

        output.to_csv(out_file, sep='\t', index=False, header=False)

def to_bed_format(cpgs):
    # Load the probe manifest file
    probe_manifest = pd.read_csv("utils/HM450.hg38.manifest.gencode.v36.probeMap", sep="\t", header=1)
    probe_manifest.columns = ["probe_ID", "gene", "chrom", "chromStart", "chromEnd", "strand"]

    # Load methylation data in the right forme
    cpgs = cpgs.T
    cpgs = cpgs.reset_index()
    cpgs = cpgs.rename(columns={'index': 'Composite Element REF'})

    # Merge the methylation data with the probe manifest on the Probe_ID (CpG ID)
    merged_data = pd.merge(cpgs, probe_manifest, left_on="Composite Element REF", right_on="probe_ID")
    
    # Select relevant columns for BED format
    bed_cpgs = merged_data[["chrom", "chromStart", "chromEnd", "Composite Element REF", "strand"]]
    bed_cpgs.rename(columns={"Composite Element REF": "name"})
    return bed_cpgs

def filter_cpgs(cpgs, filtering_genes_file, keep_promoters, keep_bodies, promoter_length, promoter_overlap):
    # Read the file line by line into a list
    with open(filtering_genes_file, 'r') as file:
        filtering_genes = [line.strip() for line in file]
    
    filtered_cpgs = pd.DataFrame()

    dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
    
    # Get coordinates for genes in filtering_genes
    coordinates_df = get_gene_coordinates(filtering_genes, keep_promoters=keep_promoters, keep_bodies=keep_bodies, dataset = dataset, promoter_length = promoter_length, promoter_overlap = promoter_overlap)
    bed_cpgs = to_bed_format(cpgs)

    result = []

    for _, cpg_row in bed_cpgs.iterrows():
        # Extract relevant CpG region details
        cpg_chrom = cpg_row['chrom']
        cpg_start = cpg_row['chromStart']
        cpg_end = cpg_row['chromEnd']
        cpg_ref = cpg_row['Composite Element REF']
        
        for _, gene_row in coordinates_df.iterrows():
            # Extract gene region details
            gene_chrom = gene_row['chrom']
            gene_start = gene_row['chromStart']
            gene_end = gene_row['chromEnd']
            gene_id = gene_row['ensembl_gene_id']
            
            # Check if the CpG region lies within the gene region
            if cpg_chrom == gene_chrom:
                if cpg_start >= gene_start and cpg_end <= gene_end:
                    result.append({'Composite Element REF': cpg_ref, 'ensembl_gene_id': gene_id})

    # Convert the result to a dataframe
    if not result:
        print("Error: None of the given CpGs are located in the region of the genes specified in filtering_genes_file.")
        sys.exit(1)

    cpgs_in_genes = pd.DataFrame(result)

    # Extract the list of CpGs from cpgs_in_genes
    elements_to_filter = cpgs_in_genes['Composite Element REF'].tolist()
    filtered_cpgs = cpgs[[col for col in elements_to_filter if col in cpgs.columns]]

    return filtered_cpgs, cpgs_in_genes

def aggregate_cpgs(cpgs, genes_list, keep_promoters, keep_bodies, promoter_length, promoter_overlap):
    filtered_cpgs, cpgs_in_genes = filter_cpgs(cpgs, genes_list, keep_promoters, keep_bodies, promoter_length, promoter_overlap)
    
    # To test aggregate_cpgs()
    # cpgs_in_genes.loc[1, "ensembl_gene_id"] = cpgs_in_genes.loc[0, "ensembl_gene_id"]

    # Create dictionary with genes in Ensembl ID as keys and CpGs, that are in the DNA region of the Gene, as values
    gene_dict = cpgs_in_genes.groupby('ensembl_gene_id')['Composite Element REF'].apply(list).to_dict()
    cpg_columns_to_drop = []

    # Iterate over the keys directly
    for gene in gene_dict:
        if len(gene_dict[gene]) > 1:
            to_aggregate = gene_dict[gene]
            # Sum the relevant columns and create a new column in the cpgs DataFrame
            cpgs['cpgs_aggregated_for_' + str(gene)] = cpgs[to_aggregate].mean(axis=1)
            cpg_columns_to_drop.extend(to_aggregate)

    cpgs.drop(columns=cpg_columns_to_drop, inplace=True)
    print(cpgs)

    return cpgs

def save_n_grnboost_runs_to_file(nruns):
    with open('output/nruns.txt', 'w') as file:
        file.write(str(nruns))

def get_gene_coordinates(genes_list, dataset, keep_promoters, keep_bodies, promoter_length, promoter_overlap):
    # Query for basic gene information
    attributes = ['ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'strand']
    gene_info = dataset.query(attributes=attributes) #'link_ensembl_gene_id'
    filtered_gene_info = gene_info[gene_info['Gene stable ID'].isin(genes_list)]
    results = []
    for _, row in filtered_gene_info.iterrows():
            ensembl_id = row['Gene stable ID']
            chrom = row['Chromosome/scaffold name']
            start = row['Gene start (bp)']
            end = row['Gene end (bp)']
            strand = row['Strand']
         
            #promoter_start, promoter_end = start, start  # Initialize promoter region
            if keep_promoters:
                # Calculate promoter region based on strand
                if strand == 1:  # Positive strand
                    promoter_start = max(0, start - promoter_length)
                    promoter_end = start + promoter_overlap
                else:  # Negative strand
                    promoter_start = end - promoter_overlap
                    promoter_end = end + promoter_length

            if keep_promoters and keep_bodies:
                # If both promoters and bodies are to be included, extend from promoter start to gene end
                region_start = promoter_start
                region_end = end
                results.append({
                    'ensembl_gene_id': ensembl_id,
                    'region_type': 'promoter_and_body',
                    'chrom': 'chr' + str(chrom),
                    'chromStart': region_start,
                    'chromEnd': region_end
                })
            else:
                if keep_promoters:
                    results.append({
                        'ensembl_gene_id': ensembl_id,
                        'region_type': 'promoter',
                        'chrom': 'chrom' + str(chrom),
                        'chromStart': promoter_start,
                        'chromEnd': promoter_end
                    })

                if keep_bodies:
                    results.append({
                        'ensembl_gene_id': ensembl_id,
                        'region_type': 'body',
                        'chrom': 'chrom' + str(chrom),
                        'chromStart': start,
                        'chromEnd': end
                    })
    return pd.DataFrame(results)

def main():
    sys.path.insert(0,'/arboreto_added')

    # Argument parser for command-line execution
    parser = argparse.ArgumentParser(description="Run GRNBoost2 analysis with CpG and gene data to get the Methylation Factor Identifiers.")
    parser.add_argument('--cpgs_file', type=str, required=True, help='Path to the CpGs file (TSV format)')
    parser.add_argument('--genes_file', type=str, required=True, help='Path to the genes file (TSV format)')
    parser.add_argument('--should_filter_cpgs', action='store_true', help='Apply filtering to CpGs')
    parser.add_argument('--filtering_genes_file', type=str, help='Path to the file containing genes to filter CpGs (TSV format)')
    parser.add_argument('--filtering_keep_promoters', action='store_true', help='Keep promoter regions during filtering')
    parser.add_argument('--filtering_keep_bodies', action='store_true', help='Keep gene bodies during filtering')
    parser.add_argument('--should_aggregate_cpgs', action='store_true', help='Apply aggregation to CpGs')
    parser.add_argument('--aggregation_genes_file', type=str, help='Path to the file containing genes for aggregation (TSV format)')
    parser.add_argument('--aggregation_keep_promoters', action='store_true', help='Keep promoter regions during aggregation')
    parser.add_argument('--aggregation_keep_bodies', action='store_true', help='Keep gene bodies during aggregation')
    parser.add_argument('--output_dir', type=str, default='output', help='Directory to save the output')
    parser.add_argument('--only_coding_genes', action='store_true', help='Select only coding genes')
    parser.add_argument('--promoter_length', type=int, default=1000, help='Length in base pairs to define the promoter region') 
    parser.add_argument('--promoter_overlap', type=int, default=200, help='Length in base pairs to define the promoter overlap') 
    parser.add_argument('--nruns', type=int, default=10, help='Number of times to run GRNBoost2')

    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    try:
        # Prepare data
        print("Preparing data...")
        ex_matrix, gene_names, cpg_names = prepare_data(
            cpgs_file = args.cpgs_file,
            genes_file = args.genes_file,
            should_filter_cpgs = args.should_filter_cpgs,
            filtering_genes_file = args.filtering_genes_file,
            filtering_keep_promoters = args.filtering_keep_promoters,
            filtering_keep_bodies = args.filtering_keep_bodies,
            should_aggregate_cpgs = args.should_aggregate_cpgs,
            aggregation_genes_file = args.aggregation_genes_file,
            aggregation_keep_promoters = args.aggregation_keep_promoters,
            aggregation_keep_bodies = args.aggregation_keep_bodies,
            only_coding_genes = args.only_coding_genes,
            promoter_length = args.promoter_length,
            promoter_overlap = args.promoter_overlap
        )

        # Run GRNBoost2
        print("Running GRNBoost2...")
        run_grnboost2(ex_matrix, gene_names, cpg_names, nruns = args.nruns)
        print(f"GRNBoost2 analysis completed. Results saved in '{args.output_dir}'.")

        # Save nruns for evaluation
        save_n_grnboost_runs_to_file(args.nruns)

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
