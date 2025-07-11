#python motif.py --cpgs_file data/sample_cpgs.tsv --genes_file data/sample_genes.tsv --run_number 1

import warnings
import sys
import argparse
import os
import pandas as pd
import zmq
from pybiomart import Dataset
from arboreto_added.algo import grnboost2 as grnb2

# Suppress future and specific warnings for cleaner output
warnings.filterwarnings("ignore", category=FutureWarning)


def prepare_data(cpgs_file, genes_file,
                 should_filter_cpgs, filtering_genes_file, filtering_keep_promoters, filtering_keep_bodies,
                 should_aggregate_cpgs, aggregation_genes_file, aggregation_keep_promoters, aggregation_keep_bodies,
                 only_coding_genes, promoter_length, promoter_overlap):
    """
    Load and preprocess CpG and gene expression data.

    - Reads full CpG and gene expression matrices from TSV files.
    - Optionally filters CpGs to those in specified gene regions.
    - Optionally aggregates CpGs by gene region (mean per gene).
    - Selects only coding genes (if requested) and top-5000 variable genes.
    - Returns a combined expression matrix (samples × [genes + CpGs]) plus lists of gene and CpG names.
    """
    # Load data (all rows/columns, using first column as sample index)
    cpgs = pd.read_csv(cpgs_file, sep='\t', index_col=0).transpose()
    genes = pd.read_csv(genes_file, sep='\t', index_col=0).transpose()
    
    # Remove any axis labels for cleanliness
    cpgs.index.name = None
    cpgs.columns.name = None
    genes.index.name = None
    genes.columns.name = None

    # Drop entirely-empty columns (all-NA)
    cpgs = cpgs.dropna(axis=1, how='all')
    genes = genes.dropna(axis=1, how='all')  # remove genes with no expression

    # Load list of protein-coding gene IDs (to filter, if requested)
    # The gene_types.txt file should map transcript types to stable IDs.
    gene_types = pd.read_csv('motif_utils/gene_types.txt', sep='\t')
    coding_list = gene_types[
        gene_types['Transcript type'].isin(['protein_coding', 'protein_coding_CDS_not_defined'])
    ]['Gene stable ID'].tolist()
    # Remove version suffix (after '.') from gene IDs for matching
    trimmed_genes = [g.split('.')[0] for g in genes.columns]
    genes.columns = trimmed_genes

    # Filter to coding genes if requested
    if only_coding_genes:
        mask = genes.columns.isin(coding_list)
        genes = genes.loc[:, mask]

    # Keep only the top 5000 most-variable genes (by variance) to reduce dimensionality
    gene_variances = genes.var(axis=0, skipna=True)
    top_genes = gene_variances.nlargest(5000).index
    genes = genes.loc[:, top_genes]

    # Filter CpGs to genes' promoters/bodies if requested
    if should_filter_cpgs:
        cpgs, _ = filter_cpgs(
            cpgs, filtering_genes_file,
            filtering_keep_promoters, filtering_keep_bodies,
            promoter_length, promoter_overlap
        )

    # Aggregate CpGs by gene region if requested
    if should_aggregate_cpgs:
        cpgs = aggregate_cpgs(
            cpgs, aggregation_genes_file,
            aggregation_keep_promoters, aggregation_keep_bodies,
            promoter_length, promoter_overlap
        )

    # Combine gene expression and CpG data into one matrix
    genes.reset_index(drop=True, inplace=True)
    cpgs.reset_index(drop=True, inplace=True)
    ex_matrix = pd.concat([genes, cpgs], axis=1)

    # Get regulator (gene) and target (CpG) names for GRNBoost2
    gene_names = list(genes.columns)
    cpg_names = list(cpgs.columns)
    return ex_matrix, gene_names, cpg_names


def to_bed_format(cpgs):
    """
    Convert CpG DataFrame into BED-format (chrom, start, end, name, strand).
    Merges CpGs with a probe manifest that maps CpG IDs to genomic coordinates.
    """
    # Load a probe manifest mapping CpG probe IDs to genomic coordinates
    probe_manifest = pd.read_csv(
        "motif_utils/HM450.hg38.manifest.gencode.v36.probeMap",
        sep="\t", header=1,
        names=["probe_ID", "gene", "chrom", "chromStart", "chromEnd", "strand"]
    )
    # Reset index to turn CpG IDs into a column for merging
    cpgs_bed = cpgs.T.reset_index().rename(columns={'index': 'probe_ID'})
    # Merge on probe_ID to attach coordinates
    merged = pd.merge(cpgs_bed, probe_manifest, how="inner", left_on="probe_ID", right_on="probe_ID")
    # Select and rename columns for BED (0-based start)
    bed = merged[["chrom", "chromStart", "chromEnd", "probe_ID", "strand"]].copy()
    bed.rename(columns={"probe_ID": "name"}, inplace=True)
    return bed


def filter_cpgs(cpgs, filtering_genes_file, keep_promoters, keep_bodies, promoter_length, promoter_overlap):
    """
    Filter CpG columns to those falling within promoters or bodies of specified genes.
    Returns a filtered CpG DataFrame and a mapping of CpG->gene.
    """
    # Read gene list (one Ensembl gene ID per line)
    with open(filtering_genes_file, 'r') as file:
        genes_of_interest = [line.strip() for line in file]

    # Get genomic regions for each gene in genes_of_interest
    dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
    regions = get_gene_coordinates(
        genes_of_interest, dataset,
        keep_promoters, keep_bodies,
        promoter_length, promoter_overlap
    )

    # Convert CpG matrix to BED-format for overlap checking
    bed_cpgs = to_bed_format(cpgs)

    matches = []
    # Find overlaps: for each CpG, check if it lies in any gene region
    for _, cpg_row in bed_cpgs.iterrows():
        for _, gene_row in regions.iterrows():
            if (cpg_row['chrom'] == gene_row['chrom'] and
                cpg_row['chromStart'] >= gene_row['chromStart'] and
                cpg_row['chromEnd'] <= gene_row['chromEnd']):
                matches.append({
                    'probe_ID': cpg_row['name'],
                    'ensembl_gene_id': gene_row['ensembl_gene_id']
                })
    if not matches:
        print("Error: No CpGs found in specified gene regions.")
        sys.exit(1)

    cpgs_in_genes = pd.DataFrame(matches)
    selected_cpgs = cpgs_in_genes['probe_ID'].tolist()
    filtered_cpgs = cpgs.loc[:, cpgs.columns.intersection(selected_cpgs)]
    return filtered_cpgs, cpgs_in_genes


def aggregate_cpgs(cpgs, genes_file, keep_promoters, keep_bodies, promoter_length, promoter_overlap):
    """
    Aggregate CpG values by gene region. For genes with multiple CpGs in their region,
    replaces those CpG columns by their mean.
    """
    filtered_cpgs, mapping = filter_cpgs(
        cpgs, genes_file,
        keep_promoters, keep_bodies,
        promoter_length, promoter_overlap
    )
    # Group CpGs by gene and compute mean across each group (row-wise)
    grouped = mapping.groupby('ensembl_gene_id')['probe_ID'].apply(list).to_dict()
    cols_to_drop = []
    for gene_id, probes in grouped.items():
        if len(probes) > 1:
            agg_col = f"cpgs_agg_for_{gene_id}"
            cpgs[agg_col] = cpgs[probes].mean(axis=1)
            cols_to_drop.extend(probes)
    cpgs.drop(columns=cols_to_drop, inplace=True, errors='ignore')
    return cpgs


def get_gene_coordinates(genes_list, dataset, keep_promoters, keep_bodies, promoter_length, promoter_overlap):
    """
    Query Ensembl via pybiomart for gene start/end and strand, then define promoter/body regions.
    Returns a DataFrame with columns: ['ensembl_gene_id','chrom','chromStart','chromEnd'].
    """
    # Query the dataset for relevant fields
    attributes = ['ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'strand']
    gene_info = dataset.query(attributes=attributes)
    # Filter to requested genes
    gene_info = gene_info[gene_info['Gene stable ID'].isin(genes_list)]
    regions = []
    for _, row in gene_info.iterrows():
        gid = row['Gene stable ID']
        chrom = f"chr{row['Chromosome/scaffold name']}"
        start = int(row['Gene start (bp)'])
        end = int(row['Gene end (bp)'])
        strand = int(row['Strand'])

        if keep_promoters:
            # Calculate promoter region based on strand
            if strand == 1:  # plus strand
                prom_start = max(0, start - promoter_length)
                prom_end = start + promoter_overlap
            else:  # minus strand
                prom_start = end - promoter_overlap
                prom_end = end + promoter_length
        else:
            prom_start = None
            prom_end = None

        if keep_promoters and keep_bodies:
            regions.append({'ensembl_gene_id': gid,
                            'chrom': chrom,
                            'chromStart': prom_start,
                            'chromEnd': end})
        else:
            if keep_promoters:
                regions.append({'ensembl_gene_id': gid,
                                'chrom': chrom,
                                'chromStart': prom_start,
                                'chromEnd': prom_end})
            if keep_bodies:
                regions.append({'ensembl_gene_id': gid,
                                'chrom': chrom,
                                'chromStart': start,
                                'chromEnd': end})
    return pd.DataFrame(regions)


def run_grnboost2(ex_matrix, gene_names, cpg_names, grnboost_output_path):
    """
    Run the adapted GRNBoost2 algorithm to infer a regulatory network.
    Saves the output (gene-CpG importance) to a TSV file.
    """
    # Clean up any existing ZMQ contexts (as recommended by arboreto)
    context = zmq.Context()
    context.term()

    # Run GRNBoost2: gene_names = TF/regulators, cpg_names = targets
    network = grnb2(expression_data=ex_matrix, tf_names=gene_names, target_genes=cpg_names)

    # Save the inferred network to file
    network.to_csv(grnboost_output_path, sep='\t', index=False)

def main():
    """
    Parse command-line arguments, prepare data, and run GRNBoost2.
    """
    parser = argparse.ArgumentParser(
        description="Infer GRN from gene expression and CpG methylation data using GRNBoost2."
    )
    parser.add_argument('--cpgs_file', type=str, required=True, help='Path to CpG TSV file')
    parser.add_argument('--genes_file', type=str, required=True, help='Path to gene TSV file')
    parser.add_argument('--should_filter_cpgs', action='store_true', help='Filter CpGs by gene regions')
    parser.add_argument('--filtering_genes_file', type=str, help='Gene list for filtering CpGs')
    parser.add_argument('--filtering_keep_promoters', action='store_true', help='Include promoters in filter')
    parser.add_argument('--filtering_keep_bodies', action='store_true', help='Include gene bodies in filter')
    parser.add_argument('--should_aggregate_cpgs', action='store_true', help='Aggregate CpGs by gene regions')
    parser.add_argument('--aggregation_genes_file', type=str, help='Gene list for aggregating CpGs')
    parser.add_argument('--aggregation_keep_promoters', action='store_true', help='Include promoters in aggregation')
    parser.add_argument('--aggregation_keep_bodies', action='store_true', help='Include bodies in aggregation')
    parser.add_argument('--only_coding_genes', action='store_true', help='Retain only protein-coding genes')
    parser.add_argument('--promoter_length', type=int, default=1000, help='Upstream promoter length (bp)')
    parser.add_argument('--promoter_overlap', type=int, default=200,  help='Overlap region in promoters (bp)')
    parser.add_argument('--run_number', type=int, default=1, help='Unique number to identify this GRNBoost run during output saving.')
    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs('grnboost_output', exist_ok=True)
    grnboost_output_path = f'grnboost_output/grn_output_{args.run_number}.tsv'

    # Prepare data matrix and gene/CpG lists
    ex_matrix, gene_names, cpg_names = prepare_data(
        args.cpgs_file, args.genes_file,
        args.should_filter_cpgs, args.filtering_genes_file,
        args.filtering_keep_promoters, args.filtering_keep_bodies,
        args.should_aggregate_cpgs, args.aggregation_genes_file,
        args.aggregation_keep_promoters, args.aggregation_keep_bodies,
        args.only_coding_genes, args.promoter_length, args.promoter_overlap
    )

    # Run GRNBoost2 and save output
    print("Running GRNBoost2...")
    run_grnboost2(ex_matrix, gene_names, cpg_names, grnboost_output_path)
    print(f"Network inference complete. Result saved to grnboost_output.")

if __name__ == '__main__':
    main()