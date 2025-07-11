# Installation

To install the required environment and dependencies, use the provided setup.sh script.

```bash
bash setup.sh
unzip motif_utils/hippie.txt.zip
```

# Running MOTIF

To run MOTIF, use the following command:
```bash
conda activate motif
python motif.py --cpgs_file data/sample_cpgs.tsv --genes_file data/sample_genes.tsv --only_coding_genes --nruns 20
```

The positional arguments are:
```
--cpgs_file
Path to the file containing the CpGs data (e.g., data/sample_cpgs.tsv).
Type: str | Required: Yes

--genes_file
Path to the file containing the genes data (e.g., data/sample_genes.tsv).
Type: str | Required: Yes

--should_filter_cpgs
Set to True if you want to filter the CpGs based on the provided genes.
Type: bool | Required: No | Default: False

--filtering_genes_file
Path to the file containing a list of specific genes to filter CpGs (e.g., data/filter_genes.tsv).
Type: str | Required: No (only required if --should_filter_cpgs is set)

--filtering_keep_promoters
Set to True if you want to keep CpGs that are in the promoter regions of the specified genes during filtering.
Type: bool | Required: No | Default: False

--filtering_keep_bodies
Set to True if you want to keep CpGs that are in the gene bodies of the specified genes during filtering.
Type: bool | Required: No | Default: False

--should_aggregate_cpgs
Set to True if you want to aggregate the CpGs data for the provided genes.
Type: bool | Required: No | Default: False

--aggregation_genes_file
Path to the file containing a list of genes for aggregation (e.g., data/aggregate_genes.tsv).
Type: str | Required: No (only required if --should_aggregate_cpgs is set)

--aggregation_keep_promoters
Set to True if you want to aggregate CpGs that are in the promoter regions of the specified genes.
Type: bool | Required: No | Default: False

--aggregation_keep_bodies
Set to True if you want to aggregate CpGs that are in the gene bodies of the specified genes.
Type: bool | Required: No | Default: False

--output_dir
Path to the directory where the output files will be saved.
Type: str | Required: No | Default: output

--only_coding_genes
Set to True if you only want to use protein-coding genes for the analysis.
Type: bool | Required: No | Default: False

--promoter_length
Length in base pairs to define the promoter region.
Type: int | Required: No | Default: 1000

--promoter_overlap
Length in base pairs to define the overlap region between the promoter and the gene body.
Type: int | Required: No | Default: 200

--nruns
Number of times to run the GRNBoost2 algorithm.
Type: int | Required: No | Default: 10
```

# Evaluating MOTIF

For a large-scale empirical evaluation of MOTIF, please follow the instructions given here: https://github.com/juliawal/motif-eval.

# Citing MOTIF

Please cite MOTIF as follows:
