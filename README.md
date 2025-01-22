# Installation

To install the required environment and dependencies, use the provided setup.sh script.

```bash
bash setup.sh
unzip utils/hippie.txt.zip
```

# Running MOTIF

To run MOTIF, use the following command:
```bash
conda activate motif
python motif.py --cpgs_file data/sample_cpgs.tsv --genes_file data/sample_genes.tsv --should_filter_cpgs --filtering_genes_file utils/sample_filtering_genes.txt --filtering_keep_promoters --filtering_keep_bodies --should_aggregate_cpgs --aggregation_genes_file utils/sample_filtering_genes.txt --aggregation_keep_promoters --aggregation_keep_bodies --only_coding_genes --nruns 10
```

The positional arguments are:
```
[1] path to the file containing the CpGs data (e.g., data/sample_cpgs.tsv)
[2] path to the file containing the genes data (e.g., data/sample_genes.tsv).
[3] set to True if you want to filter the CpGs by the provided genes
[4] path to a list of specific genes to filter the CpGs for
[5] set to True if you want to filter the CpGs in that are the promoters of the given genes
[6] set to True if you want to filter the CpGs in that are the bodies of the given genes
[7] set to True if you want to aggregate the CpGs data for the provided genes
[8] set to True if you want to aggregate the CpGs that are the promoters of the given genes
[9] set to True if you want to aggregate the CpGs in that are in the bodies of the given genes
[10] path to save the outputfile
[11] set to True if you only want to use protein-coding genes
[12] give the Length in base pairs to define the promoter region, default is 1000
[13] give the Length in base pairs to define the overlap region between promoter and gene body, default is 200
[14] number of times to run GRNBoost2 as integer value

### Command-Line Arguments Description

1. **`--cpgs_file`**  
   Path to the file containing the CpGs data (e.g., `data/sample_cpgs.tsv`).  
   **Type**: `str` | **Required**: Yes

2. **`--genes_file`**  
   Path to the file containing the genes data (e.g., `data/sample_genes.tsv`).  
   **Type**: `str` | **Required**: Yes

3. **`--should_filter_cpgs`**  
   Set to `True` if you want to filter the CpGs based on the provided genes.  
   **Type**: `bool` | **Required**: No | **Default**: `False`

4. **`--filtering_genes_file`**  
   Path to the file containing a list of specific genes to filter CpGs (e.g., `data/filter_genes.tsv`).  
   **Type**: `str` | **Required**: No (only required if `--should_filter_cpgs` is set)

5. **`--filtering_keep_promoters`**  
   Set to `True` if you want to keep CpGs that are in the promoter regions of the specified genes during filtering.  
   **Type**: `bool` | **Required**: No | **Default**: `False`

6. **`--filtering_keep_bodies`**  
   Set to `True` if you want to keep CpGs that are in the gene bodies of the specified genes during filtering.  
   **Type**: `bool` | **Required**: No | **Default**: `False`

7. **`--should_aggregate_cpgs`**  
   Set to `True` if you want to aggregate the CpGs data for the provided genes.  
   **Type**: `bool` | **Required**: No | **Default**: `False`

8. **`--aggregation_genes_file`**  
   Path to the file containing a list of genes for aggregation (e.g., `data/aggregate_genes.tsv`).  
   **Type**: `str` | **Required**: No (only required if `--should_aggregate_cpgs` is set)

9. **`--aggregation_keep_promoters`**  
   Set to `True` if you want to aggregate CpGs that are in the promoter regions of the specified genes.  
   **Type**: `bool` | **Required**: No | **Default**: `False`

10. **`--aggregation_keep_bodies`**  
    Set to `True` if you want to aggregate CpGs that are in the gene bodies of the specified genes.  
    **Type**: `bool` | **Required**: No | **Default**: `False`

11. **`--output_dir`**  
    Path to the directory where the output files will be saved.  
    **Type**: `str` | **Required**: No | **Default**: `output`

12. **`--only_coding_genes`**  
    Set to `True` if you only want to use protein-coding genes for the analysis.  
    **Type**: `bool` | **Required**: No | **Default**: `False`

13. **`--promoter_length`**  
    Length in base pairs to define the promoter region.  
    **Type**: `int` | **Required**: No | **Default**: `1000`

14. **`--promoter_overlap`**  
    Length in base pairs to define the overlap region between the promoter and the gene body.  
    **Type**: `int` | **Required**: No | **Default**: `200`

15. **`--nruns`**  
    Number of times to run the GRNBoost2 algorithm.  
    **Type**: `int` | **Required**: No | **Default**: `10`

```

# Evaluating MOTIF

For a large-scale empirical evaluation of MOTIF, please follow the instructions given here: https://github.com/juliawal/motif-eval.

# Citing MOTIF

Please cite MOTIF as follows:
