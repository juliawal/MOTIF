# MOTIF

MOTIF aims to identify upstream regulators of aberrant DNA methylation by inferring gene–CpG regulatory networks with a modified GRNBoost2 algorithm and aggregating multiple runs for robust associations.

## Installation

1. Clone the repository and enter its root directory:
   ```bash
   git clone https://github.com/yourusername/motif.git
   cd motif
   ```
2. Run the setup script to create a Conda environment and install dependencies:
   ```bash
   bash setup.sh
   ```
3. Unzip the probe manifest file used for mapping genomic coordinates:
   ```bash
   unzip motif_utils/hippie.txt.zip
   ```
---

## Usage

First, activate the environment:
`conda activate motif´

---

### Run MOTIF

#### Required Arguments

| Argument       | Description                                   | Type  |
| -------------- | --------------------------------------------- | ----- |
| `--cpgs_file`  | Path to a TSV file with CpG methylation data. | `str` |
| `--genes_file` | Path to a TSV file with gene expression data. | `str` |

#### Optional Arguments

| Argument                       | Description                                                    | Type   | Default |
| ------------------------------ | -------------------------------------------------------------- | ------ | ------- |
| `--should_filter_cpgs`         | Filter CpGs to gene regions if set.                            | `flag` | *False* |
| `--filtering_genes_file`       | File with Ensembl gene IDs (one per line) for filtering CpGs.  | `str`  | *–*     |
| `--filtering_keep_promoters`   | Keep promoter CpGs during filtering.                           | `flag` | *False* |
| `--filtering_keep_bodies`      | Keep gene body CpGs during filtering.                          | `flag` | *False* |
| `--should_aggregate_cpgs`      | Aggregate CpGs by region (e.g., mean per gene) if set.         | `flag` | *False* |
| `--aggregation_genes_file`     | File with Ensembl gene IDs (one per line) for CpG aggregation. | `str`  | *–*     |
| `--aggregation_keep_promoters` | Include promoter CpGs in aggregation.                          | `flag` | *False* |
| `--aggregation_keep_bodies`    | Include gene body CpGs in aggregation.                         | `flag` | *False* |
| `--only_coding_genes`          | Keep only protein-coding genes during preprocessing.           | `flag` | *False* |
| `--promoter_length`            | Number of base pairs upstream to define promoter regions.      | `int`  | `1000`  |
| `--promoter_overlap`           | Overlap in base pairs between promoters and gene bodies.       | `int`  | `200`   |
| `--run_number`                 | Identifier for GRNBoost2 output files.                         | `int`  | `1`     |

Example usage:
```bash
python motif.py --cpgs_file data/sample_cpgs.tsv --genes_file data/sample_genes.tsv
```

Results are saved in the `grnboost_output/` directory. The primary output file is `grnbboost_output/grn_output_<run_number>.tsv`.

---

## Aggregate edge weights of runs

### Single Aggregation `aggregate_weights.py`

Aaggregate GRNBoost2 outputs using various aggregation strategies.

#### Optional Arguments

| Argument           | Type             | Default | Description |
|--------------------|------------------|---------|-------------|
| `--nruns`          | int or `'all'`   | `'all'` | How many GRNBoost2 runs to aggregate. `'all'` loads all `.tsv` files from the input directory. |
| `--normalize`      | flag             | `False` | If set, normalize weights by total CpG site weight to avoid bias toward highly connected CpGs. |
| `--group_by`       | flag             | `False` | If set, group by gene and sum weights. Keeps only one entry per gene (top CpG). |
| `--method`         | str              | `'mean'` | Aggregation method:<br>• `'mean'` – average weight across runs<br>• `'borda'` – rank-based voting system<br>• `'freq'` – frequency-weighted importance combining mean, max weights and frequency|
| `--alpha`          | float            | `1.0`   | Used with `'freq'` method. Controls the strength of frequency weighting. Higher values emphasize consistency across runs. |
| `--beta`           | float            | `0.5`   | Used with `'freq'` method. Balances between max and mean weight:<br>• `0` = use only mean<br>• `1` = use only max<br>• `0.5` = equal mix |

Example usage:
```bash
python aggregate_weights.py --nruns 5 --normalize --method freq --alpha 2 --beta 0.3
```

Result is saved to `aggregation_output/aggreagted_weights.tsv`

---

### Grid Search Aggregation `aggregate_weights_grid_search.py`

Run **grid search** across multiple aggregation parameter combinations.

| Parameter             | Type   | Description |
|------------------------|--------|-------------|
| `should_normalize`     | bool   | Whether to normalize weights by CpG sums. |
| `should_group_by_gene` | bool   | Whether to sum weights per gene and keep only one entry per gene. |
| `method`               | str    | One of `'mean'`, `'borda'`, `'freq'`. See above. |
| `alpha`                | float  | Only used if `method == 'freq'`. Controls frequency weighting strength. |
| `beta`                 | float  | Only used if `method == 'freq'`. Balances max and mean weights. |

Example usage:
```bash
python aggregate_weights_grid_search.py
```

Each parameter combination and the resulting aggregations are saved to:

- `aggregation_output/params_<index>.tsv`
- `aggregation_output/aggregated_weights_<index>.tsv`

---

### Evaluation evaluate.py

Evaluate the Genes that are involved in the Gene-CpG edges with the highsest aggregated weights and compute **mean p-values** indicating their proximity to known regulators using selected graph proximity measurs.

Example usage:
```bash
python evaluate.py
```

---

### Grid Search Evaluation evaluate_grid_search.py

Automates evaluation of multiple aggregated weights files from grid search.

#### Parameters inside `evaluate_grid_search.py`

Example usage:
```bash
python evaluate_grid_search.py --results results/grid_search/ --output evaluation_report.txt
```

Results are saved to `grid_search_results.py`

## Citation

If you use MOTIF in your research, please cite:

> Wallnig, J. et al. MOTIF: Inferring gene–CpG regulatory networks using adapted GRNBoost2.

