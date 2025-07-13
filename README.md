# MOTIF: Methylation Factor Identifiers

MOTIF is a modular pipeline for inferring gene–CpG regulatory networks from matched expression and methylation data, aggregating multiple runs across parameter grids, and evaluating biological relevance via network-proximity p-values. It aims to identify upstream regulators of aberrant DNA methylation by adapting GRNBoost2 and ensuring robust associations through aggregation of multiple runs.

---

## Repository Structure

```
MOTIF/
├── README.md           
├── environment.yml               - Conda environment specification
│
├── src/
│   ├── motif_core/               - core inference and aggregation
│   │   ├── __init__.py
│   │   ├── infer_grn.py          - `infer_grn()` using modified GRNBoost2
│   │   └── aggregation.py        - `aggregate_single()`, `aggregate_grid()`
│   │ 
│   ├──motif_eval/                - evaluation routines
│   │   ├── __init__.py
│   │   └── evaluation.py         - `run_single_eval()`, `run_grid_eval()`
│   │   
│   └──motif_grnboost/
│      ├── __init__.py
│      ├── algo.py                - from arboreto, adapted
│      ├── core.py                - from arboreto, adapted
│      └── utils.py               - from arboreto
│
├── notebooks/                    - exploratory analyses and visualizations
│   ├── motif.ipynb               - run inference pipeline motif
│   ├── evaluation.ipynb          - calculate p-values for aggregated results from motif
│   ├── parameter_tuning.ipynb    - grid-search p-value analysis
│   └── gene_enrichement.ipynb    - functional enrichment of top-ranked genes
│
└── data/                         - example input templates
```

---

## Installation

Follow these steps to set up the MOTIF pipeline locally without a full package install:

1. **Clone** the repository and enter its root directory:
   ```bash
   git clone https://github.com/your-org/MOTIF.git
   cd MOTIF
   ```

2. **Create** and **activate** the Conda environment:
   ```bash
   conda env create -f environment.yml
   conda activate motif
   ```

3. **Install** the package in editable mode:
   ```bash
   pip install -e .
   ```

---

## Core Pipeline Usage

### 1. Inference

Infer gene-CpG regulatory networks using a modified GRNBoost2 algorithm over multiple random seeds:

- **`infer_grn()`** loads expression and methylation matrices and runs the adapted GRNBoost2 for the given seed.  
- Saves gene-CpGs weight table to `results/inference/grn_output_<seed>`.

Example usage:
```bash
python src/motif_core/grn_inference.py --cpgs_file data/sample_cpgs.tsv --genes_file data/sample_genes.tsv --seed 1

```

#### Required Arguments

| Argument       | Description                                   | Type  |
| -------------- | --------------------------------------------- | ----- |
| `--cpgs_file`  | Path to a TSV file with CpG methylation data. | `str` |
| `--genes_file` | Path to a TSV file with gene expression data. | `str` |
| `--seed`       | Seed for GRNBoost2 run.                       | `int` |

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

Result: `results/grn_inference/grn_with_seed_<seed>.tsv`.

### 2. Aggregation

Aggregate per-seed networks across runs.

#### Single Aggregation

Example usage:
```bash
python src/motif_core/aggregation.py --group_by --method freq --alpha 2 --beta 0.3
```

**Config options (`params/single_agg.yaml`):**

| Argument           | Type             | Default | Description |
|--------------------|------------------|---------|-------------|
| `--nruns`          | int or `'all'`   | `'all'` | How many GRNBoost2 runs to aggregate. `'all'` loads all `.tsv` files from the input directory. |
| `--normalize`      | flag             | `False` | If set, normalize weights by total CpG site weight to avoid bias toward highly connected CpGs. |
| `--group_by`       | flag             | `False` | If set, group by gene and sum weights. Keeps only one entry per gene (top CpG). |
| `--method`         | str              | `'mean'` | Aggregation methods. |
| `--alpha`          | float            | `1.0`   | Used with `'freq'` method. Controls the strength of frequency weighting. Higher values emphasize consistency across runs. |
| `--beta`           | float            | `0.5`   | Used with `'freq'` method. Balances between max and mean weight:<br>• `0` = use only mean<br>• `1` = use only max<br>• `0.5` = equal mix |

- **Aggregation Methods:**  
  - `mean`: average weight across runs.  
  - `borda`: rank-based aggregation.  
  - `freq`: frequency-weighted importance combining mean, max, and frequency; controlled by `alpha` and `beta`.

Result: `results/aggregation/aggregated_weights.tsv`.

#### Grid Search Aggregation

Run grid search over multiple aggregation parameter combinations:

Example usage:
```bash
python src/motif_core/aggregation_grid_search.py
```

Each parameter combination and the resulting aggregations are saved to:
- `results/aggregation_grid_search/params_<parameter_setting_number>.tsv`
- `results/aggregation_grid_search/aggregated_weights_<parameter_setting_number.tsv`

### 3. Evaluation

Assess biological relevance using network-proximity p-values against known regulators (e.g., DNMT/TET genes).

#### Single Evaluation

Computes empirical p-value for top genes of `results/aggregation/aggregated_weights.tsv`:

Example usage:
```bash
python src/motif_eval/evaluation.py
```

#### Grid Search Evaluation

Run evaluation over all aggregated outputs in 
- `results/aggregation_grid_search/params_<parameter_setting_number>.tsv`
- `results/aggregation_grid_search/aggregated_weights_<parameter_setting_number.tsv`
and generate a table of configurations × p-values for top genes of each grid setting:

Example usage:
```bash
python src/motif_eval/evaluation_grid_search.py
```

Results are saved to `results/evaluation_grid_search/grid_search_results.py`

---

## Run Inference Pipeline

add text
```bash
jupyter lab notebooks/motif.ipynb
```

- ..
- ..
- ..

---

## Run Evaluation

add text
```bash
jupyter lab notebooks/evaluation.ipynb
```

- ..
- ..
- ..

---

## Parameter Tuning & Visualization

Explore which parameter settings yield the most significant p-values and robust gene regulators:

```bash
jupyter lab notebooks/parameter_tuning.ipynb
```

- Analyze p-value distributions, consistency across seeds, and generate publication-quality plots.

---

## Gene Enrichement Analysis

Analyze functional enrichment of top-ranked genes identified across multiple aggregation settings:

```bash
jupyter lab notebooks/gene_enrichement_analysis.ipynb
```

- Loads aggregated top genes from results/grid_aggregation/.
- Performs enrichment against GO, KEGG, and Reactome.
- Produces interactive plots and tables of significant terms.

---

## Citation

If you use MOTIF in your research, please cite:

> **Julia Wallnig et al.** “MOTIF: Methylation Factor Identifiers—A Network-Based Approach to Identify Upstream Regulators of Aberrant DNA Methylation in Cancer.” *preprint*, 2025.


