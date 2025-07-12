# MOTIF: Methylation Factor Identifiers

MOTIF is a modular pipeline for inferring gene–CpG regulatory networks from matched expression and methylation data, aggregating multiple runs across parameter grids, and evaluating biological relevance via network-proximity p-values. It aims to identify upstream regulators of aberrant DNA methylation by adapting GRNBoost2 and ensuring robust associations through aggregation of multiple runs.

---

## Repository Structure

```
MOTIF/
├── README.md           
├── environment.yml               - Conda environment specification
├── setup.py                      - pip-installable via `pip install -e .`
│
├── src/
│   ├── motif_core/               - core inference and aggregation
│   │   ├── __init__.py
│   │   ├── infer_grn.py          - `infer_grn()` using modified GRNBoost2
│   │   └── aggregation.py        - `aggregate_single()`, `aggregate_grid()`
│   └── motif_eval/               - evaluation routines
│       ├── __init__.py
│       └── evaluation.py         - `run_single_eval()`, `run_grid_eval()`
│
├── scripts/                     
│   ├── run_core.py               - inference and aggregation entrypoints
│   └── run_eval.py               - evaluation entrypoint
│
├── notebooks/                    - exploratory analyses and visualizations
│   ├── parameter_tuning.ipynb    - grid-search p-value analysis
│   └── gene_enrichement.ipynb    - functional enrichment of top-ranked genes
│
└── data/                         - example input templates
```

---

## Installation

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

```bash
python scripts/run_core.py \  
  infer \  
  --expr data/sample_expression.tsv \  
  --meth data/sample_methylation.tsv \  
  --outdir results/inference/ \  
  --seeds 0 42 123 2021 777
```

- **`infer_grn()`** loads expression and methylation matrices and runs GRNBoost2 per seed.  
- Saves per-run weight tables to `results/inference/`.

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

### 2. Aggregation

Aggregate per-seed networks across parameter settings.

#### Single Aggregation

```bash
python scripts/run_core.py \  
  aggregate \  
  --indir results/inference/ \  
  --outdir results/aggregation/ \  
  --config params/single_agg.yaml
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

Result: `aggregation_output/aggregated_weights.tsv`.

#### Grid Search Aggregation

Run grid search over multiple aggregation parameter combinations:

```bash
python scripts/run_core.py \  
  aggregate_grid \  
  --indir results/inference/ \  
  --outdir results/grid_aggregation/ \  
  --grid params/grid.yaml
```

Grid YAML example (`params/grid.yaml`):
```yaml
normalization: [false, true]
group_by: [false, true]
method: [mean, borda, freq]
freq:
  alpha: [0.0, 0.5, 1.0, 4.0, 7.0, 10.0]
  beta: [0.0, 0.5, 1.0, 4.0, 7.0, 10.0]
```

- Writes per-combination aggregated tables and parameter logs to `results/grid_aggregation/`.

### 3. Evaluation

Assess biological relevance using network-proximity p-values against known regulators (e.g., DNMT/TET genes).

#### Single Evaluation

```bash
python scripts/run_eval.py \  
  single \  
  --agg-file results/aggregation/aggregated_weights.tsv \  
  --out results/evaluation/pval_single.json
```

- Computes empirical p-values for top genes.

#### Grid Search Evaluation

```bash
python scripts/run_eval.py \  
  grid \  
  --agg-dir results/grid_aggregation/ \  
  --out results/evaluation/pvals_grid.tsv
```

- Run evaluation over all aggregated outputs and generate a table of configurations × p-values for top genes of each grid setting.

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


