# MOTIF: Methylation Factor Identifiers

MOTIF is a modular pipeline for uncovering upstream regulators of aberrant DNA methylation. It prepares matched expression and methylation matrices, runs an adapted GRNBoost2 inference multiple times with different seeds, and aggregates the resulting gene–CpG networks into a single integrated network using user‑configurable aggregation methods. Genes with the highest aggregated edge weights are then prioritized and evaluated by their network‐proximity to known methylation enzymes (DNMT1, DNMT3A/B, TET1/2/3).

---

## Repository Structure

```
MOTIF/
├── README.md
├── environment.yml                     # Conda environment specification
├── requirements.txt                    # pip dependencies
├── pyproject.toml                      # editable install
├── src/
│   ├── motif_core/
│   │   ├── grn_inference.py            # load data and run modified GRNBoost2
│   │   ├── aggregation.py              # aggregate per-seed networks
│   │   └── aggregation_grid_search.py  # parameterized aggregation grid search
│   ├── motif_eval/
│   │   ├── evaluation.py               # single-run network-proximity p-values
│   │   └── evaluation_grid_search.py   # evaluate aggregation grid search results
│   ├── motif_utils/                    # utility modules for aggregation and evaluation
│   │   ├── aggregation_utils.py        # helper functions for aggregation workflows
│   │   └── evaluation_utils.py         # helper functions for evaluation workflows
│   └── motif_grnboost/                 # adapted Arboreto code
│       ├── algo.py
│       ├── core.py
│       └── utils.py
├── notebooks/                          # Jupyter Notebooks for exploratory analyses and visualizations
│   ├── motif.ipynb                     # execute the full inference pipeline and inspect results
│   ├── evaluation.ipynb                # compute network‑proximity p‑values for aggregated results
│   ├── parameter_tuning.ipynb          # grid‑search analysis of aggregation parameters
│   └── gene_enrichment.ipynb           # functional enrichment of top candidate genes
├── data/                               # example input templates
├── results/                            # output directory
└── resources/                          # auxiliary files
```

---

## Installation

Follow these steps to set up the MOTIF pipeline locally:

```bash
git clone https://github.com/your-org/MOTIF.git
cd MOTIF
conda env create -f environment.yml
conda activate motif
pip install -e .
```

> Ensure you unzip the probe manifest before running:
```bash
unzip resources/hippie.txt.zip -d resources/
```

> The repository includes a `.pybiomart.sqlite` file in `resources/`, allowing import of `src` modules directly.

---

## Usage

### 1. Network Inference

Infer gene-CpG regulatory networks using a modified GRNBoost2 algorithm over multiple random seeds:

- **`grn_inference()`** loads expression and methylation matrices and runs the adapted GRNBoost2 for the given seed.
- Saves gene-CpGs weight table to `results/grn_inference/grn_with_seed_<seed>`.

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

#### Optional Arguments

| Argument            | Type           | Default                 | Description                                                                                                         |
| ------------------- | -------------- | ----------------------- | ------------------------------------------------------------------------------------------------------------------- |
| `--nruns`           | int or `'all'` | `'all'`                 | How many GRNBoost2 runs to aggregate. `'all'` loads all `.tsv` files from the input directory.                      |
| `--normalize`       | flag           | `False`                 | If set, normalize weights by total CpG site weight to avoid bias toward highly connected CpGs.                      |
| `--group_by`        | flag           | `False`                 | If set, group by gene and sum weights. Keeps only one entry per gene (top CpG).                                     |
| `--method`          | str            | `'mean'`                | Aggregation method to use: `mean`, `borda`, or `freq`.                                                              |
| `--alpha`           | float          | `1.0`                   | Used with `freq` method. Controls strength of frequency weighting; higher values emphasize consistency across runs. |
| `--beta`            | float          | `0.5`                   | Used with `freq` method. Balances between max and mean weight: `0` = mean only, `1` = max only, `0.5` = equal mix.  |
| `--param_grid_file` | str            | internal parameter grid | JSON file specifying combinations of aggregation parameters to run (overrides default internal grid).               |

- **Aggregation Methods:**
  - `mean`: average weight across runs.
  - `borda`: rank-based aggregation.
  - `freq`: frequency-weighted importance combining mean, max, and frequency; controlled by `alpha` and `beta`.

Result: `results/aggregation/aggregated_weights.tsv`.

#### Grid Search Aggregation

Run grid search over multiple aggregation parameter combinations.

You can supply a parameter grid as a JSON file.

```bash
python src/motif_core/aggregation_grid_search.py --param_grid_file resources/my_params.json
```

If `--param_grid_file` is not provided, a default internal parameter grid will be used.

Each parameter combination and the resulting aggregations are saved to:

- `results/aggregation_grid_search/params_<parameter_setting_number>.tsv`
- `results/aggregation_grid_search/aggregated_weights_<parameter_setting_number>.tsv`

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
- `results/aggregation_grid_search/aggregated_weights_<parameter_setting_number>.tsv`

Example usage:
```bash
python src/motif_eval/evaluation_grid_search.py
```

Results are saved to `results/evaluation_grid_search/grid_search_results.py`

---

## Run Inference Pipeline

```bash
jupyter lab notebooks/motif.ipynb
```

- Loads gene/CpG matrices
- Runs inference over seeds
- Saves GRNs to `results/inference/`

---

## Run Evaluation

```bash
jupyter lab notebooks/evaluation.ipynb
```

- Loads aggregated results
- Computes p-values for top genes
- Visualizes significance and gene overlaps

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

- Loads aggregated top genes from results/grid\_aggregation/.
- Performs enrichment against GO, KEGG, and Reactome.
- Produces interactive plots and tables of significant terms.

---

## Citation

If you use MOTIF in your research, please cite:

> **Julia Wallnig et al.** “MOTIF: Methylation Factor Identifiers—A Network-Based Approach to Identify Upstream Regulators of Aberrant DNA Methylation in Cancer.” *preprint*, 2025.

