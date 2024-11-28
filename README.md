# Installation

To install the required environment and dependencies, use the provided setup.sh script.

```bash
bash setup.sh
```

# Running MOTIF

To run MOTIF, use the following command:
```bash
conda activate motif
python motif.py --cpgs_file data/sample_cpgs.tsv --genes_file data/sample_genes.tsv
```
The positional arguments are:
```
[1] path to the file containing the CpGs data (e.g., data/sample_cpgs.tsv)
[2] path to the file containing the genes data (e.g., data/sample_genes.tsv).
[3] set to True if you want to filter the CpGs by the provided genes, otherwise False (default is False)
[4] path to a list of specific genes to filter the CpGs for
[5] set to True if you want to aggregate the CpGs data, otherwise False (default is False)
```

# Evaluating MOTIF

For a large-scale empirical evaluation of ROBUST, please follow the instructions given here: https://github.com/juliawal/motif-eval.

# Citing MOTIF

Please cite ROBUST as follows:
