### Preprocess data

There are a number of scripts involved in this stage.

- `expressionSummarizer.r` prepares expression data columns and Id values
- `geneIdConverter.r` uses a Bioconductor annotation package to unify gene symbols
- `prepare_mutation.py` prepares mutation data by removing useless columns and unwanted rows.
- `preprocessing.py` this code generates a summary for each cancer that will be further used to aggregate in the next stage

by running each python script with a `--help` flag you will learn more about input args.