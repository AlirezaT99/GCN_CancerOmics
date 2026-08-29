# Stage 2 — Integration

Turns the per-cancer binary masks from stage 1 into the two-view numeric feature matrices that MOGONET
consumes, unions every cancer type into a single feature space, and splits train/test.

```bash
python integrate_matrices.py \
    --data-dir      ./data \
    --matrices-path ./output \
    --out-dir       ICGC \
    --min-samples   10 \
    --train-split   0.70
```

## What it does

1. Reads every stage-1 matrix in `--matrices-path` and derives the cancer label from each filename.
2. Drops cancer types with fewer than `--min-samples` donors.
3. Fills the binary mask with real values, producing two views:
   - **View 1 (mutation)** — the number of mutations per `(donor, gene)`, counted from
     `symbol_mutation.tsv`.
   - **View 2 (expression)** — the `normalized_expression_value` per `(donor, gene)`, read from
     `expression_data.tsv`.
4. Reindexes every cancer type onto the union of all gene columns, zero-filling absent genes, and tags
   each row with its cancer label.
5. Splits rows into train and test, and writes MOGONET's expected file layout.

## Inputs

| Path | From |
| --- | --- |
| `<matrices-path>/result-<Cancer>.tsv` | stage 1 — but see the naming caveat below |
| `<data-dir>/<Cancer>/symbol_mutation.tsv` | `prepare_mutation.py` |
| `<data-dir>/<Cancer>/expression_data.tsv` | `summarize_expression.R` |

## Outputs

Written to `--out-dir`, all headerless and index-free, exactly as MOGONET's loader expects:

| File | Contents |
| --- | --- |
| `1_tr.csv` / `1_te.csv` | mutation view, train / test |
| `2_tr.csv` / `2_te.csv` | expression view, train / test |
| `labels_tr.csv` / `labels_te.csv` | cancer type per row |
| `1_featname.csv` / `2_featname.csv` | gene names per column, one per view |

The script prints the number of classes at the end — MOGONET needs it as `num_class` in stage 3.

## Arguments

| Argument | Default | Notes |
| --- | --- | --- |
| `--data-dir` | `/PROJECTS/Taj/1_PreprocessData/data` | **must be overridden** |
| `--matrices-path` | `/PROJECTS/Taj/1_PreprocessData/output` | **must be overridden**, see below |
| `--out-dir` | `ICGC` | |
| `--min-samples` | `10` | minimum donors for a cancer type to be included |
| `--train-split` | `0.70` | proportion of rows assigned to train |
