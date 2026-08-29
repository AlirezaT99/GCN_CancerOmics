# Stage 1 — Preprocessing

Cleans the raw ICGC downloads, harmonises gene identifiers onto HGNC symbols, and reduces each cancer
type to the set of `(donor, gene)` pairs that have **both** mutation and expression data.

## Run order

The R scripts must run before the Python scripts. `build_sample_matrices.py` reads the
`converted_genes.csv` that `convert_gene_ids.R` emits, and will silently drop rows for the affected
cancer types if that file is missing.

```bash
Rscript summarize_expression.R                  # all cancer types
Rscript convert_gene_ids.R                      # only where IDs need converting (see below)
python prepare_mutation.py      --cancer-type Brain
python build_sample_matrices.py --cancer-type Brain
```

## Scripts

### `summarize_expression.R`

Reduces each `exp_array.tsv` to three columns — `icgc_donor_id`, `gene_id`,
`normalized_expression_value`. For Blood and Pancreas it first joins a per-cancer
`<Cancer>_gene_mapper.csv` and strips Ensembl version suffixes.

- **In:** `data/<Cancer>/exp_array.tsv` (+ `<Cancer>_gene_mapper.csv` for Blood, Pancreas)
- **Out:** `data/<Cancer>/expression_data.tsv`

> Despite the `.tsv` extension this file is **comma**-separated: the script passes `sep="\t"` to
> `write.csv()`, which ignores it. Stage 2 reads it back with a comma delimiter, so the two are
> consistent — but do not assume the extension is accurate.

### `convert_gene_ids.R`

Maps platform-specific probe identifiers to gene symbols. Detects the scheme from the ID prefix:

| Prefix | Scheme | Cancer type | Mapping |
| --- | --- | --- | --- |
| `ENS…` | Ensembl + version | Blood | strip version, then `bitr()` |
| `NM_` / `NR_` | RefSeq | Nervous System | `bitr()` with `fromType="REFSEQ"` |
| `ILMN_` | Illumina probe | Pancreas | `illuminaHumanv4.db` |

- **In:** `exp_array.tsv` for one cancer type (path set at the top of the script)
- **Out:** `converted_genes.csv` with columns `initial_id`, `Gene`

> The input path and the Illumina branch are both edited by hand — the path is a literal near the top of
> the file, and the Illumina mapping is a commented-out line that gets swapped in for Pancreas. Run this
> once per cancer type needing conversion, moving each `converted_genes.csv` into that cancer's data
> directory.

### `prepare_mutation.py`

Filters the somatic mutation table to single-base substitutions and joins it to the gene list to attach
gene symbols.

- **In:** `data/<Cancer>/simple_somatic_mutation.open.tsv`, `data/genes_list.tsv`
- **Out:** `data/<Cancer>/symbol_mutation.tsv` (genuinely tab-separated)

Read in 1M-row chunks, since the raw mutation tables do not fit comfortably in memory.

### `build_sample_matrices.py`

Intersects the two omics for one cancer type and emits a binary donor × gene matrix marking the cells
where mutation and expression data both exist. This mask is what stage 2 fills with real values.

- **In:** `data/<Cancer>/simple_somatic_mutation.open.tsv`, `data/<Cancer>/exp_array.tsv`,
  `data/<Cancer>/converted_genes.csv` (where applicable), `data/genes_list.tsv`
- **Out:** `output/result-result-<Cancer>.csv.tsv` — see the caveat below
