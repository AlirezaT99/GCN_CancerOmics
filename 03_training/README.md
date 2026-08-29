# Stage 3 — Training and biomarker ranking

Trains [MOGONET](https://github.com/txWang/MOGONET) on the stage-2 matrices to classify tumour type,
then ranks genes by ablation importance. MOGONET is cloned rather than vendored, so this stage is a thin
driver over that project's own entry points.

## Setup

```bash
git clone https://github.com/txWang/MOGONET.git
cp -r ICGC MOGONET/          # the stage-2 --out-dir
```

A GPU is strongly recommended — the full run is 500 pretraining epochs plus 2500 training epochs, and the
feature ranking repeats training 5 times.

## Running it

**Using MOGONET directly.** Set the data folder name and `num_class` (stage 2 prints it) in each script,
then:

```bash
python main_mogonet.py      # train and save the models
python main_biomarker.py    # ablate features and rank them
```

**Using `runner.py`.** It wraps the same two calls with the hyperparameters this project used. Note it
exposes two functions and has no `__main__` block, so it is imported and called rather than executed:

```python
from runner import run_mogonet, run_biomarker

run_mogonet(num_class=8)     # substitute the class count stage 2 reported
run_biomarker(num_class=8)
```

Call these in separate processes. Both begin with `os.chdir('./mogonet')`, so calling them in sequence
from one interpreter changes directory twice and the second call fails to find its paths.

## Configuration

| Setting | Value |
| --- | --- |
| `data_folder` | `ICGC` |
| `view_list` | `[1, 2]` — view 1 mutation, view 2 expression |
| `num_epoch_pretrain` | 500 |
| `num_epoch` | 2500 |
| `lr_e_pretrain` | 1e-3 |
| `lr_e` | 5e-4 |
| `lr_c` | 1e-3 |
| repeats for ranking | 5 |

`num_class` is not stored here — it depends on how many cancer types survived stage 2's
`--min-samples` filter, and must be passed in.

## How the ranking works

MOGONET builds one graph convolutional network per omics view plus a cross-omics discovery tensor that
learns label correlations across views. Biomarkers are not read off the network weights: `cal_feat_imp`
ablates each feature in turn and measures the resulting drop in classification performance, so the
features the model depended on most rank highest. `run_biomarker` repeats this over 5 independently
trained models and `summarize_imp_feat` averages the rankings, which damps the variance of any single
training run.

## Outputs

- Trained models under `<data_folder>/models/<rep>/`
- Per-class classification metrics on the held-out test split
- A ranked feature list per omics view — the candidate biomarkers, and the actual result of the project
