# WORKFLOW-README

How to run the BNV 2-Lambda analyses (B0 -> Lambda0 Lambda0 and
B+ -> LambdaC+ Lambda0) end-to-end. This document is kept current with every
change to the pipeline (see CLAUDE.md ground rules).

## Layout

```
BaBar_BNV_2Lambda_Claude_experiment/
├── data/                     # parquet files (NOT in git)
│   ├── Background_and_signal_SP_modes_All_runs_<channel>.parquet
│   └── Data_All_runs_<channel>_BLINDED.parquet
├── BNV_pLambda/              # completed p-Lambda0 analysis -- REFERENCE ONLY
├── BNV_2Lambda_analysis/     # the new analyses (this is where work happens)
│   ├── channel_config.py     # per-channel config: single source of truth
│   ├── datasets.py           # loading parquet files, MC scaling weights
│   ├── cutflow.py            # mask/cut builders and cutflow tables
│   ├── plotting.py           # histogram creation/filling, standard plots
│   ├── dataset_statistics.csv           # run/skim statistics (from p-Lambda0)
│   ├── SP_cross_sections_and_labels.csv # SP-mode cross sections (from p-Lambda0)
│   └── notebooks/
│       └── 01_load_and_diagnostics.ipynb
└── BAD_2Lambda/              # the analysis note (BAD)
```

`<channel>` is `Lam0Lam0` or `Lam0LamC`. Everything channel-specific
(mass windows, region definitions, candidate structure, histogram
definitions) lives in `channel_config.py`; analysis code takes the config
as an argument.

## Requirements

Python 3 with: `awkward`, `pyarrow`, `numpy`, `pandas`, `matplotlib`, `hist`.

## Blinding

- Only the `_BLINDED` collision-data files are ever read.
  `datasets.load_datasets()` raises if called with `UNBLINDED=True`;
  this stays in place until unblinding is approved (Stage 7).
- Notebook 01 contains a blinding-verification cell: the count of data
  candidates in the assumed signal window must be zero.

## Running the stages

### Stage 0/1 — load data + diagnostics

```
cd BNV_2Lambda_analysis/notebooks
jupyter lab 01_load_and_diagnostics.ipynb
```

- Set `CHANNEL = 'Lam0Lam0'` or `'Lam0LamC'` in the second code cell,
  then run all cells top to bottom.
- Plots are written to `BNV_2Lambda_analysis/notebooks/plots/<channel>/`.
- Review items are listed in the final "Observations / checkpoint summary"
  cell.

*(Later stages will be added here as they are implemented.)*
