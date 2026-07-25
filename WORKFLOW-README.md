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
│   ├── purity_optimization.py # Stage 2: S/sqrt(S+B) scans, LambdaC mode-fit
│   ├── plotting.py           # histogram creation/filling, standard plots
│   ├── run_stage01.py        # Stage 0/1: load + diagnostics -> results/<channel>.yaml
│   ├── run_stage02.py        # Stage 2: purity optimization -> results/<channel>.yaml
│   ├── generate_latex_macros.py # results/<channel>.yaml -> BAD_2Lambda/generated_*.tex
│   ├── dataset_statistics.csv           # run/skim statistics (from p-Lambda0)
│   ├── SP_cross_sections_and_labels.csv # SP-mode cross sections (from p-Lambda0)
│   └── notebooks/
│       ├── 01_load_and_diagnostics.ipynb
│       └── 02_lambda_purity.ipynb
├── results/                  # results/<channel>.yaml -- single source of truth for numbers
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

### Stage 2 — Lambda0 / K_S0 / LambdaC purity optimization

Re-optimizes the $\Lambda^0$ (both channels) and $K_S^0$ (`Lam0LamC`) mass
window + flight-significance cuts with an $S/\sqrt{S+B}$ scan, and fits the
per-mode $\Lambda_c^+$ mass resolution (`Lam0LamC`) to set its per-mode mass
windows. Uses MC only -- no collision data is read.

```
cd BNV_2Lambda_analysis
python run_stage02.py --channel Lam0Lam0
python run_stage02.py --channel Lam0LamC
python run_stage02.py --channel all       # both channels
```

- Writes recommended cut values, scan/fit numbers, and the multi-candidate
  study (Lam0LamC only) to `results/<channel>.yaml` (section `stage02`).
  It does **not** modify `channel_config.py` -- the recommended values are
  printed for review and must be applied to `channel_config.py` by hand
  (or ask Claude) once accepted at the checkpoint.
- **Known issue (flagged, not yet resolved):** several of the flight-
  significance and mass-halfwidth scans return a FOM optimum at the edge
  of the scanned range rather than an interior maximum -- see the BAD
  (`data_selection.tex`, Stage 2 section) and notebook 02's checkpoint
  cell for details before trusting/applying those specific values.
- For the plots and per-mode fit overlays:
  ```
  cd BNV_2Lambda_analysis/notebooks
  jupyter lab 02_lambda_purity.ipynb
  ```
  Set `CHANNEL = 'Lam0Lam0'` or `'Lam0LamC'` in the second code cell, run
  all cells top to bottom. Plots are written to
  `BNV_2Lambda_analysis/notebooks/plots/<channel>/`.

### Refreshing the BAD after any stage

```
cd BNV_2Lambda_analysis
python generate_latex_macros.py   # results/<channel>.yaml -> BAD_2Lambda/generated_*.tex
./copy_plots_to_BAD.sh            # notebooks/plots/ -> BAD_2Lambda/figures/
cd ../BAD_2Lambda
latexmk -pdf main.tex
```

*(Later stages will be added here as they are implemented.)*
