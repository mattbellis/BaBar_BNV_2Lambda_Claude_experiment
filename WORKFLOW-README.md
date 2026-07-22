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
├── results/                  # results/<channel>.yaml -- single source of
│                             #   truth for all analysis numbers (in git)
├── BNV_pLambda/              # completed p-Lambda0 analysis -- REFERENCE ONLY
├── BNV_2Lambda_analysis/     # the new analyses (this is where work happens)
│   ├── channel_config.py     # per-channel config: single source of truth
│   ├── datasets.py           # loading parquet files, MC scaling weights
│   ├── cutflow.py            # mask/cut builders and cutflow tables
│   ├── plotting.py           # histogram creation/filling, standard plots
│   ├── results_io.py         # read/write results/<channel>.yaml
│   ├── run_stage01.py        # Stage 0/1 numbers -> results yaml
│   ├── generate_latex_macros.py  # results yaml -> BAD macros + tables
│   ├── copy_plots_to_BAD.sh  # notebook plots -> BAD_2Lambda/figures/
│   ├── dataset_statistics.csv           # run/skim statistics (from p-Lambda0)
│   ├── SP_cross_sections_and_labels.csv # SP-mode cross sections (from p-Lambda0)
│   └── notebooks/
│       └── 01_load_and_diagnostics.ipynb
└── BAD_2Lambda/              # the analysis note (BAD); build with latexmk
    ├── generated_numbers.tex     # AUTO-GENERATED -- do not edit
    ├── generated_table_*.tex     # AUTO-GENERATED -- do not edit
    └── figures/                  # plots copied in by copy_plots_to_BAD.sh
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

### Stage 0/1 — results files and BAD numbers

All numbers quoted in the BAD come from `results/<channel>.yaml`, which is
written by the per-stage scripts (never hand-typed into LaTeX). The full
documentation chain is:

```
cd BNV_2Lambda_analysis

# 1. Compute the Stage 0/1 numbers and write results/<channel>.yaml
python run_stage01.py --channel all      # or --channel Lam0Lam0 / Lam0LamC

# 2. Regenerate the LaTeX macros + tables in BAD_2Lambda/
python generate_latex_macros.py

# 3. Copy the notebook plots into the BAD figures directory
./copy_plots_to_BAD.sh

# 4. Build the BAD
cd ../BAD_2Lambda
latexmk -pdf main.tex
```

Rerun steps 1-2 whenever the analysis (cuts, config, data) changes, and
step 3 after rerunning the notebooks. The `generated_*.tex` files and
`figures/` in `BAD_2Lambda/` are committed so the BAD always builds.

*(Later stages will be added here as they are implemented.)*
