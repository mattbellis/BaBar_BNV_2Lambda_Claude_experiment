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
│   ├── pid_selector.py       # PID SelectorsMap bit decoding (ported from BNV_pLambda)
│   ├── pid_optimization.py   # Stage 3: Punzi-FOM KM-ladder PID scans
│   │                         #   + Stage 4: antibaryon-veto ladder scan
│   ├── cumulative_performance.py # Stage 4 Phase 3: cumulative Stage 1-4 cutflow
│   ├── plotting.py           # histogram creation/filling, standard plots
│   ├── run_stage01.py        # Stage 0/1: load + diagnostics -> results/<channel>.yaml
│   ├── run_stage02.py        # Stage 2: purity optimization -> results/<channel>.yaml
│   ├── run_stage03.py        # Stage 3: PID optimization -> results/<channel>.yaml
│   ├── run_stage04.py        # Stage 4: antibaryon veto -> results/<channel>.yaml
│   ├── generate_latex_macros.py # results/<channel>.yaml -> BAD_2Lambda/generated_*.tex
│   ├── dataset_statistics.csv           # run/skim statistics (from p-Lambda0)
│   ├── SP_cross_sections_and_labels.csv # SP-mode cross sections (from p-Lambda0)
│   └── notebooks/
│       ├── 01_load_and_diagnostics.ipynb
│       ├── 02_lambda_purity.ipynb
│       ├── 03_pid_optimization.ipynb
│       ├── 04_antibaryon_veto.ipynb
│       └── pid_optimization_visualization.ipynb  # mass-diagnostic companion to 03
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
- Stages 2 and 3 read no collision data at all. Stage 4 does, but only for
  the cumulative cutflow's data column, and counts it in the fit region with
  the signal box explicitly removed — see the Stage 4 section below.

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

### Stage 3 — PID selector optimization

Optimizes the KM-family PID selector (Super Loose -> Super Tight) for the
$\Lambda^0$ proton/pion daughters (both channels) and, for `Lam0LamC`, each
$\Lambda_c^+$ decay mode's own daughters, with a Punzi figure of merit
($\mathrm{FOM} = \epsilon_{\rm sig}/\sqrt{B + a/2}$, $a=4$). Uses MC only --
no collision data is read.

```
cd BNV_2Lambda_analysis
python run_stage03.py --channel Lam0Lam0
python run_stage03.py --channel Lam0LamC
python run_stage03.py --channel all       # both channels
```

- Writes the recommended selector(s), FOM/efficiency/background at that
  operating point, the boundary-hugging flag, and (Lam0LamC only) the
  multi-candidate study redone with PID applied, to `results/<channel>.yaml`
  (section `stage03`). It does **not** modify `channel_config.py` -- the
  recommended values are printed for review and must be applied to
  `channel_config.py`'s `pid` section by hand (or ask Claude) once accepted
  at the checkpoint.
- **Known result (not a bug, see the BAD Stage 3 section for details):** the
  $\Lambda^0$ PID scan is boundary-hugging at the loosest KM rung in both
  channels, and 3 of the 4 LambdaC per-mode scans are boundary-hugging in at
  least one dimension -- `results/<channel>.yaml`'s `stage03` section
  records an explicit "no PID cut at all" comparison FOM for the Lambda0
  scan to check whether the recommended (loosest) selector is still a real
  improvement over no cut. **Applied to `channel_config.py`'s `pid` section
  2026-07-25** (all recommended selectors, including the boundary-hugging
  ones -- see the mass-distribution diagnostic below for why those were
  kept rather than dropped).
- For the FOM scan plots (PID scan heatmaps/lines, multi-candidate
  comparison):
  ```
  cd BNV_2Lambda_analysis/notebooks
  jupyter lab 03_pid_optimization.ipynb
  ```
  Set `CHANNEL = 'Lam0Lam0'` or `'Lam0LamC'` in the second code cell, run
  all cells top to bottom. Plots are written to
  `BNV_2Lambda_analysis/notebooks/plots/<channel>/`.
- For a visual cross-check independent of the FOM (raw $\Lambda^0$/$\Lambda_c^+$
  mass distributions before vs. after each cut, with the S window and
  sideband bands shaded -- same peak/sideband convention as Stage 2):
  ```
  cd BNV_2Lambda_analysis/notebooks
  jupyter lab pid_optimization_visualization.ipynb
  ```
  Set `CHANNEL` in the second code cell (the $\Lambda_c^+$ per-mode section
  only runs for `Lam0LamC`; the $\Lambda^0$ sections run for either
  channel). Not tied to a `run_stage*.py` script -- diagnostic plots only,
  no `results/<channel>.yaml` writes.

### Stage 4 — antibaryon veto + cumulative Stage 1–4 performance

Vetoes a $B$ candidate when the event contains an identified antiproton that
is not one of that candidate's own tracks (the BNV signal final state is
all-baryon; SM decays producing a baryon also produce an antibaryon). The
KM-family proton selector defining "identified antiproton" is optimized with
the same Punzi FOM ($a=4$) as Stage 3, on top of the Stage 2 purity + Stage 3
PID selection. Also produces the cumulative Stage 1–4 cutflow.

```
cd BNV_2Lambda_analysis
python run_stage04.py --channel Lam0Lam0
python run_stage04.py --channel Lam0LamC
python run_stage04.py --channel all       # both channels
```

- Writes the recommended selector, FOM/efficiency/background, boundary flag,
  the comparison points (no veto at all; the p-Lambda0 analysis's
  `TightKMProtonSelection`; the reference-faithful exclusion scope; the
  p-list-only track pool; the anti-Lambda0 add-on), the per-LambdaC-mode veto
  efficiency, the cumulative cutflow, and the per-step candidate
  multiplicities to `results/<channel>.yaml` (section `stage04`). It does
  **not** modify `channel_config.py` — the recommended selector is printed for
  review and must be applied by hand (or ask Claude) to
  `channel_config.py`'s `antibaryon_veto` section once accepted at the
  checkpoint. Until then `antibaryon_veto['selector'] = None`, which means
  **no veto is in force** for any other code.
- **Unlike Stages 2 and 3, this reads the BLINDED collision file** (for the
  cumulative cutflow's data column). Data is counted in the fit region with
  the signal box **explicitly removed** (`get_fit_mask & ~get_signal_region_mask`);
  the signal box is never read. The FOM scan itself is MC-only.
- **Note the ladder direction is reversed relative to Stage 3:** a *looser*
  proton selector makes the veto fire *more* often — more background
  rejection, more signal loss.
- **Known limitation (Lam0LamC modes 2 and 3):** the K_S0's two pion tracks
  cannot be resolved from these parquet files (`K_Sd1Idx`/`K_Sd2Idx` do not
  index any collection present in the file), so they stay eligible to fire
  the veto and cost ~1% of signal efficiency in those two modes. Measured per
  mode and recorded in `results/Lam0LamC.yaml`; see the note in
  `cutflow.get_signal_b_track_slots` and the Stage 4 BAD section.
- For the scan plot, the before/after diagnostics (mES, DeltaE, and a
  Lambda0-mass null check), the cumulative cutflow table and the
  multiplicity comparison:
  ```
  cd BNV_2Lambda_analysis/notebooks
  jupyter lab 04_antibaryon_veto.ipynb
  ```
  Set `CHANNEL = 'Lam0Lam0'` or `'Lam0LamC'` in the second code cell, run all
  cells top to bottom. Plots are written to
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
