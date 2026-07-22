# STATUS — BNV 2-Lambda analyses

Snapshot of where things stand. Last updated: **2026-07-22** (Stage 0/1
documentation machinery + BAD sections drafted; awaiting human review).

## Where we are

Stage 0/1 (load data + diagnostics) is **code-complete**, and the
documentation is drafted:

- `results/<channel>.yaml` written by `run_stage01.py` (event counts,
  weights, cutflow, region definitions, blinding check, LambdaC mode
  fractions). Claude test-ran it; **Matt should re-run to verify**:
  `python run_stage01.py --channel all`.
- `generate_latex_macros.py` produces `BAD_2Lambda/generated_numbers.tex`
  (~100 macros) + generated tables (samples, cutflow, LambdaC modes,
  data skims).
- BAD sections drafted: `analysis_introduction.tex` (datasets, MC samples,
  blinding definitions + verification) and `data_selection.tex`
  (multiplicities, preselection cutflow, LambdaC decay modes incl. the
  mode-4 suppression, diagnostic distributions). Title/abstract fixed
  (B+ for the LamC channel); BAD number set to TBD. Builds clean with
  latexmk (26 pp, no undefined refs).

**Checkpoint items for Matt:** review the BAD PDF and the yaml numbers,
re-run run_stage01.py to verify reproducibility, then commit.

## What exists

- `channel_config.py` — single source of truth per channel: composite
  structure, mass windows, mES/DeltaE regions, hist definitions, LambdaC
  decay-mode labels. Channels: `Lam0Lam0` (B0 -> Lam0 Lam0), `Lam0LamC`
  (B+ -> LamC+ Lam0).
- `datasets.py` — `load_datasets` (blinded-only; `UNBLINDED=True` raises),
  MC luminosity weights, `add_derived_fields` (LambdaC flight significance,
  Lam0-from-B vs Lam0-from-LambdaC flight significance).
- `cutflow.py` — single-candidate / fit-region / signal-region / purity
  masks; `get_lambdac_decay_mode` (modes 1-4 from `LambdaCnDaus` +
  `LambdaCd1Lund`); cutflow tables.
- `plotting.py` — stacked SP-mode histograms (signal rescaled per panel to
  background peak, "arb. norm."), mode-split overlays, mES-vs-DeltaE with
  region boxes.
- `notebooks/01_load_and_diagnostics.ipynb` — full Stage 0/1 diagnostics,
  channel-parametrized, incl. blinding-verification cell.

## Key findings so far

- **Blinding**: collision files are blinded upstream; verification cell
  passes. Assumed signal window (|DeltaE| < 0.07, conservative vs. the
  ~0.05 in CLAUDE.md) — final windows still TBD per channel.
- **Data/MC**: shapes agree well; known overall scale offset of order 2x
  (also seen in the p-Lambda0 work).
- **Lam0LamC cutflow** (raw events): signal 112491 -> 57497 (single cand.)
  -> 56213 (+fit region); data 124962 -> 22531 -> 4114.
- **LambdaC decay modes** in signal MC (candidates, no cuts): mode 1
  (pKpi) 31.7%, mode 2 (pKS) 22.9%, mode 3 (pKSpipi) 22.2%, mode 4
  (Lam0pipipi) 23.2%. **The `nLambda0 == 1` single-candidate requirement
  suppresses mode 4** (its Lam0 shares the collection with the B's Lam0).
- **K_S mass**: use `K_SpreFitMass` (7.8 MeV resolution); `K_SMass` is
  post-fit mass-constrained (0.3 MeV spike).
- LambdaC has no post-fit flight-significance branch; using
  `FlightLen/FlightErr` (median ~1.6).
- Candidate multiplicities in `Lam0LamC` reach ~60 B/event (vs. <=5 in
  `Lam0Lam0`) — single-candidate treatment matters much more there.

## Open decisions (parked, revisit at the noted stage)

1. **Single-candidate requirement for Lam0LamC** (Stage 1/3-4): current
   `nLambda0 == 1` kills mode 4. Options: mode-aware Lambda0 counting, or
   `nB == 1` only. Decision deferred until we see how PID/purity cuts
   reduce *surviving* candidate multiplicity.
2. **Exact signal windows per channel** (Stage 1/2): must stay no narrower
   than the upstream blinding; verify against the blinded box.
3. **LambdaC mass window** (Stage 2): placeholder ±10 MeV; set from purity
   studies.
4. **LambdaC / K_S flight cuts** (Stage 2): thresholds TBD (K_S window on
   the pre-fit mass; Lam0 window ±3 MeV from p-Lambda0 to be re-checked).
5. **Which integrated luminosity to quote**: the dataset CSV sums to
   430.9 fb^-1 (macro \IntLumiInvFb), but CLAUDE.md (and the p-Lambda0
   abstract) quote 424.3 fb^-1 at the Y(4S) — likely on-peak vs. total.
   *Decision (Matt, 2026-07-22): keep 430.9 for now, may revisit.*

## Next steps (agreed)

1. **Documentation for Stage 0/1** (next session): set up
   `results/<channel>.yaml` + LaTeX-macro generator (single source of
   truth for numbers), then draft the BAD sections (datasets/samples,
   multiplicities, diagnostics) in `BAD_2Lambda/`, using the p-Lambda0 BAD
   as the structural template.
2. Stage 1 formal checkpoint / commit, then Stage 2 (Lambda purity).
