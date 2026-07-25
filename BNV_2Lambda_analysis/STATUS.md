# STATUS — BNV 2-Lambda analyses

Snapshot of where things stand. Last updated: **2026-07-24** (Stage 2
reviewed by Matt; non-boundary-hugging recommended values applied to
`channel_config.py`. Ready for Stage 3.)

## Where we are

Stage 2 (Lambda0/K_S0/LambdaC purity optimization) is **complete and
checked in**. Matt reviewed notebook 02 and the BAD section; of the six
Stage 2 cuts, the two that showed a genuine interior FOM maximum plus the
four LambdaC per-mode windows were applied to `channel_config.py`
(2026-07-24):

- **Applied:** Lam0LamC's Lambda0 flight-significance cut (25 -> 18);
  Lam0LamC's K_S0 mass window (placeholder +/-10 MeV -> +/-26 MeV); all
  four LambdaC per-mode mass windows (mass_window_nsigma=3 around each
  mode's fitted mu/sigma -- see `_lambdac_mode_window(...)` calls in
  `channel_config.py`).
- **NOT applied (boundary-hugging scans, see "Known issue" below):**
  Lam0Lam0's Lambda0 flight cut and mass window (stay at the p-Lambda0
  defaults, 25 / +/-3 MeV); Lam0LamC's Lambda0 mass window (stays at
  +/-3 MeV); Lam0LamC's K_S0 flight cut (stays at the Stage 1 placeholder,
  5.0).
- `run_stage01.py --channel all` and `run_stage02.py --channel all`
  re-run after the config change (results/<channel>.yaml refreshed);
  `generate_latex_macros.py` re-run; BAD rebuilds clean with latexmk.
- `notebooks/02_lambda_purity.ipynb` executed for both channels (plots in
  `BAD_2Lambda/figures/`); Matt's review comments are in its checkpoint
  cell.
- BAD: new subsection in `data_selection.tex` (`sec:stage2purity`) plus a
  new table (`generated_table_lambdac_mode_fits`), with explicit
  "applied" / "not applied" callouts matching the config state above.

**Known issue -- flagged, still open (Stage 3 does not depend on it, but
resolve before trusting/using the affected cuts):** the four NOT-applied
scans above return an optimum at the *edge* of the scanned range rather
than an interior maximum -- the FOM is monotonic across the whole grid.
Most likely cause: the input ntuples are already skimmed to a narrow band
around the resonance mass, so the sideband bands run out of background
before the scan range does, biasing $S/\sqrt{S+B}$ toward the
loosest/widest setting rather than reflecting a real trade-off. The two
*applied* scans (Lam0LamC's Lambda0 flight cut and K_S0 mass window) show
clean interior maxima, so the method itself is sound -- it's specifically
these four that need either a wider scan range or a methodology change
(e.g. restricting the scan to candidates already in the mES/DeltaE fitting
region). Matt: revisit at Stage 2's "revisit" note below, or whenever it
becomes relevant (e.g. if Lam0Lam0's efficiency looks low in Stage 3+).

Stage 0/1 (load data + diagnostics) is **code-complete** and documented
(previous checkpoint):

- `results/<channel>.yaml` written by `run_stage01.py` (event counts,
  weights, cutflow, region definitions, blinding check, LambdaC mode
  fractions). Reproducibility verified (Matt re-ran notebook 01; Claude
  re-ran `run_stage01.py --channel all` after the Stage 2 config change).
- `generate_latex_macros.py` produces `BAD_2Lambda/generated_numbers.tex`
  (~100 macros) + generated tables (samples, cutflow, LambdaC modes,
  data skims).
- BAD sections drafted: `analysis_introduction.tex` (datasets, MC samples,
  blinding definitions + verification) and `data_selection.tex`
  (multiplicities, preselection cutflow, LambdaC decay modes incl. the
  mode-4 suppression, diagnostic distributions). Title/abstract fixed
  (B+ for the LamC channel); BAD number set to TBD. Builds clean with
  latexmk (26 pp, no undefined refs).

**Checkpoint items for Matt:** done -- reviewed notebook 02 and re-ran
notebook 01 himself (no significant change), left review comments in
notebook 02's checkpoint cell, and confirmed which values to apply (table
above). Claude applied them to `channel_config.py` and re-ran
`run_stage01.py`/`run_stage02.py`/`generate_latex_macros.py` to refresh
`results/<channel>.yaml` and the BAD. Remaining open item is the
boundary-scan issue (parked, see above); not a blocker for Stage 3.

## What exists

- `channel_config.py` — single source of truth per channel: composite
  structure, mass windows, mES/DeltaE regions, hist definitions, LambdaC
  decay-mode labels, Stage 2 scan ranges (`flight_scan`,
  `mass_halfwidth_scan`), `FOM_SIDEBAND_WIDTH_MULT`, `k0s` config
  (Lam0LamC), `LambdaC` per-mode `mass_windows_per_mode` +
  `mass_window_nsigma`. Channels: `Lam0Lam0` (B0 -> Lam0 Lam0), `Lam0LamC`
  (B+ -> LamC+ Lam0).
- `datasets.py` — `load_datasets` (blinded-only; `UNBLINDED=True` raises),
  MC luminosity weights, `add_derived_fields` (LambdaC flight significance,
  Lam0-from-B vs Lam0-from-LambdaC flight significance).
- `cutflow.py` — single-candidate / fit-region / signal-region masks;
  `get_lambdac_decay_mode` (modes 1-4 from `LambdaCnDaus` +
  `LambdaCd1Lund`); cutflow tables; Stage 2 purity-mask builders:
  `get_purity_mask_for_comp` (generic mass+flight), `get_lambdac_ks_info`/
  `get_lambdac_k0s_gate` (K_S0<->LambdaC linkage via daughter Lund==310),
  `get_lambdac_purity_mask` (mode window + K_S0 gate),
  `get_composite_purity_masks_per_B` (per-B daughter check, used by the
  multi-candidate study and reusable by later stages).
- `purity_optimization.py` — Stage 2 S/sqrt(S+B) scans
  (`scan_threshold_cut`, `scan_mass_halfwidth`, `sideband_subtracted_counts`,
  with `is_at_scan_boundary` flagging boundary-hugging optima) and the
  LambdaC per-mode Gaussian+linear-background mass fit (`fit_mass_peak`).
- `plotting.py` — stacked SP-mode histograms (signal rescaled per panel to
  background peak, "arb. norm."), mode-split overlays, mES-vs-DeltaE with
  region boxes, Stage 2 `plot_fom_scan` / `plot_mass_fit`.
- `run_stage01.py`, `run_stage02.py` — write `results/<channel>.yaml`.
- `notebooks/01_load_and_diagnostics.ipynb` — full Stage 0/1 diagnostics,
  channel-parametrized, incl. blinding-verification cell.
- `notebooks/02_lambda_purity.ipynb` — Stage 2 scans/fits, cross-check
  plot, multi-candidate study, checkpoint cell. MC only.

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
- **Stage 2 FOM scans**: several boundary-hugging optima (see "Known
  issue" above) -- likely an upstream mass-skim depleting the sidebands.
  The one clean interior result (Lam0LamC's Lambda0 flight-significance
  cut, recommended 18, FOM 321) sits close to the p-Lambda0 default of 25,
  which is reassuring for the method.
- **LambdaC per-mode mass resolution** (signal MC, `K_S0` gate applied):
  all four modes fit cleanly (Gaussian + linear background) with
  sigma = 4.2/5.0/4.3/5.2 MeV for modes 1-4 respectively -- modes 2 and 4
  (fewer/different daughters) are visibly wider than 1 and 3. Fitted means
  are ~1.5-1.7 MeV below the PDG value (2284.8-2285.0 vs. 2286.46) in all
  four modes -- small, consistent, not yet investigated.
- **Multi-candidate study** (Lam0LamC, after all Stage 2 purity cuts):
  16.3% of signal events have 0 good B candidates, 75.4% have exactly 1,
  8.3% still have >1 -- a large reduction from the Stage 1 baseline
  (~40% single-candidate for modes 1/2, ~20% for mode 3, 0% for mode 4),
  though still no candidate-selection policy applied.

## Open decisions (parked, revisit at the noted stage)

1. **Single-candidate requirement for Lam0LamC** (Stage 1/3-4): current
   `nLambda0 == 1` kills mode 4. Options: mode-aware Lambda0 counting, or
   `nB == 1` only. Decision deferred until we see how PID cuts further
   reduce *surviving* candidate multiplicity (Stage 2's multi-candidate
   numbers above are a first look, still 8.3% multi-candidate).
2. **Exact signal windows per channel** (Stage 1/2): must stay no narrower
   than the upstream blinding; verify against the blinded box.
3. **Stage 2 FOM-scan boundary issue** (Stage 2, still open): 4 of 6
   scanned cuts hit the edge of their scan range (see "Known issue"
   above) and were NOT applied to `channel_config.py`; they stay at their
   Stage 1 defaults. Needs a decision (wider scan range vs.
   fit-region-restricted rescan vs. other) whenever it becomes a priority
   -- not a Stage 3 blocker, since Stage 3 (PID) is a separate cut applied
   on top of whatever purity cuts are in force.
4. **Which integrated luminosity to quote**: the dataset CSV sums to
   430.9 fb^-1 (macro \IntLumiInvFb), but CLAUDE.md (and the p-Lambda0
   abstract) quote 424.3 fb^-1 at the Y(4S) — likely on-peak vs. total.
   *Decision (Matt, 2026-07-22): keep 430.9 for now, may revisit.*

## Next steps (agreed)

1. Stage 3: PID selector optimization (Punzi figure of merit) for the
   proton (both channels' Lambda0 -> p pi-) and, for Lam0LamC, the
   LambdaC daughters (p, K, pi per mode) and the K_S0 pions (K_S0 decays
   weak/neutral, no PID needed on the K_S0 itself, but its pi+pi-
   daughters may still benefit from a loose PID veto -- TBD). See "Stage
   3 pickup notes" below.

## Stage 3 pickup notes (PID selector optimization)

Enough to start Stage 3 cold, without re-deriving from scratch:

- **PID branches exist** in both channels' parquet files:
  `pSelectorsMap`, `KSelectorsMap`, `piSelectorsMap`, `muSelectorsMap`,
  `eSelectorsMap` -- one per charged-track collection, each an integer
  bitmap (one bit per named selector, e.g.
  `TightLHProtonSelection`/`VeryTightKMProtonSelection` for protons,
  `TightBDTKaonMicroSelection` for kaons, etc.).
- **Reference decoder**: `BNV_pLambda/myPIDselector.py` (`PIDselector`
  class) maps bit position -> selector name per particle type
  (`particle='p'/'K'/'pi'/...`). Not yet ported into
  `BNV_2Lambda_analysis/`; porting it (channel-agnostic, so it belongs
  next to `datasets.py`/`cutflow.py`, not duplicated) is the natural
  first step of Stage 3.
- **Reference PID-optimization notebooks** (read for the FOM/scan
  pattern, not the exact code -- they predate this channel-parametrized
  framework): `BNV_pLambda/PID_study_and_plots_for_BAD.ipynb`,
  `PID_selector_function_explainer.ipynb`, `trying_out_PID_selector.ipynb`.
- **Reference Punzi FOM** (the actual PID-selector one; verified by
  reading the notebook, not guessed): `BNV_pLambda/PID_study_and_plots_for_BAD.ipynb`,
  `fom = sig_eff / sqrt(bkg998 + bkg1005 + a/2)` with `a = 4` (Punzi 2003
  form, eff / (a/2 + sqrt(B)) -- appropriate for optimizing toward an
  upper limit rather than a discovery significance, per CLAUDE.md).
  `sig_eff` there is signal MC efficiency relative to the pre-PID sample;
  background is a weighted sum of specific continuum SP modes (998, 1005)
  using ad hoc Run1-era scale factors in the reference -- in this
  framework, use `datasets.get_scaling_weights` instead of copying those
  factors (same substitution already made for Stage 2's background
  weighting). NOTE: `babar_analysis_tools.py`'s `punzi_fom_nn` is a
  *different* stage's reference (the DeltaE-sideband MLP-output Punzi FOM
  for Stage 5, not Stage 3) -- don't reuse it here by mistake.
  Note the Stage 2 boundary-hugging issue above: if a PID-threshold scan
  runs into the same skimmed-sideband problem, the fix path is the same
  (widen range / restrict to fit region) -- but PID scans are typically
  over a discrete selector ladder (Very Loose -> ... -> Very Tight), not
  a continuous variable, so this may not apply the same way.
- **What Stage 3 selects PID for**: the proton in every Lambda0 -> p pi-
  (both channels); for `Lam0LamC`, additionally the LambdaC daughters,
  which differ **by decay mode** (`cutflow.get_lambdac_decay_mode`):
  mode 1 (p K- pi+) needs p+K+pi selectors, mode 2 (p K_S0) needs only a
  p selector (K_S0's pi+pi- are from a V0, typically not PID-cut), mode 3
  (p K_S0 pi+pi-) needs p (+ loose pi), mode 4 (Lambda0 pi+pi+pi-) needs
  the *inner* Lambda0's proton (already covered) + pi selectors on the
  three pions.
- **Architecture expectation**: following the existing pattern (Stage 2's
  `get_purity_mask_for_comp` / `get_lambdac_purity_mask`), PID cuts should
  be new candidate-level mask builders in `cutflow.py`, parametrized by
  channel config (new PID selector choices + thresholds probably belong
  in a new `channel_config.py` section, e.g. `pid` or per-composite
  `pid_selector`), reusing `get_composite_purity_masks_per_B` (or
  extending it) for the per-B combination. `purity_optimization.py`'s
  scan/FOM helpers are Stage-2-specific (S/sqrt(S+B)); Stage 3 needs a
  Punzi-FOM analog -- probably a new module or an addition to
  `purity_optimization.py` (name TBD, e.g. `pid_optimization.py`, if the
  scan mechanics differ enough to not share code cleanly).
- **Follow CLAUDE.md rule 4**: PID selector choices (which map/threshold
  per particle, per mode) are exactly the kind of nontrivial design
  decision to propose and agree with Matt before implementing, not to
  decide unilaterally.
