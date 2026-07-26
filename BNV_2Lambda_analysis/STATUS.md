# STATUS — BNV 2-Lambda analyses

Snapshot of where things stand. Last updated: **2026-07-25** (Stage 3 PID
selector optimization reviewed by Matt; all recommended selectors, including
the boundary-hugging ones, applied to `channel_config.py`. Ready for Stage 4.)

## Where we are

Stage 3 (PID selector optimization, Punzi FOM) is **complete and checked
in**. Matt reviewed notebook 03, the BAD section, and a new mass-distribution
diagnostic notebook (see below), and decided to apply every recommended
selector, including the boundary-hugging ones (2026-07-25):

- `pid_selector.py` ports `BNV_pLambda/myPIDselector.py`'s bit -> selector-name
  mapping (channel-agnostic), validated bit-for-bit against the reference
  decoder on real signal MC. Decoding itself is simplified to plain vectorized
  bitwise ops (`(bits >> i) & 1`) rather than the reference's binary-string/
  decimal-digit encoding trick -- verified equivalent, not just assumed.
- New PID candidate masks in `cutflow.py`: `get_lambda0_pid_mask` (proton +
  pion KM-ladder cuts, reused for both channels' Lambda0 and for Lam0LamC
  mode 4's inner Lambda0) and `get_lambdac_pid_mask` (mode-dispatched: LambdaC
  daughter-slot-to-particle assignment per mode was verified empirically
  against `LambdaCdNLund` on signal MC, not assumed -- see "Key findings"
  below). Both reuse `get_composite_purity_masks_per_B` unmodified for the
  per-B combination (pass a `{'Lambda0': ..., 'LambdaC': ...}` candidate-mask
  dict, same call as Stage 2 uses for purity).
- `pid_optimization.py`: Punzi FOM (`sig_eff/sqrt(bkg + a/2)`, a=4, verified
  against `BNV_pLambda/PID_study_and_plots_for_BAD.ipynb`) and KM-ladder grid
  scans (`scan_lambda0_pid`, `scan_lambdac_mode_pid`), counting background
  directly in the mES/DeltaE signal region from luminosity-weighted MC (all
  `BACKGROUND_SP_MODES`, not sideband-subtracted -- different from Stage 2's
  method, since this is an MC-truth estimate). Boundary-hugging is flagged by
  reusing `purity_optimization.is_at_scan_boundary` on the ladder's integer
  rung index; a second, PID-specific check flags any recommended point whose
  background estimate rests on very few (<5) raw MC candidates (none did).
- `run_stage03.py --channel all` writes `results/<channel>.yaml` (section
  `stage03`): recommended selector(s) per target, FOM/efficiency/background,
  boundary flag, an explicit "no PID cut at all" FOM comparison for the
  Lambda0 scan, and (Lam0LamC) the multi-candidate study redone with PID
  applied. **Result:** the Lambda0 scan is boundary-hugging at the loosest
  KM rung in *both* channels (though genuinely better than no cut at all);
  3 of 4 LambdaC per-mode scans are boundary-hugging in at least one
  dimension, but mode 2 gives a clean interior optimum and mode 1's
  proton/kaon dimensions do too (only its pion dimension hugs the boundary).
  See the BAD Stage 3 section for full discussion.
- `notebooks/03_pid_optimization.ipynb` executed for both channels (plots in
  `BAD_2Lambda/figures/`); checkpoint cell filled in with the findings above.
- **New: `notebooks/pid_optimization_visualization.ipynb`** (requested by
  Matt at the checkpoint, not tied to a `run_stage*.py` script) -- a visual
  cross-check independent of the Punzi FOM: raw $\Lambda^0$/$\Lambda_c^+$
  mass distributions before vs. after each cut, with the S window and
  contiguous sideband bands shaded (`plotting.plot_mass_cut_diagnostic`,
  same peak/sideband convention as Stage 2's `FOM_SIDEBAND_WIDTH_MULT`).
  Covers the 4 LambdaC-mode PID cuts (`Lam0LamC` only) and the Lambda0
  flight-significance + PID cuts (either channel), each shown for signal MC
  and weighted background MC side by side. This is what settled the
  boundary-hugging cases: LambdaC mode 1's proton/kaon cuts visibly gut the
  background under the peak; the Lambda0 selector's effect is visibly small
  (matching its small FOM margin) but not negligible, so it was kept rather
  than dropped.
- BAD: new subsection in `data_selection.tex` (`sec:stage3pid`) plus a new
  table (`generated_table_pid_selectors`), updated post-review with explicit
  "applied" callouts (matching Stage 2's pattern) and a note on the
  diagnostic-notebook cross-check that justified keeping the boundary cases.
- **Applied to `channel_config.py`'s `pid` section (2026-07-25):**
  `lambda0_selector = {'p': SuperLooseKMProtonSelection, 'pi':
  SuperLooseKMPionMicroSelection}` (both channels); Lam0LamC's
  `lambdac_selector_per_mode`: mode 1 `{p: Tight, K: VeryLoose, pi:
  SuperLoose}`, mode 2 `{p: Loose}`, mode 3 `{p: SuperTight, pi:
  SuperLoose}`, mode 4 `{pi: SuperLoose}` -- i.e. every recommended value
  from the scan, including the boundary-hugging ones.

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
  `mass_window_nsigma`, `pid` section (Lambda0 selector + LambdaC per-mode
  selectors, all `None` pending Stage 3 review). Channels: `Lam0Lam0`
  (B0 -> Lam0 Lam0), `Lam0LamC` (B+ -> LamC+ Lam0).
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
  multi-candidate study and reusable by later stages). Stage 3 PID masks:
  `get_lambda0_pid_mask` (proton+pion KM-selector cut, reused for LambdaC
  mode 4's inner Lambda0), `get_lambdac_pid_mask` (mode-dispatched, fixed
  daughter slots per mode -- verified empirically, see "Key findings").
- `purity_optimization.py` — Stage 2 S/sqrt(S+B) scans
  (`scan_threshold_cut`, `scan_mass_halfwidth`, `sideband_subtracted_counts`,
  with `is_at_scan_boundary` flagging boundary-hugging optima) and the
  LambdaC per-mode Gaussian+linear-background mass fit (`fit_mass_peak`).
- `plotting.py` — stacked SP-mode histograms (signal rescaled per panel to
  background peak, "arb. norm."), mode-split overlays, mES-vs-DeltaE with
  region boxes, Stage 2 `plot_fom_scan` / `plot_mass_fit`, Stage 3
  `plot_pid_ladder_scan_1d`/`_2d` (FOM heatmaps/lines) and
  `plot_mass_cut_diagnostic` (raw mass before/after a cut, S window +
  sidebands shaded -- independent visual check of any purity/PID cut).
- `pid_selector.py` — ported PID SelectorsMap bit decoder (`SELECTORS`
  bit->name mapping per particle, `KM_LADDER` per-particle Super Loose ->
  Super Tight ladder, `bits_for_hypothesis`/`passes_selector`).
- `pid_optimization.py` — Stage 3 Punzi FOM (`punzi_fom`), per-B evaluation
  (`evaluate_combo`), KM-ladder grid scans (`scan_lambda0_pid`,
  `scan_lambdac_mode_pid`), and boundary/low-stats flagging
  (`best_from_ladder_scan`, reusing `purity_optimization.is_at_scan_boundary`).
- `run_stage01.py`, `run_stage02.py`, `run_stage03.py` — write
  `results/<channel>.yaml`.
- `notebooks/01_load_and_diagnostics.ipynb` — full Stage 0/1 diagnostics,
  channel-parametrized, incl. blinding-verification cell.
- `notebooks/02_lambda_purity.ipynb` — Stage 2 scans/fits, cross-check
  plot, multi-candidate study, checkpoint cell. MC only.
- `notebooks/03_pid_optimization.ipynb` — Stage 3 KM-ladder PID scans
  (Lambda0 both channels; LambdaC per mode), no-PID-baseline sanity check,
  multi-candidate study with PID applied, checkpoint cell. MC only.
- `notebooks/pid_optimization_visualization.ipynb` — mass-distribution
  diagnostic companion to notebook 03 (not a `run_stage*.py`-backed stage):
  peak/sideband-shaded Lambda0 and per-mode LambdaC mass, before vs. after
  each applied cut, signal MC and weighted background MC side by side.

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
- **LambdaC daughter slots are fixed per mode** (checked empirically on
  ~2500 signal candidates spanning all 4 modes, not assumed): mode 1
  d1=p/d2=K/d3=pi; mode 2 d1=p/d2=K_S0; mode 3 d1=p/d2=K_S0/d3,d4=pi; mode 4
  d1=Lambda0/d2,d3,d4=pi. No Lund-slot search needed for LambdaC's own
  daughters (unlike the K_S0 gate, which does need one).
- **Stage 3 PID scans**: the Lambda0 proton x pion KM-ladder scan is
  boundary-hugging at the loosest rung (SuperLoose) in *both* channels --
  but an explicit "no PID at all" check confirms SuperLoose is still a real,
  if modest, improvement (FOM 0.364 vs. 0.229 for Lam0Lam0; 0.082 vs. 0.081
  for Lam0LamC, a much smaller margin). Of the 4 LambdaC per-mode scans,
  mode 2 gives a clean interior optimum (Loose proton) and mode 1's
  proton/kaon dimensions do too (Tight proton, Very Loose kaon) -- only its
  pion dimension hugs the boundary; modes 3 and 4 are boundary-hugging in
  every dimension, with mode 3 notable for pulling in opposite directions
  (proton wants tightest, pion wants loosest). None of the recommended
  points rest on <5 raw background-MC candidates (range 6-333), so this
  isn't a low-stats artifact. Redoing the multi-candidate study with these
  PID cuts applied on top of Stage 2 purity (Lam0LamC): good-B fractions
  move from (16.3%, 75.4%, 8.3%) to (32.4%, 62.6%, 5.0%) for (0, 1, >1) --
  PID removes a real chunk of candidates outright, mostly via the LambdaC
  per-mode cuts.

## Open decisions (parked, revisit at the noted stage)

1. **Single-candidate requirement for Lam0LamC** (Stage 1/3-4): current
   `nLambda0 == 1` kills mode 4. Options: mode-aware Lambda0 counting, or
   `nB == 1` only. Stage 3's (now-applied) PID cuts reduce the multi-candidate
   fraction from 8.3% to ~5.0% -- an improvement, not a resolution. Decision
   still deferred, now until after the antibaryon-veto stage.
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

1. Stage 4: antibaryon veto.
