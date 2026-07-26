# STATUS — BNV 2-Lambda analyses

Snapshot of where things stand. Last updated: **2026-07-26** (Stage 4
antibaryon veto **run and documented**; awaiting Matt's review decisions.)

## >>> PICK UP HERE (read this first after a break) <<<

Stage 4 is **done and run**. `run_stage04.py --channel all`,
`generate_latex_macros.py` and `latexmk` have all been executed;
`results/<channel>.yaml` section `stage04` is populated and the BAD builds
clean (exit 0, no undefined refs, 38 pp). **Nothing is committed** -- the
working tree holds all of Stage 4 plus earlier uncommitted edits.

**Nothing has been applied to `channel_config.py`**:
`antibaryon_veto['selector']` is still `None`, i.e. **no veto is in force**
for any other code.

**Four decisions are waiting on Matt** (details + numbers below):

1. **Which veto selector to apply per channel.** Scan says
   `TightKMProtonSelection` (Lam0Lam0) and `LooseKMProtonSelection`
   (Lam0LamC). Caveat: the Lam0LamC optimum is *shallow and non-monotonic*
   (FOM 0.3579 at Loose, dips to 0.3530 at Tight, back up to 0.3552 at
   SuperTight) -- so adopting the p-Lambda0 `Tight` for BOTH channels, for
   cross-analysis consistency, costs only 1.4% of FOM in Lam0LamC and
   nothing at all in Lam0Lam0 (where Tight *is* the optimum). That is a
   defensible alternative to the literal scan output.
2. **Anti-Lambda0 add-on**: recommend leaving OFF. Measured: no FOM gain in
   either channel (Lam0LamC 0.3556 vs 0.3579; Lam0Lam0 bit-identical, i.e.
   fully redundant).
3. **Single-candidate requirement** (open decision 1 below). This was parked
   "until after the antibaryon-veto stage" -- that has now happened, so it
   is ready to decide. The multi-candidate fraction is no longer the main
   problem (4.4%); the 39.7% **zero**-good-B fraction is.
4. **Stale `stage02` numbers** -- re-run `run_stage02.py --channel all`
   (see item 5 under "Open decisions").

**Also not done:** `notebooks/04_antibaryon_veto.ipynb` has been written but
**never executed**, so it has no output cells and its plots
(`antibaryon_veto_scan.png`, `..._mes_diagnostic.png`,
`..._deltae_diagnostic.png`, `..._lambda0_mass_diagnostic.png`) do not exist
yet and are not in the BAD. Run it for both channels, then
`./copy_plots_to_BAD.sh`. (VS Code: use "Revert File" first if it is open.)

### Stage 4 final numbers (from `results/<channel>.yaml`, section `stage04`)

Veto scan (Punzi, a=4, on top of Stage 2 purity + Stage 3 PID). Neither
optimum is boundary-hugging and neither rests on low MC stats:

| | Lam0Lam0 | Lam0LamC |
|---|---|---|
| recommended | `TightKMProtonSelection` | `LooseKMProtonSelection` |
| signal eff | 0.9096 | 0.8934 |
| weighted bkg | 1.75 | 4.23 |
| FOM | 0.4698 | 0.3579 |
| FOM, no veto at all | 0.3920 | 0.2654 |
| FOM at p-Lambda0 `Tight` | 0.4698 (same point) | 0.3530 |

Full ladders (selector, eff, bkg, FOM):

- Lam0Lam0: SuperLoose .8886/1.75/.4589 | VeryLoose .8958/1.75/.4626 |
  Loose .9058/1.75/.4678 | **Tight .9096/1.75/.4698** |
  VeryTight .9121/2.00/.4561 | SuperTight .9131/2.25/.4429.
  Background is *flat* at 1.75 over the first four rungs, so the FOM there
  is driven purely by signal efficiency; it only turns over once bkg rises.
- Lam0LamC: SuperLoose .8650/4.23/.3465 | VeryLoose .8772/4.23/.3514 |
  **Loose .8934/4.23/.3579** | Tight .8987/4.48/.3530 |
  VeryTight .9028/4.48/.3546 | SuperTight .9042/4.48/.3552. Note the
  non-monotonicity -- the "optimum" is a shallow local max.

Cumulative Stage 1-4 (signal/bkg = B candidates in the signal box; data =
fit region EXCLUDING the signal box):

| Step | Lam0Lam0 sig / eff / bkg / data | Lam0LamC sig / eff / bkg / data |
|---|---|---|
| 0 presel + single cand | 42658 / 1.000 / 272.09 / 7831 | 53801 / 1.000 / 142.32 / 4114 |
| 1 + fit region | 42658 / 1.000 / 272.09 / 7831 | 53801 / 1.000 / 142.32 / 4114 |
| 2 + Stage 2 purity | 30212 / 0.708 / 16.99 / 437 | 43398 / 0.807 / 12.72 / 412 |
| 3 + Stage 3 PID | 28024 / 0.657 / 4.51 / 151 | 38303 / 0.712 / 2.75 / 81 |
| 4 + Stage 4 veto | 25498 / **0.598** / **1.75** / 59 | 34425 / **0.640** / **1.49** / 34 |

(Step 1 is a no-op for all three columns *by construction* -- the signal box
is a subset of the fit region, and the data column already has the fit
region imposed. Kept only so the pipeline order is complete.)

Per-event good-B fractions (0 / 1 / >1), Stage 2/3/4 convention -- no
single-candidate cut, no region cut:

- Lam0Lam0: 0.305/0.695/0.000 -> 0.356/0.644/0.000 -> **0.414/0.586/0.000**
- Lam0LamC: 0.232/0.700/0.068 -> 0.324/0.626/0.050 -> **0.397/0.559/0.044**

Per-LambdaC-mode veto efficiency (quantifies the K_S0 caveat below):
mode 1 **0.9004**, mode 2 0.8880, mode 3 0.8910, mode 4 **0.8934**.
Modes 1 and 4 have no K_S0 and are the clean comparison.

## Where we are

Stage 4 (antibaryon veto + cumulative Stage 1-4 performance) is **complete
and run** -- numbers and pending decisions are in the "PICK UP HERE" block
above. Implementation detail follows.

- **Veto definition**: per-B-candidate (the p-Lambda0 reference is per-event
  and explicitly assumes one B/event; Lam0LamC reaches ~60). A candidate is
  vetoed if the event holds a track that (a) is not one of *that*
  candidate's tracks, (b) passes the chosen KM proton selector, and (c) has
  charge opposite to that candidate's signal baryon charge, read off its own
  Lambda0 proton daughter (`sign(Lambda0d1Lund)`) so B/Bbar flip
  automatically. Exclusion scope is **all** of the candidate's tracks, not
  just the proton daughters as in the reference -- the signal Lambda0 pions
  carry exactly the antiproton charge, so leaving them eligible lets the
  signal veto itself.
- `cutflow.py`: `get_signal_b_track_slots` (mode-dispatched daughter walk to
  TRK indices, reusing the Stage 3 slot map), `count_signal_b_tracks`
  (validation handle), `get_antibaryon_veto_mask`, `get_anti_lambda0_veto_mask`
  (add-on, measured not adopted), plus `get_pid_candidate_masks` /
  `get_pid_mask_per_B` (factored out of run_stage03's inline logic).
  The veto counts the event's opposite-charge selector-passing tracks once and
  subtracts each candidate's own (with explicit de-duplication), rather than
  materializing an O(nB x nTRK) per-B/per-track boolean.
- `pid_optimization.py`: `scan_antibaryon_veto` (reuses `evaluate_combo`
  unchanged) and `veto_efficiency_by_lambdac_mode`.
- `cumulative_performance.py` (new): `build_pipeline_masks`,
  `cumulative_cutflow`, `candidate_multiplicities`.
- `run_stage04.py --channel all`, `notebooks/04_antibaryon_veto.ipynb`, BAD
  subsection `sec:stage4veto` + two new generated tables
  (`generated_table_antibaryon_veto`,
  `generated_table_cumulative_cutflow_<channel>`).

**Structural facts verified empirically, not assumed** (both channels, full
MC files):

- `pSelectorsMap` is indexed by **track** (`len == nTRK`). The `p` hypothesis
  list is a strict subset of tracks and **24% of tracks passing
  SuperLooseKMProtonSelection lie outside it**, so the eligible-track pool is
  a real choice. Using the `p` list alone is **degenerate** -- it holds only
  ~2 candidates/event (essentially the signal protons, all of signal-baryon
  charge), so it contains almost no opposite-charge track and the veto
  approaches never firing. The per-track map is used, as in the reference.
- `sign(TRKLund)` is the electric charge: agrees with
  `sign(pLund[pTrkIdx])` in 207562/207563 proton entries (one isolated
  ntuple inconsistency, not a pattern).
- Signal baryon sign is consistent within a B candidate: both Lambda0s share
  a sign (Lam0Lam0) and `sign(LambdaCLund)` matches the Lambda0's proton
  (Lam0LamC, including mode 4's inner Lambda0) -- **0 mismatches** either
  way. So no separate same-baryon-number requirement is needed.
- Track walk validated: distinct resolvable tracks per B candidate are
  exactly 4 (Lam0Lam0) and 5/3/5/7 for LambdaC modes 1/2/3/4, with zero
  deviation over 368345 candidates. `run_stage04.py` asserts on this.

**Known limitation (Lam0LamC modes 2 and 3) -- flagged, low impact:** the
K_S0's two pion tracks **cannot be resolved** from these parquet files.
`K_Sd1Idx`/`K_Sd2Idx` index no collection present in the file: they exceed
`npi` for 94% of K_S0 candidates, and in TRK space the implied daughter
charges are opposite only 55% of the time (random = 50%). Where both indices
do land inside the pi collection (5.6%) the charges are opposite 100% of the
time and d1 is always positive, matching `K_Sd1Lund = +211` -- so the
intended index space *is* a pion list, just not the one stored here.
Kinematic matching to TRK is not a clean substitute (daughters are refit at
the displaced vertex; <1 MeV match only 45% of the time, consistent with
combinatorics). Those two pions therefore stay eligible to fire the veto in
modes 2/3. **Measured cost ~1%**: veto efficiency 0.900/0.888/0.891/0.893
for modes 1/2/3/4 -- mode 4 has no K_S0 and sits at 0.893, so the deficit is
small and partly confounded with track multiplicity anyway.

**Found while doing Stage 4 -- needs Matt's attention:** the `stage02`
multi-candidate numbers stored in `results/Lam0LamC.yaml` (16.3/75.4/8.3%)
are **stale** -- they predate the current `channel_config.py`. Re-running
`run_stage02.multi_candidate_numbers` against today's config gives
**23.2/70.0/6.8%** (verified bit-for-bit against Stage 4's own Stage-2 row,
distributions identical). The `stage03` numbers (32.4/62.6/5.0%) *are*
current and are reproduced exactly. Fix: re-run
`python run_stage02.py --channel all` and `generate_latex_macros.py`.

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
  Stage 4 additions: `scan_antibaryon_veto` (1D KM proton-ladder scan,
  reusing `evaluate_combo` unchanged) and `veto_efficiency_by_lambdac_mode`.
- `cumulative_performance.py` — Stage 4 Phase 3: `build_pipeline_masks`
  (the whole selection in pipeline order), `cumulative_cutflow` (signal/bkg
  in the signal box, data in fit-region-minus-signal-box), and
  `candidate_multiplicities` (Stage 2/3 convention: no single-candidate cut,
  no region cut, so it stays comparable across stages).
- `run_stage01.py`, `run_stage02.py`, `run_stage03.py`, `run_stage04.py` —
  write `results/<channel>.yaml`. `run_stage04.py` is the only one that
  reads collision data (fit region minus signal box only).
- `notebooks/01_load_and_diagnostics.ipynb` — full Stage 0/1 diagnostics,
  channel-parametrized, incl. blinding-verification cell.
- `notebooks/02_lambda_purity.ipynb` — Stage 2 scans/fits, cross-check
  plot, multi-candidate study, checkpoint cell. MC only.
- `notebooks/03_pid_optimization.ipynb` — Stage 3 KM-ladder PID scans
  (Lambda0 both channels; LambdaC per mode), no-PID-baseline sanity check,
  multi-candidate study with PID applied, checkpoint cell. MC only.
- `notebooks/04_antibaryon_veto.ipynb` — Stage 4 veto ladder scan,
  comparison points, per-mode efficiency, before/after diagnostics (mES,
  DeltaE, plus a Lambda0-mass null check), cumulative cutflow and
  multiplicity tables, checkpoint cell. **Written but never executed** —
  see the resume checklist at the end of this file.
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
   fraction from 8.3% to ~5.0% -- an improvement, not a resolution. Stage 4's
   veto takes it further, to ~4.4% (0/1/>1 = 39.7/55.9/4.4% at the
   recommended selector), still without resolving it -- and the zero-good-B
   fraction is now the dominant term, which is the real cost to weigh
   against relaxing the requirement. **This was deferred "until after the
   antibaryon-veto stage", which has now happened: it is ready to decide.**
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

5. **Stale `stage02` multi-candidate numbers** (found during Stage 4, see
   above): `results/Lam0LamC.yaml`'s stage02 section predates the current
   `channel_config.py`. Needs `run_stage02.py --channel all` re-run. Not a
   Stage 4 blocker -- Stage 4 recomputes its own Stage-2 row -- but the BAD
   quotes the stale numbers until it is fixed.

## Next steps (agreed)

1. Matt reviews the Stage 4 checkpoint: notebook 04, the BAD
   `sec:stage4veto` section, and `results/<channel>.yaml` section `stage04`.
   Decisions needed: (a) which veto selector to apply per channel, or
   whether to take the p-Lambda0 value for cross-analysis consistency;
   (b) confirm the anti-Lambda0 add-on stays off (measured: costs efficiency,
   no FOM gain); (c) whether to re-run Stage 2 to refresh the stale
   multi-candidate numbers; (d) the single-candidate requirement.
2. Stage 5: MLP classifier.

## Resume checklist (2026-07-26)

What is already done, so it is not repeated: `run_stage04.py --channel all`,
`generate_latex_macros.py`, and `latexmk` have all been run. Do NOT re-run
them unless something changes.

To finish the Stage 4 checkpoint:

```
# 1. Execute notebook 04 for BOTH channels (it has never been run --
#    no output cells, and its four plots do not exist yet).
cd BNV_2Lambda_analysis/notebooks
jupyter lab 04_antibaryon_veto.ipynb     # set CHANNEL, run all, for each channel

# 2. Refresh the stale stage02 numbers (independent of the above)
cd ..
python run_stage02.py --channel all

# 3. Only after 1 and/or 2:
python generate_latex_macros.py
./copy_plots_to_BAD.sh
cd ../BAD_2Lambda && latexmk -pdf main.tex
```

Then apply the accepted veto selector to `channel_config.py`'s
`antibaryon_veto['selector']` per channel (currently `None` = no veto), and
add explicit "applied"/"not applied" callouts to the BAD Stage 4 section,
matching the pattern used at the Stage 2 and Stage 3 checkpoints.

### Uncommitted working tree (nothing has been committed)

New files (untracked): `cumulative_performance.py`, `run_stage04.py`,
`notebooks/04_antibaryon_veto.ipynb`,
`BAD_2Lambda/generated_table_antibaryon_veto.tex`,
`BAD_2Lambda/generated_table_cumulative_cutflow_{Lam0Lam0,Lam0LamC}.tex`,
plus five Stage 2/3 diagnostic figures.

Modified: `cutflow.py`, `channel_config.py`, `pid_optimization.py`,
`generate_latex_macros.py`, `results/{Lam0Lam0,Lam0LamC}.yaml`,
`BAD_2Lambda/data_selection.tex`, `BAD_2Lambda/generated_numbers.tex`,
`WORKFLOW-README.md`, `STATUS.md`, and several Stage 2/3 figures.

The commit for this checkpoint has not been made -- ask before committing,
and never push (CLAUDE.md).

### Gotcha worth remembering

`BAD_2Lambda/data_selection.tex` was briefly broken by a careless global
`sed 's/LamLamZero/LamLam/g'`: the string `LamLamZero` occurs *inside*
pre-existing macro names (`\LamLamLamZeroPidFom` contains it at offset 3),
so four Stage 2/3 macros were silently mangled. Already fixed and verified
(zero undefined `LamLam*` macros across all `.tex`), but do not run
unanchored `sed` over that file again -- scope any replacement to the
specific macro names.
