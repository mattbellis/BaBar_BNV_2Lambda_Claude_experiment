"""
Cumulative selection performance through Stages 1-4 (Phase 3 of Stage 4).

Puts the whole chain of cuts currently in force -- preselection /
single-candidate, mES/DeltaE fit region, Stage 2 purity, Stage 3 PID, Stage 4
antibaryon veto -- into one table per channel, so the compounded effect is
visible in one place rather than spread across four stage checkpoints.

Counting conventions (stated explicitly because they are easy to misread):

- Signal MC and background MC are counted as B CANDIDATES inside the
  mES/DeltaE SIGNAL region, matching the Stage 3 Punzi FOM convention
  (pid_optimization.evaluate_combo). Because the signal region is a subset of
  the fit region, the "+ fit region" step is a no-op for these two columns by
  construction -- it is kept in the table so the pipeline order is complete
  and so its effect on the data column is visible.
- Background MC is luminosity-weighted (datasets.get_scaling_weights over
  channel_config.BACKGROUND_SP_MODES).
- COLLISION DATA is counted in the fit region with the signal box REMOVED
  (get_fit_mask & ~get_signal_region_mask). BLINDING: the signal box is never
  read. The upstream _BLINDED files already have it removed, but the mask is
  applied explicitly here as well so the column means what its label says and
  cannot be silently reinterpreted later.

Event-level counts are reported alongside candidate counts because the
single-candidate requirement is an event cut while everything after it is a
candidate cut; quoting only one of the two makes the first row misleading.
"""

import awkward as ak
import numpy as np
import pandas as pd

import cutflow
from channel_config import BACKGROUND_SP_MODES, SIGNAL_SP_MODE, DATA_SP_MODE

################################################################################
def build_pipeline_masks(data, config, veto_selector=None, include_veto=True):
    """
    The cumulative selection, in pipeline order.

    Returns a list of dicts, each {'name', 'per_b', 'event'}, where 'per_b'
    is the cumulative per-B-candidate mask up to and including that step and
    'event' is the cumulative event-level mask.

    veto_selector: the Stage 4 selector to use for the final step. Defaults
    to the channel config's 'antibaryon_veto'/'selector', which is None (=
    veto not applied) until the Stage 4 checkpoint is accepted -- so pass the
    scan's recommendation explicitly to see the veto's effect before it is
    adopted.
    """
    veto_cfg = config.get('antibaryon_veto', {}) or {}
    if veto_selector is None:
        veto_selector = veto_cfg.get('selector')

    steps = []

    # 1. Preselection / single-candidate requirement (an EVENT cut)
    event_mask = cutflow.get_single_candidate_mask(data, config)
    per_b = ak.ones_like(data['BpostFitMes'], dtype=bool) & event_mask
    steps.append({'name': 'preselection + single candidate',
                  'per_b': per_b, 'event': event_mask})

    # 2. mES/DeltaE fit region
    per_b = per_b & cutflow.get_fit_mask(data, config)
    steps.append({'name': '+ fit region',
                  'per_b': per_b, 'event': ak.any(per_b, axis=-1)})

    # 3. Stage 2 purity (mass windows, flight significance, K_S0 gate)
    per_b = per_b & cutflow.get_composite_purity_masks_per_B(data, config)
    steps.append({'name': '+ Stage 2 purity',
                  'per_b': per_b, 'event': ak.any(per_b, axis=-1)})

    # 4. Stage 3 PID selectors
    per_b = per_b & cutflow.get_pid_mask_per_B(data, config)
    steps.append({'name': '+ Stage 3 PID',
                  'per_b': per_b, 'event': ak.any(per_b, axis=-1)})

    # 5. Stage 4 antibaryon veto
    if include_veto:
        per_b = per_b & cutflow.get_antibaryon_veto_mask(
            data, config, selector=veto_selector,
            scope=veto_cfg.get('scope', 'all_signal_tracks'),
            pool=veto_cfg.get('pool', 'all_tracks'))
        label = ('+ Stage 4 antibaryon veto' if veto_selector is not None
                 else '+ Stage 4 antibaryon veto (NOT APPLIED: selector=None)')
        steps.append({'name': label, 'per_b': per_b, 'event': ak.any(per_b, axis=-1)})

    return steps

################################################################################
def cumulative_cutflow(data_sp, config, weights, data_collision=None,
                       veto_selector=None):
    """
    Cumulative cutflow table (one row per pipeline step) as a DataFrame.

    Columns: signal-MC candidates/events in the signal box with efficiency
    relative to the previous step and to the Stage 1 baseline (row 0);
    luminosity-weighted background-MC yield in the signal box with rejection
    relative to the baseline; and collision-data candidates in the fit region
    EXCLUDING the signal box (see the module docstring on blinding).

    data_collision may be None (MC-only run); the data columns are then NaN.
    """
    steps = build_pipeline_masks(data_sp, config, veto_selector=veto_selector)

    signal_box = cutflow.get_signal_region_mask(data_sp, config)
    spmode = data_sp['spmode']
    is_sig = (spmode == SIGNAL_SP_MODE)

    data_steps, data_region = None, None
    if data_collision is not None:
        data_steps = build_pipeline_masks(data_collision, config, veto_selector=veto_selector)
        # Fit region MINUS the signal box -- the blinded region is never counted.
        data_region = (cutflow.get_fit_mask(data_collision, config) &
                       ~cutflow.get_signal_region_mask(data_collision, config))

    rows = []
    n_sig_base = n_bkg_base = None

    for i, step in enumerate(steps):
        sel = step['per_b'] & signal_box

        n_sig = int(ak.sum(sel & is_sig))
        n_sig_events = int(ak.sum(ak.any(sel & is_sig, axis=-1)))

        bkg_weighted, bkg_raw = 0.0, 0
        for sp in BACKGROUND_SP_MODES:
            n_raw = int(ak.sum(sel & (spmode == sp)))
            bkg_weighted += n_raw * weights[str(sp)]
            bkg_raw += n_raw

        if i == 0:
            n_sig_base, n_bkg_base = n_sig, bkg_weighted

        row = {
            'step': i,
            'name': step['name'],
            'n_signal_candidates': n_sig,
            'n_signal_events': n_sig_events,
            'sig_eff_rel_prev': (n_sig / rows[-1]['n_signal_candidates']
                                 if i > 0 and rows[-1]['n_signal_candidates'] > 0 else 1.0),
            'sig_eff_rel_baseline': (n_sig / n_sig_base) if n_sig_base else 0.0,
            'bkg_weighted': float(bkg_weighted),
            'bkg_raw_mc_candidates': bkg_raw,
            'bkg_rejection_rel_baseline': (1.0 - bkg_weighted / n_bkg_base) if n_bkg_base else 0.0,
        }

        if data_steps is not None:
            dsel = data_steps[i]['per_b'] & data_region
            row['n_data_candidates_fit_excl_signal'] = int(ak.sum(dsel))
            row['n_data_events_fit_excl_signal'] = int(ak.sum(ak.any(dsel, axis=-1)))
        else:
            row['n_data_candidates_fit_excl_signal'] = np.nan
            row['n_data_events_fit_excl_signal'] = np.nan

        rows.append(row)

    return pd.DataFrame(rows)

################################################################################
def _distribution(counts):
    vals, n = np.unique(counts, return_counts=True)
    return {int(v): int(c) for v, c in zip(vals, n)}

################################################################################
def candidate_multiplicities(data_sp, config, veto_selector=None):
    """
    Per-event surviving-candidate multiplicities on signal MC, after Stage 2
    purity, after +Stage 3 PID, and after +Stage 4 antibaryon veto.

    IMPORTANT -- this deliberately does NOT use build_pipeline_masks. It
    reproduces the convention of Stage 2's and Stage 3's multi-candidate
    studies (run_stage02.multi_candidate_numbers /
    run_stage03.multi_candidate_numbers_with_pid): no single-candidate
    requirement and no mES/DeltaE region cut. Both exclusions matter:
    applying the pipeline's nB == 1 cut would force the multi-candidate
    fraction to 0 by construction and make the whole study vacuous, which is
    precisely the question the parked single-candidate decision (STATUS.md
    open decision 1) is asking about. The endpoint is therefore directly
    comparable to the Stage 2 (16.3/75.4/8.3%) and Stage 3 (32.4/62.6/5.0%)
    zero/one/multi good-B fractions.

    Numbers only -- no candidate-selection policy is applied or recommended.
    """
    mask_sig = data_sp['spmode'] == SIGNAL_SP_MODE
    n_events = int(ak.sum(mask_sig))

    veto_cfg = config.get('antibaryon_veto', {}) or {}

    # Candidate-level masks, accumulated in stage order (same construction as
    # run_stage03.multi_candidate_numbers_with_pid).
    _, purity_masks = cutflow.get_all_composite_purity_masks(data_sp, config)
    pid_masks = cutflow.get_pid_candidate_masks(data_sp, config)

    with_pid = {name: purity_masks[name] & pid_masks[name] for name in purity_masks}

    veto_per_b = cutflow.get_antibaryon_veto_mask(
        data_sp, config, selector=veto_selector,
        scope=veto_cfg.get('scope', 'all_signal_tracks'),
        pool=veto_cfg.get('pool', 'all_tracks'))

    stages = [
        ('Stage 2 purity', purity_masks, None),
        ('+ Stage 3 PID', with_pid, None),
        ('+ Stage 4 antibaryon veto', with_pid, veto_per_b),
    ]

    is_lam0lamc = 'LambdaC' in config['composites']
    if is_lam0lamc:
        mode = cutflow.get_lambdac_decay_mode(data_sp)
        mode_per_b = cutflow.get_lambdac_decay_mode_per_B(data_sp, config)
    idxvar_lam0_from_b = dict(config['B_daughters'])['Lambda0']

    out = {}
    for i, (name, cand_masks, extra_per_b) in enumerate(stages):
        mask_b = cutflow.get_composite_purity_masks_per_B(data_sp, config, cand_masks)
        if extra_per_b is not None:
            mask_b = mask_b & extra_per_b

        n_good_b = ak.to_numpy(ak.sum(mask_b[mask_sig], axis=1))

        entry = {
            'step': i,
            'name': name,
            'n_signal_events': n_events,
            'n_good_B_distribution': _distribution(n_good_b),
            'frac_zero_good_B': float(np.mean(n_good_b == 0)) if n_events else 0.0,
            'frac_one_good_B': float(np.mean(n_good_b == 1)) if n_events else 0.0,
            'frac_multi_good_B': float(np.mean(n_good_b > 1)) if n_events else 0.0,
        }

        # Lambda0 candidates, split by the role they play: the B's direct
        # daughter vs. the one embedded in a mode-4 LambdaC.
        lam0_mask = cand_masks['Lambda0']
        from_b_ref = cutflow.indices_to_booleans(data_sp[idxvar_lam0_from_b], lam0_mask)
        entry['n_good_Lambda0_from_B_distribution'] = _distribution(
            ak.to_numpy(ak.sum(lam0_mask[mask_sig] & from_b_ref[mask_sig], axis=1)))

        if is_lam0lamc:
            lamc_mask = cand_masks['LambdaC']
            entry['n_good_LambdaC_distribution'] = _distribution(
                ak.to_numpy(ak.sum(lamc_mask[mask_sig], axis=1)))

            from_lamc_idx = data_sp['LambdaCd1Idx'][mode == 4]
            from_lamc_ref = cutflow.indices_to_booleans(from_lamc_idx, lam0_mask)
            entry['n_good_Lambda0_from_LambdaC_distribution'] = _distribution(
                ak.to_numpy(ak.sum(lam0_mask[mask_sig] & from_lamc_ref[mask_sig], axis=1)))

            good_b_sig = mask_b & mask_sig
            entry['n_good_B_by_lambdac_mode'] = {
                int(m): int(ak.sum(good_b_sig & (mode_per_b == m)))
                for m in sorted(config['lambdac_modes'])}
            entry['n_good_LambdaC_by_mode'] = {
                int(m): int(ak.sum(lamc_mask & mask_sig & (mode == m)))
                for m in sorted(config['lambdac_modes'])}

        out[i] = entry

    return out
