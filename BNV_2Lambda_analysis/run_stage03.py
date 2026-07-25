#!/usr/bin/env python
"""
Stage 3: PID selector optimization (Punzi figure of merit).

Uses MC only (signal + luminosity-weighted background) -- no collision data
is read, same as Stage 2. For each channel this:

1. Scans the KM-family selector ladder (pid_selector.KM_LADDER) for the
   Lambda0 proton/pion daughters, on top of the current Stage 1/2 selection
   (single-candidate + fit region + purity cuts, per channel_config.py).
2. (Lam0LamC only) Scans the KM ladder(s) relevant to each LambdaC decay
   mode's own daughters, independently per mode, with the Lambda0 scan's
   recommendation (step 1) already fixed and reused for mode 4's inner
   Lambda0 (sequential optimization, same spirit as Stage 2's K_S0-before-
   LambdaC ordering).
3. (Lam0LamC only) Redoes the multi-candidate study with the Stage 3
   recommended PID cuts folded in on top of Stage 2's purity cuts -- numbers
   only, no candidate-selection policy (see STATUS.md open decision 1).

This does NOT modify channel_config.py: it prints and records the
*recommended* selector choices (results/<channel>.yaml, section 'stage03');
the current channel_config.py 'pid' section (all None) remains in force for
other code until a human reviews the Stage 3 checkpoint and applies them.

Usage:
    python run_stage03.py --channel Lam0Lam0
    python run_stage03.py --channel Lam0LamC
    python run_stage03.py --channel all
"""

import argparse

import awkward as ak
import numpy as np

from channel_config import get_channel_config, BACKGROUND_SP_MODES, SIGNAL_SP_MODE
import cutflow
import datasets
import pid_optimization as po
import pid_selector
import results_io

################################################################################
def summarize_ladder_result(df, best, ladders):
    """Summary dict (no full scan table) for one ladder scan, for results_io."""
    particles = sorted(ladders.keys())
    return {
        'ladders': {p: list(ladders[p]) for p in particles},
        'recommended_selector': {p: best['best_row'][f'{p}_selector'] for p in particles},
        'fom_at_recommended': float(best['best_row']['fom']),
        'sig_eff_at_recommended': float(best['best_row']['sig_eff']),
        'bkg_weighted_at_recommended': float(best['best_row']['bkg_weighted']),
        'bkg_raw_mc_at_recommended': int(best['best_row']['bkg_raw_mc_candidates']),
        'at_scan_boundary': {k: bool(v) for k, v in best['at_scan_boundary'].items()},
        'any_at_scan_boundary': bool(best['any_at_scan_boundary']),
        'low_mc_stats_at_recommended': bool(best['low_mc_stats_at_best']),
    }

################################################################################
def multi_candidate_numbers_with_pid(data_sp, config, lambda0_selector, lambdac_selector_per_mode=None):
    """
    Stage 2's multi-candidate study (see run_stage02.multi_candidate_numbers),
    redone with the Stage 3 recommended PID selectors folded in on top of the
    Stage 2 purity cuts. Reports B, LambdaC, and Lambda0 (Lam0LamC: split by
    from-B vs from-LambdaC) candidate multiplicities, and (Lam0LamC only) the
    good-B / good-LambdaC breakdown by decay mode. Numbers only -- no
    candidate-selection policy (Phase 3 of the Stage 3 plan).
    """
    mask_sig = data_sp['spmode'] == SIGNAL_SP_MODE

    _, purity_masks = cutflow.get_all_composite_purity_masks(data_sp, config)

    lambda0_pid = cutflow.get_lambda0_pid_mask(
        data_sp, p_selector=lambda0_selector.get('p'), pi_selector=lambda0_selector.get('pi'))
    combined = dict(purity_masks)
    combined['Lambda0'] = purity_masks['Lambda0'] & lambda0_pid

    is_lam0lamc = 'LambdaC' in config['composites']
    if is_lam0lamc:
        lambdac_pid = cutflow.get_lambdac_pid_mask(data_sp, lambdac_selector_per_mode or {}, lambda0_pid)
        combined['LambdaC'] = purity_masks['LambdaC'] & lambdac_pid

    mask_b = cutflow.get_composite_purity_masks_per_B(data_sp, config, combined)
    n_good_b = ak.to_numpy(ak.sum(mask_b[mask_sig], axis=1))
    n_events = len(n_good_b)
    vals, counts = np.unique(n_good_b, return_counts=True)

    out = {
        'n_signal_events': int(n_events),
        'n_good_B_distribution': {int(v): int(c) for v, c in zip(vals, counts)},
        'frac_zero_good_B': float(np.sum(n_good_b == 0)) / n_events if n_events else 0.0,
        'frac_one_good_B': float(np.sum(n_good_b == 1)) / n_events if n_events else 0.0,
        'frac_multi_good_B': float(np.sum(n_good_b > 1)) / n_events if n_events else 0.0,
    }

    if not is_lam0lamc:
        return out

    # LambdaC candidate multiplicity, and good-LambdaC / good-B counts by mode
    n_good_lamc = ak.to_numpy(ak.sum(combined['LambdaC'][mask_sig], axis=1))
    out['n_good_LambdaC_distribution'] = {
        int(v): int(c) for v, c in zip(*np.unique(n_good_lamc, return_counts=True))}

    mode = cutflow.get_lambdac_decay_mode(data_sp)
    mode_per_b = cutflow.get_lambdac_decay_mode_per_B(data_sp, config)
    good_b_sig = mask_b & mask_sig
    good_lamc_sig = combined['LambdaC'] & mask_sig

    modes = sorted(config['lambdac_modes'].keys())
    out['n_good_B_by_lambdac_mode'] = {int(m): int(ak.sum(good_b_sig & (mode_per_b == m))) for m in modes}
    out['n_good_LambdaC_by_mode'] = {int(m): int(ak.sum(good_lamc_sig & (mode == m))) for m in modes}

    # Lambda0 candidate multiplicity, split by where it's used: the B's
    # direct daughter (from-B) vs. mode 4's LambdaC-embedded daughter
    # (from-LambdaC) -- reuses cutflow.indices_to_booleans, the same utility
    # get_lambdac_k0s_gate-style linkage functions build on.
    idxvar_b_lambda0 = dict(config['B_daughters'])['Lambda0']
    from_b_idx = data_sp[idxvar_b_lambda0]
    from_b_referenced = cutflow.indices_to_booleans(from_b_idx, combined['Lambda0'])
    n_good_lambda0_from_b = ak.to_numpy(ak.sum(combined['Lambda0'][mask_sig] & from_b_referenced[mask_sig], axis=1))
    out['n_good_Lambda0_from_B_distribution'] = {
        int(v): int(c) for v, c in zip(*np.unique(n_good_lambda0_from_b, return_counts=True))}

    from_lamc_idx = data_sp['LambdaCd1Idx'][mode == 4]
    from_lamc_referenced = cutflow.indices_to_booleans(from_lamc_idx, combined['Lambda0'])
    n_good_lambda0_from_lamc = ak.to_numpy(
        ak.sum(combined['Lambda0'][mask_sig] & from_lamc_referenced[mask_sig], axis=1))
    out['n_good_Lambda0_from_LambdaC_distribution'] = {
        int(v): int(c) for v, c in zip(*np.unique(n_good_lambda0_from_lamc, return_counts=True))}

    return out

################################################################################
def run_stage03(channel):
    print(f"\n========== Stage 3 for {channel} ==========")
    config = get_channel_config(channel)

    data_sp, _ = datasets.load_datasets(channel, sp_or_data='sp')
    datasets.add_derived_fields(data_sp, config)

    weights = datasets.get_scaling_weights(BACKGROUND_SP_MODES)

    # Stage 1/2 selection currently in force (purity cuts, no PID) -- both
    # the numerator/denominator of the Punzi FOM's signal efficiency, and
    # its background count, are computed on top of this.
    baseline_mask = cutflow.get_composite_purity_masks_per_B(data_sp, config)

    out = {}

    # ------------------------------------------------------------------
    # Lambda0 PID (both channels)
    # ------------------------------------------------------------------
    l0_ladders = {'p': pid_selector.KM_LADDER['p'], 'pi': pid_selector.KM_LADDER['pi']}
    df_l0 = po.scan_lambda0_pid(data_sp, config, weights, baseline_mask,
                               l0_ladders['p'], l0_ladders['pi'])
    best_l0 = po.best_from_ladder_scan(df_l0, ['p_selector_idx', 'pi_selector_idx'])

    out['lambda0_pid'] = summarize_ladder_result(df_l0, best_l0, l0_ladders)
    lambda0_selector = out['lambda0_pid']['recommended_selector']

    # Sanity check (persisted, not just printed): does the recommended
    # selector actually beat applying no PID cut at all, or does the true
    # optimum lie outside the scanned ladder entirely? Relevant whenever the
    # scan is boundary-hugging at the loosest rung (see STATUS.md).
    no_pid_mask = ak.ones_like(baseline_mask, dtype=bool)
    no_pid_result = po.evaluate_combo(data_sp, config, weights, baseline_mask, no_pid_mask)
    out['lambda0_pid']['fom_no_pid_at_all'] = float(no_pid_result['fom'])
    out['lambda0_pid']['bkg_weighted_no_pid_at_all'] = float(no_pid_result['bkg_weighted'])

    print(f"Lambda0 PID: recommended {lambda0_selector} "
          f"(FOM {best_l0['best_row']['fom']:.3f}, sig_eff {best_l0['best_row']['sig_eff']:.3f}, "
          f"bkg {best_l0['best_row']['bkg_weighted']:.2f}; no-PID-at-all FOM {no_pid_result['fom']:.3f}); "
          f"boundary={best_l0['any_at_scan_boundary']}, low_mc_stats={best_l0['low_mc_stats_at_best']}")

    if channel != 'Lam0LamC':
        out['multi_candidate_study'] = multi_candidate_numbers_with_pid(data_sp, config, lambda0_selector)
        results_io.update_results(channel, 'stage03', out)
        return out

    # ------------------------------------------------------------------
    # LambdaC PID, per decay mode (Lam0LamC only) -- Lambda0 selector above
    # is fixed first; mode 4's inner Lambda0 reuses it (sequential order).
    # ------------------------------------------------------------------
    mode_ladders = {
        1: {'p': pid_selector.KM_LADDER['p'], 'K': pid_selector.KM_LADDER['K'], 'pi': pid_selector.KM_LADDER['pi']},
        2: {'p': pid_selector.KM_LADDER['p']},
        3: {'p': pid_selector.KM_LADDER['p'], 'pi': pid_selector.KM_LADDER['pi']},
        4: {'pi': pid_selector.KM_LADDER['pi']},
    }

    lambdac_selector_per_mode = {}
    out['lambdac_pid_per_mode'] = {}
    for m, ladders in mode_ladders.items():
        df_m = po.scan_lambdac_mode_pid(data_sp, config, weights, baseline_mask, m, ladders,
                                        lambda0_selector=lambda0_selector)
        idx_cols = [f'{p}_selector_idx' for p in ladders]
        best_m = po.best_from_ladder_scan(df_m, idx_cols)

        summary = summarize_ladder_result(df_m, best_m, ladders)
        out['lambdac_pid_per_mode'][m] = summary
        lambdac_selector_per_mode[m] = summary['recommended_selector']

        print(f"LambdaC mode {m} PID: recommended {summary['recommended_selector']} "
              f"(FOM {summary['fom_at_recommended']:.3f}, sig_eff {summary['sig_eff_at_recommended']:.3f}, "
              f"bkg {summary['bkg_weighted_at_recommended']:.2f}); "
              f"boundary={summary['any_at_scan_boundary']}, low_mc_stats={summary['low_mc_stats_at_recommended']}")

    out['multi_candidate_study'] = multi_candidate_numbers_with_pid(
        data_sp, config, lambda0_selector, lambdac_selector_per_mode)

    results_io.update_results(channel, 'stage03', out)

    print("\nRECOMMENDED channel_config.py updates (review at the Stage 3 checkpoint "
          "before applying -- see notebooks/03_pid_optimization.ipynb):")
    print(f"  pid['lambda0_selector']           = {lambda0_selector}")
    print(f"  pid['lambdac_selector_per_mode']  = {lambdac_selector_per_mode}")

    return out

################################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--channel', required=True,
                        choices=['Lam0Lam0', 'Lam0LamC', 'all'])
    args = parser.parse_args()

    channels = ['Lam0Lam0', 'Lam0LamC'] if args.channel == 'all' else [args.channel]
    for channel in channels:
        run_stage03(channel)
