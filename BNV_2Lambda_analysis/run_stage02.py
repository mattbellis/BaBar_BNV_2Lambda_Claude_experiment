#!/usr/bin/env python
"""
Stage 2: Lambda0 (and, for Lam0LamC, K_S0 and LambdaC) purity optimization.

Uses MC only (signal + luminosity-weighted background) -- no collision data
is read, per the Stage 2 plan. For each channel this:

1. Re-optimizes the Lambda0 flight-significance cut and mass window with an
   S/sqrt(S+B) scan (sideband-subtracted signal MC vs. background MC).
2. (Lam0LamC only) Optimizes the K_S0 flight-significance cut and mass
   window the same way, FIRST -- before the LambdaC windows are set, since
   the K_S0 cut gates two of the four LambdaC decay modes (sequential
   optimization; see purity_optimization.py and cutflow.get_lambdac_k0s_gate).
3. (Lam0LamC only) Fits the per-mode LambdaC mass resolution (Gaussian +
   linear background) and derives mass_window_nsigma * sigma windows per
   mode, with the K_S0 gate applied.
4. (Lam0LamC only) Reports the multi-candidate numbers (candidates per
   event, fraction with >1) after all Stage 2 purity cuts -- numbers only,
   no candidate-selection policy yet (parked for a later stage).

This does NOT modify channel_config.py: it prints and records the
*recommended* cut values (results/<channel>.yaml, section 'stage02'); the
current channel_config.py values remain in force for other code until a
human reviews the Stage 2 checkpoint and applies them.

Usage:
    python run_stage02.py --channel Lam0Lam0
    python run_stage02.py --channel Lam0LamC
    python run_stage02.py --channel all
"""

import argparse

import awkward as ak
import numpy as np

from channel_config import get_channel_config, FOM_SIDEBAND_WIDTH_MULT, SIGNAL_SP_MODE, BACKGROUND_SP_MODES
import cutflow
import datasets
import purity_optimization as po
import results_io

################################################################################
def multi_candidate_numbers(data_sp, config):
    """
    Per-event count of B candidates whose composite daughters all pass
    their purity masks (config's *current* cut values -- caller sets these
    to whatever it wants "after all purity cuts" to mean), restricted to
    signal MC, plus the fraction of events with 0 / 1 / >1 such candidates.
    """
    mask_sig = data_sp['spmode'] == SIGNAL_SP_MODE

    _, candidate_masks = cutflow.get_all_composite_purity_masks(data_sp, config)
    mask_b = cutflow.get_composite_purity_masks_per_B(data_sp, config, candidate_masks)

    n_good_b = ak.to_numpy(ak.sum(mask_b[mask_sig], axis=1))
    n_events = len(n_good_b)

    vals, counts = np.unique(n_good_b, return_counts=True)

    return {
        'n_signal_events': int(n_events),
        'n_good_B_distribution': {int(v): int(c) for v, c in zip(vals, counts)},
        'frac_zero_good_B': float(np.sum(n_good_b == 0)) / n_events if n_events else 0.0,
        'frac_one_good_B': float(np.sum(n_good_b == 1)) / n_events if n_events else 0.0,
        'frac_multi_good_B': float(np.sum(n_good_b > 1)) / n_events if n_events else 0.0,
    }

################################################################################
def run_stage02(channel):
    print(f"\n========== Stage 2 for {channel} ==========")
    config = get_channel_config(channel)

    data_sp, _ = datasets.load_datasets(channel, sp_or_data='sp')
    datasets.add_derived_fields(data_sp, config)

    weights = datasets.get_scaling_weights(BACKGROUND_SP_MODES)

    out = {}

    # ------------------------------------------------------------------
    # Lambda0 purity optimization (both channels)
    # ------------------------------------------------------------------
    comp_l0 = config['composites']['Lambda0']
    sig_l0 = po.flat_signal_arrays(data_sp, [comp_l0['mass_var'], comp_l0['flight_var']])
    bkg_l0 = po.flat_background_arrays(data_sp, [comp_l0['mass_var'], comp_l0['flight_var']], weights)

    lambda0_opt = po.optimize_flight_and_mass(
        sig_l0, bkg_l0, comp_l0['mass_var'], comp_l0['flight_var'], comp_l0['mass_pdg'],
        comp_l0['mass_window'], comp_l0['flight_scan'], comp_l0['mass_halfwidth_scan'],
        FOM_SIDEBAND_WIDTH_MULT)

    out['lambda0'] = {
        'current_flight_cut': comp_l0['flight_cut'],
        'current_mass_window': [float(x) for x in comp_l0['mass_window']],
        'recommended_flight_cut': lambda0_opt['chosen_flight_cut'],
        'recommended_mass_window': lambda0_opt['chosen_mass_window'],
        'fom_at_recommended_flight_cut': lambda0_opt['fom_flight'],
        'fom_at_recommended_mass_window': lambda0_opt['fom_mass'],
        'flight_at_scan_boundary': lambda0_opt['flight_at_scan_boundary'],
        'mass_halfwidth_at_scan_boundary': lambda0_opt['mass_halfwidth_at_scan_boundary'],
    }
    print(f"Lambda0: recommended flight cut > {lambda0_opt['chosen_flight_cut']:.1f} "
          f"(FOM {lambda0_opt['fom_flight']:.2f}); "
          f"recommended mass window {lambda0_opt['chosen_mass_window']} "
          f"(FOM {lambda0_opt['fom_mass']:.2f})")

    recommended_config = get_channel_config(channel)
    recommended_config['composites']['Lambda0']['flight_cut'] = lambda0_opt['chosen_flight_cut']
    recommended_config['composites']['Lambda0']['mass_window'] = lambda0_opt['chosen_mass_window']

    if channel != 'Lam0LamC':
        out['multi_candidate_study'] = multi_candidate_numbers(data_sp, recommended_config)
        results_io.update_results(channel, 'stage02', out)
        return out

    # ------------------------------------------------------------------
    # K_S0 purity optimization (Lam0LamC only) -- FIRST, sequential order
    # ------------------------------------------------------------------
    k0s = config['k0s']
    sig_k0s = po.flat_signal_arrays(data_sp, [k0s['mass_var'], k0s['flight_var']])
    bkg_k0s = po.flat_background_arrays(data_sp, [k0s['mass_var'], k0s['flight_var']], weights)

    k0s_opt = po.optimize_flight_and_mass(
        sig_k0s, bkg_k0s, k0s['mass_var'], k0s['flight_var'], k0s['mass_pdg'],
        k0s['mass_window'], k0s['flight_scan'], k0s['mass_halfwidth_scan'],
        FOM_SIDEBAND_WIDTH_MULT)

    out['k0s'] = {
        'current_flight_cut': k0s['flight_cut'],
        'current_mass_window': [float(x) for x in k0s['mass_window']],
        'recommended_flight_cut': k0s_opt['chosen_flight_cut'],
        'recommended_mass_window': k0s_opt['chosen_mass_window'],
        'fom_at_recommended_flight_cut': k0s_opt['fom_flight'],
        'fom_at_recommended_mass_window': k0s_opt['fom_mass'],
        'flight_at_scan_boundary': k0s_opt['flight_at_scan_boundary'],
        'mass_halfwidth_at_scan_boundary': k0s_opt['mass_halfwidth_at_scan_boundary'],
    }
    print(f"K_S0: recommended flight cut > {k0s_opt['chosen_flight_cut']:.1f} "
          f"(FOM {k0s_opt['fom_flight']:.2f}); "
          f"recommended mass window {k0s_opt['chosen_mass_window']} "
          f"(FOM {k0s_opt['fom_mass']:.2f})")

    recommended_config['k0s']['flight_cut'] = k0s_opt['chosen_flight_cut']
    recommended_config['k0s']['mass_window'] = k0s_opt['chosen_mass_window']

    # ------------------------------------------------------------------
    # LambdaC per-mode mass-resolution fit, with the recommended K_S0 gate
    # applied (sequential order: K_S0 cuts are fixed above; this only
    # changes the *input sample* to the fit, not a scan)
    # ------------------------------------------------------------------
    mask_sig = data_sp['spmode'] == SIGNAL_SP_MODE

    mode = cutflow.get_lambdac_decay_mode(data_sp)
    mass_lamc = data_sp[config['composites']['LambdaC']['mass_var']]
    k0s_gate = cutflow.get_lambdac_k0s_gate(data_sp, recommended_config)

    mode_flat = ak.to_numpy(ak.flatten(mode[mask_sig], axis=None))
    mass_flat = ak.to_numpy(ak.flatten(mass_lamc[mask_sig], axis=None))
    gate_flat = ak.to_numpy(ak.flatten(k0s_gate[mask_sig], axis=None))

    lambdac_mass_pdg = config['composites']['LambdaC']['mass_pdg']
    nsigma = config['composites']['LambdaC']['mass_window_nsigma']

    mode_fit_summary = {}
    mode_fit_summary_gated = {}
    recommended_windows = {}

    for m in sorted(config['lambdac_modes'].keys()):
        sel_m = (mode_flat == m)
        fit = po.fit_mass_peak(mass_flat[sel_m], lambdac_mass_pdg)
        mode_fit_summary[m] = {'mu': fit['mu'], 'sigma': fit['sigma'],
                               'converged': fit['converged'], 'n_candidates': fit['n_candidates']}
        recommended_windows[m] = po.nsigma_window(fit['mu'], fit['sigma'], nsigma)

        # Cross-check: same fit with the K_S0 gate applied (a no-op for
        # modes 1 and 4, which have no K_S0 daughter) -- confirms the
        # LambdaC peak position/width are not strongly coupled to the K_S0
        # window choice.
        sel_m_gated = sel_m & gate_flat
        if sel_m_gated.sum() > 20:
            fit_gated = po.fit_mass_peak(mass_flat[sel_m_gated], lambdac_mass_pdg)
            mode_fit_summary_gated[m] = {'mu': fit_gated['mu'], 'sigma': fit_gated['sigma'],
                                         'converged': fit_gated['converged'],
                                         'n_candidates': fit_gated['n_candidates']}
        else:
            mode_fit_summary_gated[m] = None

        print(f"LambdaC mode {m}: mu={fit['mu']*1000:.2f} MeV sigma={fit['sigma']*1000:.2f} MeV "
              f"(converged={fit['converged']}, n={fit['n_candidates']}) "
              f"-> recommended window {[f'{x*1000:.1f}' for x in recommended_windows[m]]} MeV")

    out['lambdac_mode_fits'] = mode_fit_summary
    out['lambdac_mode_fits_k0s_gated'] = mode_fit_summary_gated
    out['lambdac_mass_window_nsigma'] = nsigma
    out['lambdac_recommended_windows_per_mode'] = recommended_windows

    recommended_config['composites']['LambdaC']['mass_windows_per_mode'] = recommended_windows

    # ------------------------------------------------------------------
    # Multi-candidate study, after ALL Stage 2 purity cuts (recommended
    # values). Numbers only -- no candidate-selection policy at this stage.
    # ------------------------------------------------------------------
    out['multi_candidate_study'] = multi_candidate_numbers(data_sp, recommended_config)
    mc = out['multi_candidate_study']
    print(f"Multi-candidate study: {mc['frac_zero_good_B']*100:.1f}% of signal events have 0 "
          f"good B candidates, {mc['frac_one_good_B']*100:.1f}% have 1, "
          f"{mc['frac_multi_good_B']*100:.1f}% have >1 (after all Stage 2 purity cuts).")

    results_io.update_results(channel, 'stage02', out)

    print("\nRECOMMENDED channel_config.py updates (review at the Stage 2 checkpoint "
          "before applying -- see notebooks/02_lambda_purity.ipynb):")
    print(f"  composites['Lambda0']['flight_cut']              = {lambda0_opt['chosen_flight_cut']}")
    print(f"  composites['Lambda0']['mass_window']              = {lambda0_opt['chosen_mass_window']}")
    print(f"  k0s['flight_cut']                                 = {k0s_opt['chosen_flight_cut']}")
    print(f"  k0s['mass_window']                                = {k0s_opt['chosen_mass_window']}")
    print(f"  composites['LambdaC']['mass_windows_per_mode']    = {recommended_windows}")

    return out

################################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--channel', required=True,
                        choices=['Lam0Lam0', 'Lam0LamC', 'all'])
    args = parser.parse_args()

    channels = ['Lam0Lam0', 'Lam0LamC'] if args.channel == 'all' else [args.channel]
    for channel in channels:
        run_stage02(channel)
