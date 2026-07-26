#!/usr/bin/env python
"""
Stage 4: antibaryon (antiproton) veto, plus the cumulative Stage 1-4
selection-performance summary.

Physics: the BNV signal final state is all-baryon -- two protons, or two
antiprotons for the charge-conjugate B -- and contains no antibaryon, whereas
SM decays that produce a baryon also produce a compensating antibaryon. A B
candidate is vetoed if the event holds an identified antiproton that is not
one of that candidate's own tracks (cutflow.get_antibaryon_veto_mask). This
follows the B+ -> p Lambda0 reference (BNV_pLambda's
build_antiproton_antimask), generalized from its one-B-candidate-per-event
assumption to a per-candidate mask.

For each channel this:

1. Validates the signal-B track walk (cutflow.count_signal_b_tracks) against
   the expected final-state size before anything depends on it.
2. Scans the KM proton ladder for the veto with the Punzi FOM, the same
   machinery and same a=4 as Stage 3, and flags boundary-hugging optima.
   Note the ladder runs the opposite way from Stage 3: a LOOSER selector
   makes the veto fire MORE often.
3. Records the comparisons needed to judge the result: no veto at all, the
   p-Lambda0 reference selector, the reference-faithful exclusion scope, the
   p-list-only eligible-track pool, and (Lam0LamC) the per-mode efficiency
   that quantifies the unresolvable-K_S0-daughter caveat.
4. Builds the cumulative Stage 1-4 cutflow and per-step candidate
   multiplicities (cumulative_performance.py).

Steps 1-3 use MC only. Step 4 additionally reads the BLINDED collision file
and counts data candidates in the fit region with the signal box explicitly
removed -- the signal box is never read (see CLAUDE.md; blinding).

This does NOT modify channel_config.py: it prints and records the
*recommended* selector (results/<channel>.yaml, section 'stage04'). The
config's 'antibaryon_veto'/'selector' stays None -- i.e. no veto in force for
other code -- until a human accepts the Stage 4 checkpoint.

Usage:
    python run_stage04.py --channel Lam0Lam0
    python run_stage04.py --channel Lam0LamC
    python run_stage04.py --channel all
"""

import argparse

import awkward as ak
import numpy as np

from channel_config import get_channel_config, BACKGROUND_SP_MODES, SIGNAL_SP_MODE
import cumulative_performance as cp
import cutflow
import datasets
import pid_optimization as po
import pid_selector
import results_io

# The p-Lambda0 analysis's hand-picked choice (BNV_pLambda's cutflow), carried
# as an explicit comparison point rather than adopted by default.
REFERENCE_SELECTOR = 'TightKMProtonSelection'

# Expected number of RESOLVABLE distinct tracks per B candidate. Modes 2 and 3
# fall 2 short of their true final-state size (5 and 7) because the K_S0's
# pion tracks cannot be resolved from these files -- see the note in
# cutflow.get_signal_b_track_slots.
EXPECTED_TRACKS = {
    'Lam0Lam0': {None: 4},
    'Lam0LamC': {1: 5, 2: 3, 3: 5, 4: 7},
}

################################################################################
def validate_track_walk(data_sp, config):
    """
    Assert that the signal-B daughter walk yields the expected number of
    distinct tracks per candidate. The LambdaC daughter-slot map is
    mode-dependent and was established empirically; this turns it into a
    checked precondition rather than a standing assumption.
    """
    channel = config['name']
    n = cutflow.count_signal_b_tracks(data_sp, config)
    expected = EXPECTED_TRACKS[channel]

    out = {}
    if channel == 'Lam0Lam0':
        ok = bool(ak.all(n == expected[None]))
        out['all'] = {'expected': expected[None], 'ok': ok}
        if not ok:
            raise AssertionError(
                f"{channel}: signal-B track walk gave "
                f"{np.unique(ak.to_numpy(ak.flatten(n))).tolist()} distinct tracks/candidate, "
                f"expected {expected[None]}")
    else:
        mode = cutflow.get_lambdac_decay_mode_per_B(data_sp, config)
        for m, exp in expected.items():
            sel = (mode == m)
            if not bool(ak.any(sel)):
                continue
            ok = bool(ak.all(n[sel] == exp))
            out[int(m)] = {'expected': int(exp), 'ok': ok,
                           'k0s_daughters_unresolved': bool(m in (2, 3))}
            if not ok:
                raise AssertionError(
                    f"{channel} mode {m}: signal-B track walk gave "
                    f"{np.unique(ak.to_numpy(ak.flatten(n[sel]))).tolist()} distinct "
                    f"tracks/candidate, expected {exp}")

    print(f"  track-walk validation passed: {out}")
    return out

################################################################################
def summarize_veto_scan(df, best, ladder):
    """Summary dict (no full scan table) for the veto ladder scan."""
    row = best['best_row']
    return {
        'ladder': list(ladder),
        'recommended_selector': row['p_selector'],
        'fom_at_recommended': float(row['fom']),
        'sig_eff_at_recommended': float(row['sig_eff']),
        'bkg_weighted_at_recommended': float(row['bkg_weighted']),
        'bkg_raw_mc_at_recommended': int(row['bkg_raw_mc_candidates']),
        'at_scan_boundary': {k: bool(v) for k, v in best['at_scan_boundary'].items()},
        'any_at_scan_boundary': bool(best['any_at_scan_boundary']),
        'low_mc_stats_at_recommended': bool(best['low_mc_stats_at_best']),
        'per_rung': [
            {'selector': r['p_selector'], 'sig_eff': float(r['sig_eff']),
             'bkg_weighted': float(r['bkg_weighted']), 'fom': float(r['fom'])}
            for _, r in df.iterrows()
        ],
    }

################################################################################
def run_stage04(channel):
    print(f"\n========== Stage 4 for {channel} ==========")
    config = get_channel_config(channel)

    data_sp, _ = datasets.load_datasets(channel, sp_or_data='sp')
    datasets.add_derived_fields(data_sp, config)
    weights = datasets.get_scaling_weights(BACKGROUND_SP_MODES)

    out = {}

    # ------------------------------------------------------------------
    # 1. Validate the exclusion set before anything depends on it
    # ------------------------------------------------------------------
    out['track_walk_validation'] = validate_track_walk(data_sp, config)

    # ------------------------------------------------------------------
    # 2. Punzi scan over the KM proton ladder
    #
    # Baseline is the selection in force before the veto: Stage 2 purity AND
    # the Stage 3 PID cuts (now applied in channel_config.py), so the FOM
    # measures the veto's own marginal effect -- the same convention Stage 3
    # used when it baselined on Stage 2.
    # ------------------------------------------------------------------
    baseline_mask = (cutflow.get_composite_purity_masks_per_B(data_sp, config) &
                     cutflow.get_pid_mask_per_B(data_sp, config))

    ladder = pid_selector.KM_LADDER['p']
    df = po.scan_antibaryon_veto(data_sp, config, weights, baseline_mask, ladder)
    best = po.best_from_ladder_scan(df, ['p_selector_idx'])

    out['veto_scan'] = summarize_veto_scan(df, best, ladder)
    recommended = out['veto_scan']['recommended_selector']

    print(f"  veto: recommended {recommended} "
          f"(FOM {best['best_row']['fom']:.4f}, sig_eff {best['best_row']['sig_eff']:.4f}, "
          f"bkg {best['best_row']['bkg_weighted']:.2f}); "
          f"boundary={best['any_at_scan_boundary']}, "
          f"low_mc_stats={best['low_mc_stats_at_best']}")

    # ------------------------------------------------------------------
    # 3. Comparisons the checkpoint needs to judge the recommendation
    # ------------------------------------------------------------------
    no_veto = po.evaluate_combo(data_sp, config, weights, baseline_mask,
                                ak.ones_like(baseline_mask, dtype=bool))
    out['veto_scan']['fom_no_veto_at_all'] = float(no_veto['fom'])
    out['veto_scan']['bkg_weighted_no_veto_at_all'] = float(no_veto['bkg_weighted'])
    out['veto_scan']['recommended_beats_no_veto'] = bool(
        out['veto_scan']['fom_at_recommended'] > no_veto['fom'])

    ref_row = df[df['p_selector'] == REFERENCE_SELECTOR].iloc[0]
    out['veto_scan']['reference_selector'] = REFERENCE_SELECTOR
    out['veto_scan']['fom_at_reference_selector'] = float(ref_row['fom'])
    out['veto_scan']['sig_eff_at_reference_selector'] = float(ref_row['sig_eff'])
    out['veto_scan']['bkg_weighted_at_reference_selector'] = float(ref_row['bkg_weighted'])

    print(f"  no veto at all: FOM {no_veto['fom']:.4f} (bkg {no_veto['bkg_weighted']:.2f}); "
          f"p-Lambda0 reference {REFERENCE_SELECTOR}: FOM {ref_row['fom']:.4f}")

    # Alternative definitions, evaluated at the recommended selector so the
    # design choices are measured rather than asserted.
    notes = {
        'reference_scope_proton_daughters_only':
            ('p-Lambda0-reference exclusion scope: only the proton-hypothesis '
             'daughters are exempt, so the signal Lambda0 pions -- which carry '
             'exactly the antiproton charge -- can fake the veto.'),
        'pool_p_list_only':
            ('DEGENERATE, recorded as the justification for pool=all_tracks, not '
             'as a competing option: the "p" hypothesis collection holds only '
             '~2 candidates/event (essentially the signal protons themselves, '
             'all of signal-baryon charge), so it contains almost no '
             'opposite-charge track and the veto approaches never firing.'),
    }
    alt = {}
    for label, kwargs in (
            ('reference_scope_proton_daughters_only', {'scope': 'proton_daughters_only'}),
            ('pool_p_list_only', {'pool': 'p_list'})):
        mask = cutflow.get_antibaryon_veto_mask(data_sp, config, selector=recommended, **kwargs)
        r = po.evaluate_combo(data_sp, config, weights, baseline_mask, mask)
        alt[label] = {'sig_eff': float(r['sig_eff']), 'bkg_weighted': float(r['bkg_weighted']),
                      'fom': float(r['fom']), 'note': notes[label]}
        print(f"  alt [{label}]: sig_eff {r['sig_eff']:.6f}, FOM {r['fom']:.6f}")
    out['alternative_definitions'] = alt

    # Optional anti-Lambda0 add-on: measured, not adopted.
    anti_lam = cutflow.get_anti_lambda0_veto_mask(data_sp, config)
    combined = cutflow.get_antibaryon_veto_mask(data_sp, config, selector=recommended) & anti_lam
    r = po.evaluate_combo(data_sp, config, weights, baseline_mask, combined)
    out['anti_lambda0_addon'] = {
        'sig_eff': float(r['sig_eff']), 'bkg_weighted': float(r['bkg_weighted']),
        'fom': float(r['fom']),
        'note': ('measured, not adopted; expected near-redundant with the '
                 'single-candidate requirement'),
    }
    print(f"  anti-Lambda0 add-on: sig_eff {r['sig_eff']:.4f}, FOM {r['fom']:.4f}")

    # Per-mode efficiency: quantifies the unresolvable-K_S0-daughter caveat
    if 'LambdaC' in config['composites']:
        out['veto_efficiency_by_lambdac_mode'] = po.veto_efficiency_by_lambdac_mode(
            data_sp, config, baseline_mask, recommended)
        for m, v in out['veto_efficiency_by_lambdac_mode'].items():
            flag = ' (K_S0 daughters unresolved)' if v['k0s_daughters_unresolved'] else ''
            print(f"    mode {m}: veto eff {v['efficiency']:.4f}{flag}")

    # ------------------------------------------------------------------
    # 4. Cumulative Stage 1-4 performance (reads BLINDED collision data;
    #    the signal box is explicitly excluded -- see cumulative_performance)
    # ------------------------------------------------------------------
    _, data_col = datasets.load_datasets(channel, sp_or_data='col')
    datasets.add_derived_fields(data_col, config)

    df_cum = cp.cumulative_cutflow(data_sp, config, weights, data_collision=data_col,
                                   veto_selector=recommended)
    out['cumulative_cutflow'] = df_cum.to_dict(orient='records')

    print("\n  Cumulative cutflow (signal/bkg in signal box; data in fit region "
          "EXCLUDING signal box):")
    print(df_cum.to_string(index=False,
                           float_format=lambda v: f"{v:.4g}"))

    out['candidate_multiplicities'] = cp.candidate_multiplicities(
        data_sp, config, veto_selector=recommended)

    final = out['candidate_multiplicities'][max(out['candidate_multiplicities'])]
    print(f"\n  Endpoint good-B fractions (0 / 1 / >1): "
          f"{final['frac_zero_good_B']:.3f} / {final['frac_one_good_B']:.3f} / "
          f"{final['frac_multi_good_B']:.3f}")

    results_io.update_results(channel, 'stage04', out)

    print(f"\nRECOMMENDED channel_config.py update for {channel} (review at the "
          f"Stage 4 checkpoint before applying -- see notebooks/04_antibaryon_veto.ipynb):")
    print(f"  antibaryon_veto['selector'] = {recommended!r}")

    return out

################################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--channel', required=True,
                        choices=['Lam0Lam0', 'Lam0LamC', 'all'])
    args = parser.parse_args()

    channels = ['Lam0Lam0', 'Lam0LamC'] if args.channel == 'all' else [args.channel]
    for channel in channels:
        run_stage04(channel)
