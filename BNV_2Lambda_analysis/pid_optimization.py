"""
Stage 3 PID selector optimization (Punzi figure of merit) for the BNV
2-Lambda analyses.

Scans the KM-family selector ladder (pid_selector.KM_LADDER) for each PID
target (Lambda0's p/pi daughters; for Lam0LamC, additionally the LambdaC
daughters, which differ by decay mode -- see cutflow.get_lambdac_pid_mask),
optimizing the Punzi figure of merit (verified against
BNV_pLambda/PID_study_and_plots_for_BAD.ipynb, not guessed):

    fom = sig_eff / sqrt(bkg + a/2),   a = 4

sig_eff is the signal-MC efficiency of the PID cut *on top of* the current
Stage 1/2 selection (single-candidate + fit region + purity cuts, whatever
is presently in channel_config.py) -- i.e. this measures the PID cut's own
marginal effect, not efficiency from scratch. bkg is the luminosity-weighted
count of background MC (all of channel_config.BACKGROUND_SP_MODES, via
datasets.get_scaling_weights) surviving the same selection, directly in the
mES/DeltaE signal region -- no sideband subtraction, unlike Stage 2: this is
an MC-truth background estimate, not a data-driven one (Stage 5/6 will do
the data-driven sideband estimate). Both counts are per-B-candidate, summed
across all candidates in all events (no single-candidate requirement),
matching the reference notebook's convention.

The PID ladder is discrete (6 KM rungs, loosest -> tightest), unlike Stage
2's continuous scans, but the same "boundary-hugging optimum" concern
applies: is_at_scan_boundary from purity_optimization.py is reused directly
by scanning each ladder rung as an integer index. A second, PID-specific
failure mode is also checked: at a tight enough selector, the weighted
background estimate can collapse toward zero on very few raw MC candidates,
inflating the FOM from a noisy denominator rather than real background
rejection -- flagged per grid point via `low_mc_stats` (raw background MC
candidates in the signal region below `low_stats_floor`).
"""

import itertools

import awkward as ak
import numpy as np
import pandas as pd

import cutflow
import pid_selector
import purity_optimization as po
from channel_config import BACKGROUND_SP_MODES, SIGNAL_SP_MODE

PUNZI_A = 4.0

################################################################################
def punzi_fom(sig_eff, bkg_weighted, a=PUNZI_A):
    denom = bkg_weighted + a / 2.0
    return sig_eff / np.sqrt(denom) if denom > 0 else 0.0

################################################################################
def _pass_all(data, config, name):
    """All-True candidate mask shaped like composite `name`'s collection."""
    return ak.ones_like(data[config['composites'][name]['mass_var']], dtype=bool)

################################################################################
def _per_b_mask(data, config, candidate_overrides):
    """
    get_composite_purity_masks_per_B, but with only some composites' masks
    overridden (candidate_overrides); every other composite in this
    channel's config passes all its candidates (unconstrained), so the
    per-B AND only reflects the composite(s) actually being scanned.
    """
    candidate_masks = {name: _pass_all(data, config, name) for name in config['composites']}
    candidate_masks.update(candidate_overrides)
    return cutflow.get_composite_purity_masks_per_B(data, config, candidate_masks)

################################################################################
def evaluate_combo(data_sp, config, weights, baseline_mask, per_b_pid_mask, low_stats_floor=5):
    """
    Punzi FOM for one PID selector combination.

    baseline_mask: per-B mask (aligned with the B collection) with the
    current Stage 1/2 (purity, no PID) cuts applied, optionally further
    restricted (e.g. to one LambdaC decay mode) by the caller -- both the
    signal-efficiency numerator/denominator and the background count are
    computed on top of this, so the FOM measures the PID cut's marginal
    effect.
    per_b_pid_mask: per-B mask for the PID selector combination being
    evaluated (see cutflow.get_lambda0_pid_mask / get_lambdac_pid_mask +
    get_composite_purity_masks_per_B).
    """
    signal_box = cutflow.get_signal_region_mask(data_sp, config)
    spmode = data_sp['spmode']

    denom_mask = baseline_mask & signal_box & (spmode == SIGNAL_SP_MODE)
    numer_mask = denom_mask & per_b_pid_mask

    n_denom = float(ak.sum(denom_mask))
    n_numer = float(ak.sum(numer_mask))
    sig_eff = n_numer / n_denom if n_denom > 0 else 0.0

    bkg_weighted = 0.0
    bkg_raw = 0
    for spmode_val in BACKGROUND_SP_MODES:
        mask_bkg = baseline_mask & per_b_pid_mask & signal_box & (spmode == spmode_val)
        n_raw = int(ak.sum(mask_bkg))
        bkg_weighted += n_raw * weights[str(spmode_val)]
        bkg_raw += n_raw

    return {
        'sig_eff': sig_eff, 'n_sig_numer': int(n_numer), 'n_sig_denom': int(n_denom),
        'bkg_weighted': float(bkg_weighted), 'bkg_raw_mc_candidates': bkg_raw,
        'fom': punzi_fom(sig_eff, bkg_weighted),
        'low_mc_stats': bool(bkg_raw < low_stats_floor),
    }

################################################################################
def scan_lambda0_pid(data_sp, config, weights, baseline_mask,
                     p_ladder=None, pi_ladder=None):
    """
    2D grid scan over the Lambda0 proton/pion KM ladders (both channels).
    Reused as-is for Lam0LamC's LambdaC mode 4, whose inner Lambda0 shares
    this same selector choice (see scan_lambdac_mode_pid).
    """
    p_ladder = p_ladder or pid_selector.KM_LADDER['p']
    pi_ladder = pi_ladder or pid_selector.KM_LADDER['pi']

    rows = []
    for ip, p_sel in enumerate(p_ladder):
        for ipi, pi_sel in enumerate(pi_ladder):
            lambda0_pid = cutflow.get_lambda0_pid_mask(data_sp, p_selector=p_sel, pi_selector=pi_sel)
            per_b_pid = _per_b_mask(data_sp, config, {'Lambda0': lambda0_pid})

            result = evaluate_combo(data_sp, config, weights, baseline_mask, per_b_pid)
            result.update({'p_selector': p_sel, 'p_selector_idx': ip,
                           'pi_selector': pi_sel, 'pi_selector_idx': ipi})
            rows.append(result)

    return pd.DataFrame(rows)

################################################################################
def scan_lambdac_mode_pid(data_sp, config, weights, baseline_mask, mode, ladders,
                          lambda0_selector=None):
    """
    Grid scan over the KM ladder(s) relevant to one LambdaC decay mode
    (Lam0LamC only): `ladders` is {particle: ladder_list} for just the
    particles this mode needs (e.g. {'p':.., 'K':.., 'pi':..} for mode 1;
    {'p':..} for mode 2; {'p':.., 'pi':..} for mode 3; {'pi':..} for mode 4).

    Restricted to candidates whose LambdaC daughter is this decay mode
    (cutflow.get_lambdac_decay_mode_per_B), so each mode's ladder is
    optimized independently of the others.

    lambda0_selector: {'p': name_or_None, 'pi': name_or_None} -- the
    (already chosen, or provisional) Lambda0 selector, reused for mode 4's
    inner Lambda0 daughter (not re-scanned here; see scan_lambda0_pid).
    """
    lambda0_selector = lambda0_selector or {}
    lambda0_pid_for_inner = cutflow.get_lambda0_pid_mask(
        data_sp, p_selector=lambda0_selector.get('p'), pi_selector=lambda0_selector.get('pi'))

    mode_per_b = cutflow.get_lambdac_decay_mode_per_B(data_sp, config)
    mode_mask = (mode_per_b == mode)
    baseline_this_mode = baseline_mask & mode_mask

    particles = sorted(ladders.keys())
    grids = [ladders[p] for p in particles]

    rows = []
    for combo in itertools.product(*grids):
        sel_this_mode = dict(zip(particles, combo))
        lamc_pid = cutflow.get_lambdac_pid_mask(data_sp, {mode: sel_this_mode}, lambda0_pid_for_inner)
        per_b_pid = _per_b_mask(data_sp, config, {'LambdaC': lamc_pid})

        result = evaluate_combo(data_sp, config, weights, baseline_this_mode, per_b_pid)
        result['mode'] = mode
        for particle, sel_name in sel_this_mode.items():
            result[f'{particle}_selector'] = sel_name
            result[f'{particle}_selector_idx'] = ladders[particle].index(sel_name)
        rows.append(result)

    return pd.DataFrame(rows)

################################################################################
# Stage 4: antibaryon-veto selector scan
#
# Same Punzi FOM, same ladder machinery, same boundary/low-stats flagging as
# the Stage 3 PID scans above -- the veto is, mechanically, one more KM-ladder
# cut expressed as a per-B mask, so it reuses evaluate_combo unchanged and
# produces a DataFrame that best_from_ladder_scan and
# plotting.plot_pid_ladder_scan_1d consume without modification.
#
# The ladder runs the OPPOSITE way from Stage 3, which is worth keeping in
# mind when reading the scan: a LOOSER proton selector tags more tracks as
# antiprotons, so the veto fires more often -- more background rejection AND
# more signal loss. Tightening the selector moves toward "no veto at all".
################################################################################
def scan_antibaryon_veto(data_sp, config, weights, baseline_mask,
                         ladder=None, scope=None, pool=None):
    """
    1D scan over the KM proton ladder for the antibaryon veto.

    scope/pool default to the channel's 'antibaryon_veto' config; pass them
    explicitly to produce the alternative-definition comparisons reported at
    the checkpoint (reference-faithful scope, p-list-only pool).
    """
    veto_cfg = config.get('antibaryon_veto', {})
    scope = scope if scope is not None else veto_cfg.get('scope', 'all_signal_tracks')
    pool = pool if pool is not None else veto_cfg.get('pool', 'all_tracks')
    ladder = ladder or pid_selector.KM_LADDER['p']

    rows = []
    for i, sel in enumerate(ladder):
        veto_mask = cutflow.get_antibaryon_veto_mask(data_sp, config, selector=sel,
                                                     scope=scope, pool=pool)

        result = evaluate_combo(data_sp, config, weights, baseline_mask, veto_mask)
        result.update({'p_selector': sel, 'p_selector_idx': i,
                       'scope': scope, 'pool': pool})
        rows.append(result)

    return pd.DataFrame(rows)

################################################################################
def veto_efficiency_by_lambdac_mode(data_sp, config, baseline_mask, selector,
                                    scope=None, pool=None):
    """
    Signal-MC veto efficiency split by LambdaC decay mode (Lam0LamC only).

    Reported because the exclusion set is not equally complete across modes:
    modes 2 and 3 contain a K_S0 whose two pion tracks cannot be resolved
    from these files, so those pions stay eligible to fire the veto (see
    cutflow.get_signal_b_track_slots). Modes 1 and 4 have no K_S0 and are the
    clean comparison -- this function is what turns that caveat into a
    number in results/Lam0LamC.yaml instead of a footnote.
    """
    veto_cfg = config.get('antibaryon_veto', {})
    scope = scope if scope is not None else veto_cfg.get('scope', 'all_signal_tracks')
    pool = pool if pool is not None else veto_cfg.get('pool', 'all_tracks')

    signal_box = cutflow.get_signal_region_mask(data_sp, config)
    base = baseline_mask & signal_box & (data_sp['spmode'] == SIGNAL_SP_MODE)

    keep = cutflow.get_antibaryon_veto_mask(data_sp, config, selector=selector,
                                            scope=scope, pool=pool)
    mode = cutflow.get_lambdac_decay_mode_per_B(data_sp, config)

    out = {}
    for m in sorted(config['lambdac_modes']):
        b = base & (mode == m)
        n_den = int(ak.sum(b))
        n_num = int(ak.sum(b & keep))
        out[int(m)] = {
            'n_before': n_den, 'n_after': n_num,
            'efficiency': float(n_num / n_den) if n_den > 0 else 0.0,
            'k0s_daughters_unresolved': bool(m in (2, 3)),
        }
    return out

################################################################################
def best_from_ladder_scan(df, ladder_idx_cols, fom_col='fom'):
    """
    Best-FOM row from a ladder scan, plus a per-dimension boundary flag
    (reusing purity_optimization.is_at_scan_boundary on each ladder's
    integer rung-index column) and the row's own low_mc_stats flag.

    ladder_idx_cols: list of the '<particle>_selector_idx' column names
    scanned in this grid (one or more dimensions).
    """
    best = po.best_from_scan(df, fom_col)

    boundary_flags = {}
    for col in ladder_idx_cols:
        boundary_flags[col] = bool(po.is_at_scan_boundary(df, col, best))

    return {
        'best_row': best,
        'at_scan_boundary': boundary_flags,
        'any_at_scan_boundary': any(boundary_flags.values()),
        'low_mc_stats_at_best': bool(best.get('low_mc_stats', False)),
    }
################################################################################
