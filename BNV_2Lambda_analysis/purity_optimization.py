"""
Stage 2 purity-cut optimization for the BNV 2-Lambda analyses.

Two independent things live here:

1. S/sqrt(S+B) scans (sideband-subtracted signal MC vs. luminosity-weighted
   background MC) used to choose the Lambda0 and K_S0 mass-window and
   flight-significance cuts. The peak/sideband convention mirrors the
   B+ -> p Lambda0 reference study (BNV_pLambda/Lambda_purity_studies.ipynb):
   sideband bands are contiguous with the peak window and have the same
   width as the peak (channel_config.FOM_SIDEBAND_WIDTH_MULT = 1.0).

2. A Gaussian-plus-linear-background fit used to extract the per-mode
   LambdaC mass resolution, which sets the (mode-dependent) mass windows --
   those are NOT chosen from an S/sqrt(S+B) scan, per the Stage 2 plan.

All scan/fit functions take FLAT arrays (one entry per candidate); flattening
and any selection that should already be applied (e.g. a chosen flight cut,
when scanning the mass window) is the caller's responsibility -- see
run_stage02.py.
"""

import numpy as np
import pandas as pd
import awkward as ak
from scipy.optimize import curve_fit

from channel_config import FOM_SIDEBAND_WIDTH_MULT, BACKGROUND_SP_MODES, SIGNAL_SP_MODE

################################################################################
def _to_numpy(x):
    if x is None:
        return None
    return ak.to_numpy(x) if isinstance(x, ak.Array) else np.asarray(x)

################################################################################
def sideband_subtracted_counts(mass, weight, lo, hi, sideband_mult=FOM_SIDEBAND_WIDTH_MULT):
    """
    (S, B) for one mass-peak window [lo, hi]:

        B = average of the two sideband bands (each of width
            sideband_mult * (hi - lo), contiguous with the peak edges)
        S = weighted peak count - B

    weight=None means unit weight (used for signal MC, which is filled with
    weight 1 elsewhere too).
    """
    mass = _to_numpy(mass)
    weight = np.ones_like(mass) if weight is None else _to_numpy(weight)

    width = hi - lo
    sb_width = sideband_mult * width

    peak = (mass > lo) & (mass < hi)
    sb_lo = (mass > lo - sb_width) & (mass <= lo)
    sb_hi = (mass >= hi) & (mass < hi + sb_width)

    npeak = float(weight[peak].sum())
    nsblo = float(weight[sb_lo].sum())
    nsbhi = float(weight[sb_hi].sum())

    B = 0.5 * (nsblo + nsbhi)
    S = npeak - B
    return S, B

################################################################################
def fom_s_over_sqrt_s_plus_b(S, B):
    return S / np.sqrt(S + B) if (S + B) > 0 else 0.0

################################################################################
def scan_threshold_cut(cutvar_sig, mass_sig, cutvar_bkg, mass_bkg, weight_bkg,
                       cut_values, mass_lo, mass_hi, sideband_mult=FOM_SIDEBAND_WIDTH_MULT):
    """
    Scan a 'greater than' threshold on cutvar (e.g. flight significance) at
    a fixed mass-peak window, computing S/sqrt(S+B) at each cut value.
    """
    cutvar_sig = _to_numpy(cutvar_sig)
    mass_sig = _to_numpy(mass_sig)
    cutvar_bkg = _to_numpy(cutvar_bkg)
    mass_bkg = _to_numpy(mass_bkg)
    weight_bkg = _to_numpy(weight_bkg)

    rows = []
    for cut in cut_values:
        sel_sig = cutvar_sig > cut
        sel_bkg = cutvar_bkg > cut

        S, _ = sideband_subtracted_counts(mass_sig[sel_sig], None, mass_lo, mass_hi, sideband_mult)
        _, B = sideband_subtracted_counts(mass_bkg[sel_bkg], weight_bkg[sel_bkg], mass_lo, mass_hi, sideband_mult)

        rows.append({
            'cut': float(cut), 'S': S, 'B': B,
            'fom': fom_s_over_sqrt_s_plus_b(S, B),
            'sig_eff': float(sel_sig.sum()) / len(cutvar_sig) if len(cutvar_sig) else 0.0,
        })

    return pd.DataFrame(rows)

################################################################################
def scan_mass_halfwidth(mass_sig, mass_bkg, weight_bkg, halfwidths, mass_center,
                        sideband_mult=FOM_SIDEBAND_WIDTH_MULT):
    """
    Scan the mass-window half-width around mass_center (any flight cut
    already applied by the caller), computing S/sqrt(S+B) at each width.
    """
    mass_sig = _to_numpy(mass_sig)
    mass_bkg = _to_numpy(mass_bkg)
    weight_bkg = _to_numpy(weight_bkg)

    rows = []
    for hw in halfwidths:
        lo, hi = mass_center - hw, mass_center + hw
        S, _ = sideband_subtracted_counts(mass_sig, None, lo, hi, sideband_mult)
        _, B = sideband_subtracted_counts(mass_bkg, weight_bkg, lo, hi, sideband_mult)
        rows.append({'halfwidth': float(hw), 'S': S, 'B': B,
                     'fom': fom_s_over_sqrt_s_plus_b(S, B)})

    return pd.DataFrame(rows)

################################################################################
def best_from_scan(df, fom_col='fom'):
    """Row (as a dict) with the maximum FOM in a scan DataFrame."""
    return df.loc[df[fom_col].idxmax()].to_dict()

################################################################################
def is_at_scan_boundary(df, xcol, best_row, edge_frac=0.02):
    """
    True if best_row[xcol] is (within edge_frac of the range) the first or
    last value scanned in df[xcol] -- a sign the true optimum may lie
    outside the scanned range (or that the FOM is monotonic across it and
    the 'optimum' is not a real interior maximum). Callers should surface
    this rather than silently trust a boundary-hugging recommendation.
    """
    xmin, xmax = df[xcol].min(), df[xcol].max()
    span = xmax - xmin
    if span <= 0:
        return False
    x = best_row[xcol]
    return (x - xmin) <= edge_frac * span or (xmax - x) <= edge_frac * span

################################################################################
# LambdaC per-mode mass-resolution fit
################################################################################
def _gauss_plus_linear(x, amp, mu, sigma, a, b):
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + a + b * x

################################################################################
def fit_mass_peak(mass_values, center_guess, sigma_guess=0.01, fit_halfrange=0.05, nbins=60):
    """
    Bin `mass_values` (flat array, e.g. one LambdaC decay mode's candidates)
    and fit a Gaussian-plus-linear-background model in
    [center_guess - fit_halfrange, center_guess + fit_halfrange].

    Returns a dict with the fit parameters plus the binned data and fit
    curve (for plotting). Falls back to a clipped mean/std of the windowed
    data (converged=False) if the fit does not converge -- the notebook
    should flag this for a mode rather than silently trust it.
    """
    mass_values = _to_numpy(mass_values)

    lo, hi = center_guess - fit_halfrange, center_guess + fit_halfrange
    x = mass_values[(mass_values > lo) & (mass_values < hi)]

    counts, edges = np.histogram(x, bins=nbins, range=(lo, hi))
    centers = 0.5 * (edges[:-1] + edges[1:])

    edge_bkg_guess = float(np.mean(np.concatenate([counts[:5], counts[-5:]]))) if nbins >= 10 else 0.0
    p0 = [max(counts.max() - edge_bkg_guess, 1.0), center_guess, sigma_guess, edge_bkg_guess, 0.0]

    try:
        popt, pcov = curve_fit(
            _gauss_plus_linear, centers, counts, p0=p0,
            bounds=([0, lo, 1e-4, 0, -np.inf], [np.inf, hi, fit_halfrange, np.inf, np.inf]),
            maxfev=10000)
        perr = np.sqrt(np.diag(pcov))
        result = {'mu': float(popt[1]), 'sigma': float(abs(popt[2])),
                  'mu_err': float(perr[1]), 'sigma_err': float(perr[2]),
                  'amp': float(popt[0]), 'a': float(popt[3]), 'b': float(popt[4]),
                  'converged': True}
    except Exception as exc:
        print(f"WARNING: Gaussian fit failed ({exc}); falling back to a "
              f"windowed mean/std estimate")
        result = {'mu': float(np.mean(x)) if len(x) else center_guess,
                  'sigma': float(np.std(x)) if len(x) else sigma_guess,
                  'mu_err': None, 'sigma_err': None, 'amp': None, 'a': None, 'b': None,
                  'converged': False}

    result['n_candidates'] = int(len(x))
    result['bin_centers'] = centers
    result['bin_counts'] = counts
    result['fit_curve'] = (_gauss_plus_linear(centers, result['amp'], result['mu'],
                                              result['sigma'], result['a'], result['b'])
                           if result['converged'] else None)

    return result

################################################################################
def nsigma_window(mu, sigma, nsigma):
    return [float(mu - nsigma * sigma), float(mu + nsigma * sigma)]

################################################################################
# Data-array helpers (candidate-level, flattened, aligned across variables)
################################################################################
def flat_signal_arrays(data_sp, varnames):
    """Flattened candidate-level arrays for signal MC, aligned across varnames."""
    sub = data_sp[data_sp['spmode'] == SIGNAL_SP_MODE]
    return {v: ak.to_numpy(ak.flatten(sub[v], axis=None)) for v in varnames}

################################################################################
def flat_background_arrays(data_sp, varnames, weights):
    """
    Flattened, luminosity-weighted candidate-level arrays for all background
    SP modes present in data_sp, aligned across varnames, plus a 'weight'
    array (one entry per candidate, repeating that mode's scaling weight).
    """
    per_var = {v: [] for v in varnames}
    wts = []

    for spmode in BACKGROUND_SP_MODES:
        sub = data_sp[data_sp['spmode'] == spmode]
        if len(sub) == 0:
            continue
        n = None
        for v in varnames:
            flat = ak.to_numpy(ak.flatten(sub[v], axis=None))
            per_var[v].append(flat)
            n = len(flat)
        wts.append(np.full(n, weights[str(spmode)]))

    out = {v: (np.concatenate(a) if a else np.array([])) for v, a in per_var.items()}
    out['weight'] = np.concatenate(wts) if wts else np.array([])
    return out

################################################################################
# Sequential 2-step optimization: flight cut first (nominal mass window),
# then the mass half-width (chosen flight cut applied)
################################################################################
def optimize_flight_and_mass(sig, bkg, mass_var, flight_var, mass_pdg, nominal_window,
                             flight_scan_cfg, mass_scan_cfg, sideband_mult=FOM_SIDEBAND_WIDTH_MULT):
    """
    sig, bkg: dicts from flat_signal_arrays / flat_background_arrays,
    already containing mass_var and flight_var (and 'weight' for bkg).

    Returns a dict with both scan DataFrames and the recommended
    (flight_cut, mass_window, FOM at each).
    """
    lo0, hi0 = nominal_window
    flight_cuts = np.arange(flight_scan_cfg['lo'], flight_scan_cfg['hi'], flight_scan_cfg['step'])

    df_flight = scan_threshold_cut(sig[flight_var], sig[mass_var],
                                   bkg[flight_var], bkg[mass_var], bkg['weight'],
                                   flight_cuts, lo0, hi0, sideband_mult)
    best_flight = best_from_scan(df_flight)
    chosen_flight_cut = best_flight['cut']

    sel_sig = sig[flight_var] > chosen_flight_cut
    sel_bkg = bkg[flight_var] > chosen_flight_cut

    halfwidths = np.arange(mass_scan_cfg['lo'], mass_scan_cfg['hi'], mass_scan_cfg['step'])
    df_mass = scan_mass_halfwidth(sig[mass_var][sel_sig], bkg[mass_var][sel_bkg], bkg['weight'][sel_bkg],
                                  halfwidths, mass_pdg, sideband_mult)
    best_mass = best_from_scan(df_mass)
    chosen_hw = best_mass['halfwidth']

    flight_at_boundary = is_at_scan_boundary(df_flight, 'cut', best_flight)
    mass_at_boundary = is_at_scan_boundary(df_mass, 'halfwidth', best_mass)
    if flight_at_boundary:
        print(f"WARNING: flight-cut FOM optimum ({chosen_flight_cut}) is at the edge of the "
              f"scanned range [{flight_scan_cfg['lo']}, {flight_scan_cfg['hi']}) -- the true "
              f"optimum may lie outside this range, or the FOM may simply be monotonic across "
              f"it in this sample. Do not trust this value without checking the scan plot.")
    if mass_at_boundary:
        print(f"WARNING: mass-halfwidth FOM optimum ({chosen_hw}) is at the edge of the "
              f"scanned range [{mass_scan_cfg['lo']}, {mass_scan_cfg['hi']}) -- same caveat "
              f"as above (check for an upstream skim cut depleting the sidebands, in "
              f"particular).")

    return {
        'df_flight': df_flight, 'df_mass': df_mass,
        'chosen_flight_cut': float(chosen_flight_cut),
        'chosen_mass_window': [float(mass_pdg - chosen_hw), float(mass_pdg + chosen_hw)],
        'fom_flight': float(best_flight['fom']), 'fom_mass': float(best_mass['fom']),
        'flight_at_scan_boundary': bool(flight_at_boundary),
        'mass_halfwidth_at_scan_boundary': bool(mass_at_boundary),
    }
################################################################################
