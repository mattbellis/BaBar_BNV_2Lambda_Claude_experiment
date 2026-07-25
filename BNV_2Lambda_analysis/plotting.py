"""
Histogramming and plotting functions for the BNV 2-Lambda analyses.

Histograms are hist.Hist objects with a regular axis for the variable and
growing string-category axes for SP mode and cut, following the p-Lambda0
reference analysis.

Plots are written to plots/<channel>/ (created on demand).
"""

import os

import awkward as ak
import numpy as np
import pandas as pd

import matplotlib.pylab as plt

import hist
from hist import Hist

from channel_config import SP_MODE_INFO, BACKGROUND_SP_MODES, SIGNAL_SP_MODE, DATA_SP_MODE

################################################################################
def plot_dir(config, topdir='plots'):
    """Return (and create if needed) the plot directory for a channel."""
    path = os.path.join(topdir, config['name'])
    os.makedirs(path, exist_ok=True)
    return path

################################################################################
def create_empty_histograms(hist_defs):
    """
    Create empty Hist objects (variable x SP mode x cut, weighted storage)
    from a {varname: {nbins, lo, hi, label}} dictionary.
    """
    all_hists = {}
    for var in hist_defs.keys():
        h = Hist.new.Reg(hist_defs[var]["nbins"], hist_defs[var]["lo"], hist_defs[var]["hi"],
                         name='var', label=f"{hist_defs[var]['label']}") \
                 .StrCat([], name="SP", label="SP modes", growth=True) \
                 .StrCat([], name="cuts", label="Cuts", growth=True) \
                 .Weight()

        all_hists[var] = h

    return all_hists

################################################################################
def fill_histograms(data, all_hists, dcuts, spmodes=None, weights=None, subset=None, verbose=False):
    """
    Fill the histograms for each variable / SP mode / cut combination.

    data:      awkward array (MC or collision)
    all_hists: dict of Hist objects from create_empty_histograms
    dcuts:     dict of cuts (see cutflow.build_diagnostic_cuts); the event
               masks are used, and jagged variables are flattened after
               the event selection
    spmodes:   list of SP-mode strings to fill (default: all in the array)
    weights:   dict {spmode: weight} (default: weight 1 for everything)
    subset:    list of variable names to fill (default: all in all_hists)

    Returns a DataFrame recording the number of entries per fill.
    """
    if spmodes is None:
        spmodes = list(np.unique(np.array(data['spmode'].to_list())))

    if weights is None:
        weights = {str(spmode): 1.0 for spmode in spmodes}

    if subset is None:
        subset = list(all_hists.keys())

    df_dict = {'var': [], 'cut': [], 'spmode': [], 'n': []}

    for key in subset:
        if verbose:
            print(key)

        # Skip variables not present in this file
        if key not in data.fields:
            print(f"WARNING: {key} not found in the array; skipping")
            continue

        for spmode in spmodes:
            weight = weights[str(spmode)]
            mask_sp = data['spmode'] == spmode

            for cutname in dcuts.keys():
                cut = dcuts[cutname]['event']

                x = data[key][cut & mask_sp]

                if len(x) <= 0:
                    continue

                # Flatten jagged (candidate-level) variables
                try:
                    float(x[0])
                except (TypeError, ValueError):
                    x = ak.flatten(x)

                n = len(x)
                if n > 0:
                    all_hists[key].fill(var=x, SP=str(spmode), cuts=f"{cutname}", weight=weight)

                df_dict['var'].append(key)
                df_dict['cut'].append(cutname)
                df_dict['spmode'].append(str(spmode))
                df_dict['n'].append(n)

    return pd.DataFrame.from_dict(df_dict)

################################################################################
def plot_stacked(all_hists, var, cut='0', bkg_spmodes=None, sig_spmode=SIGNAL_SP_MODE,
                 data_spmode=DATA_SP_MODE, overlay_data=True, overlay_signal=True,
                 sig_norm='match_peak', sig_norm_label='(arb. norm.)', logy=False,
                 ax=None):
    """
    Stacked background MC for one variable/cut, with signal MC and collision
    data overlaid.

    The signal MC normalization on these plots is arbitrary: by default
    (sig_norm='match_peak') it is rescaled, per panel, so that its peak bin
    matches the peak of the summed stacked background. Pass sig_norm=None
    to draw the signal with the weights it was filled with.

    Only SP modes actually present in the histogram are drawn.
    """
    h = all_hists[var]

    if bkg_spmodes is None:
        bkg_spmodes = BACKGROUND_SP_MODES

    # Only categories that were actually filled
    filled = list(h.axes['SP'])
    bkg_spmodes = [s for s in bkg_spmodes if s in filled]

    if ax is None:
        plt.figure(figsize=(6, 4))
        ax = plt.gca()
    plt.sca(ax)

    if len(bkg_spmodes) > 0:
        hstack = h[:, bkg_spmodes, cut].stack('SP')
        # Relabel with human-readable names
        labels = [SP_MODE_INFO.get(s, {}).get('label', s) for s in bkg_spmodes]
        hstack[:].project('var').plot(stack=True, histtype='fill', label=labels)

    if overlay_signal and sig_spmode in filled:
        hsig = h[:, sig_spmode, cut].project('var')

        # Rescale the signal for visibility (its normalization is arbitrary
        # on these diagnostic plots): peak bin -> peak of the summed background
        if sig_norm == 'match_peak' and len(bkg_spmodes) > 0:
            bkg_max = h[:, bkg_spmodes, cut].project('var').values().max()
            sig_max = hsig.values().max()
            if sig_max > 0 and bkg_max > 0:
                hsig = hsig * (bkg_max / sig_max)

        hsig.plot(
            histtype='step', color=SP_MODE_INFO[sig_spmode]['color'], linewidth=1.5,
            label=f"{SP_MODE_INFO[sig_spmode]['label']} {sig_norm_label}")

    if overlay_data and data_spmode in filled:
        h[:, data_spmode, cut].project('var').plot(
            histtype='errorbar', color='black', markersize=4,
            label=SP_MODE_INFO[data_spmode]['label'])

    plt.legend(loc='best', fontsize=8)
    plt.xlabel(plt.gca().get_xlabel(), fontsize=14)
    plt.ylabel('Entries / bin')
    if logy:
        plt.yscale('log')

################################################################################
def plot_stacked_grid(all_hists, vars=None, cut='0', ncols=4, save=False,
                      config=None, extra_tag="", **kwargs):
    """
    Grid of plot_stacked panels, one per variable.
    If save is True, writes plots/<channel>/diagnostics_cut<cut><extra_tag>.png
    """
    if vars is None:
        vars = list(all_hists.keys())

    nrows = int(np.ceil(len(vars) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.2 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for ax_extra in axes[len(vars):]:
        ax_extra.set_visible(False)

    for ax, var in zip(axes, vars):
        plot_stacked(all_hists, var, cut=cut, ax=ax, **kwargs)

    plt.tight_layout()

    if save and config is not None:
        path = plot_dir(config)
        outfilename = f"diagnostics_cut{cut}{extra_tag}.png"
        plt.savefig(os.path.join(path, outfilename), dpi=150)
        print(f"Saved {os.path.join(path, outfilename)}")

################################################################################
def plot_split_by_mode(x, mode, hist_def, mode_labels, density=True, logy=False,
                       ax=None, title=None):
    """
    Overlay the distribution of `x` for each category in `mode` (e.g. the
    LambdaC decay modes). `x` and `mode` must be flat, aligned arrays
    (flatten candidate-level jagged arrays first).

    density=True normalizes each mode's histogram to unit area, which is
    usually what you want when the modes have very different yields.
    """
    x = np.asarray(ak.to_numpy(x))
    mode = np.asarray(ak.to_numpy(mode))

    if ax is None:
        plt.figure(figsize=(6, 4))
        ax = plt.gca()
    plt.sca(ax)

    bins = np.linspace(hist_def['lo'], hist_def['hi'], hist_def['nbins'] + 1)

    for m, label in mode_labels.items():
        xm = x[mode == m]
        if len(xm) == 0:
            continue
        plt.hist(xm, bins=bins, histtype='step', linewidth=1.5,
                 density=density, label=f"{m}: {label}")

    plt.xlabel(hist_def['label'], fontsize=14)
    plt.ylabel('arb. units' if density else 'Entries / bin')
    if logy:
        plt.yscale('log')
    if title is not None:
        plt.title(title)
    plt.legend(loc='best', fontsize=8)

################################################################################
def plot_mes_vs_DeltaE(mes, DeltaE, config, draw_signal_region=False, draw_sidebands=False,
                       draw_fit_region=False, bins=100, zoom=False, title=None,
                       tag=None, verbose=False):
    """
    2D histogram of mES vs DeltaE for (flattened) candidate arrays, with
    optional region boxes drawn from the channel's region definitions.

    Returns a dict of the counts in the signal / fit / sideband regions.

    BLINDING NOTE: for collision data this should only ever be used with the
    _BLINDED files; the signal-region count for data must be 0 -- notebook 01
    uses this as a verification of the blinding.
    """
    region_definitions = config['region_definitions']

    meslo, meshi = region_definitions['fitting MES']
    DeltaElo, DeltaEhi = -0.5, 0.5

    if zoom:
        DeltaElo, DeltaEhi = region_definitions['fitting DeltaE']

    h = Hist(
        hist.axis.Regular(bins, meslo, meshi, name="mes", label=r"$M_{ES}$ [GeV/c$^2$]", flow=True),
        hist.axis.Regular(bins, DeltaElo, DeltaEhi, name="de", label=r"$\Delta E$ [GeV]", flow=True),
    )

    h.fill(mes, DeltaE)

    h.plot2d(cmap='plasma')

    def draw_box(xlo, xhi, ylo, yhi, fmt='w--'):
        plt.plot([xlo, xlo, xhi, xhi, xlo], [ylo, yhi, yhi, ylo, ylo], fmt, linewidth=2)

    sigmeslo, sigmeshi = region_definitions['signal MES']
    sigDeltaElo, sigDeltaEhi = region_definitions['signal DeltaE']

    if draw_signal_region:
        draw_box(sigmeslo, sigmeshi, sigDeltaElo, sigDeltaEhi)

    if draw_sidebands:
        sbmeslo, sbmeshi = region_definitions['sideband MES']
        draw_box(sbmeslo, sbmeshi, *region_definitions['sideband 1 DeltaE'], fmt='y--')
        draw_box(sbmeslo, sbmeshi, *region_definitions['sideband 2 DeltaE'], fmt='y--')

    if draw_fit_region:
        fitDeltaElo, fitDeltaEhi = region_definitions['fitting DeltaE']
        plt.plot([meslo, meshi], [fitDeltaElo, fitDeltaElo], 'w:')
        plt.plot([meslo, meshi], [fitDeltaEhi, fitDeltaEhi], 'w:')

    plt.xlabel(plt.gca().get_xlabel(), fontsize=14)
    plt.ylabel(plt.gca().get_ylabel(), fontsize=14)
    if title is not None:
        plt.title(title)

    # Region counts
    signal_mask = (mes > sigmeslo) & (mes < sigmeshi) & (DeltaE > sigDeltaElo) & (DeltaE < sigDeltaEhi)

    fitDeltaElo, fitDeltaEhi = region_definitions['fitting DeltaE']
    fit_mask = (mes > meslo) & (mes < meshi) & (DeltaE > fitDeltaElo) & (DeltaE < fitDeltaEhi)

    sbmeslo, sbmeshi = region_definitions['sideband MES']
    sb1lo, sb1hi = region_definitions['sideband 1 DeltaE']
    sb2lo, sb2hi = region_definitions['sideband 2 DeltaE']
    sideband1_mask = (mes > sbmeslo) & (mes < sbmeshi) & (DeltaE > sb1lo) & (DeltaE < sb1hi)
    sideband2_mask = (mes > sbmeslo) & (mes < sbmeshi) & (DeltaE > sb2lo) & (DeltaE < sb2hi)

    d = {}
    d['nsig'] = int(ak.sum(signal_mask))
    d['nfit'] = int(ak.sum(fit_mask))
    d['nside1'] = int(ak.sum(sideband1_mask))
    d['nside2'] = int(ak.sum(sideband2_mask))

    if verbose:
        print(f"signal: {d['nsig']}   fit: {d['nfit']}  "
              f"sb1: {d['nside1']}  sb2: {d['nside2']}")

    if tag is not None:
        path = plot_dir(config)
        plt.savefig(os.path.join(path, f"mes_vs_de_{tag}.png"), dpi=150)

    return d
################################################################################

################################################################################
# Stage 2: purity-optimization scan / fit plots
################################################################################
def plot_fom_scan(df, xcol, chosen_x=None, xlabel='cut value', title=None, ax=None):
    """
    Standard 3-panel view of a purity_optimization scan DataFrame (S, B vs.
    the scanned variable; and the FOM itself), with the chosen operating
    point marked.
    """
    if ax is None:
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.5))

    ax[0].plot(df[xcol], df['S'], 'o-', label='S (sig. MC, sideband-sub.)')
    ax[0].plot(df[xcol], df['B'], 's-', label='B (bkg MC, sideband est.)')
    ax[0].set_xlabel(xlabel)
    ax[0].set_ylabel('weighted counts')
    ax[0].legend(fontsize=8)

    ax[1].plot(df[xcol], df['fom'], 'o-', color='tab:green')
    ax[1].set_xlabel(xlabel)
    ax[1].set_ylabel(r'FOM $= S/\sqrt{S+B}$')

    if 'sig_eff' in df.columns:
        ax[2].plot(df[xcol], df['sig_eff'], 'o-', color='tab:purple')
        ax[2].set_xlabel(xlabel)
        ax[2].set_ylabel('signal efficiency')
    else:
        ax[2].plot(df['B'], df['S'], 'o-', color='tab:gray')
        ax[2].set_xlabel('B')
        ax[2].set_ylabel('S')

    if chosen_x is not None:
        for a in ax:
            a.axvline(chosen_x, color='k', linestyle='--', linewidth=1)

    if title is not None:
        plt.suptitle(title)
    plt.tight_layout()

################################################################################
def plot_mass_fit(fit, hist_def=None, title=None, ax=None):
    """
    Binned mass distribution with the Gaussian-plus-linear-background fit
    overlaid, from purity_optimization.fit_mass_peak's return dict.
    """
    if ax is None:
        plt.figure(figsize=(5, 3.5))
        ax = plt.gca()
    plt.sca(ax)

    width = fit['bin_centers'][1] - fit['bin_centers'][0]
    ax.bar(fit['bin_centers'], fit['bin_counts'], width=width, alpha=0.5, label='signal MC')

    if fit['converged']:
        ax.plot(fit['bin_centers'], fit['fit_curve'], 'r-', linewidth=1.5,
                label=fr"$\mu$={fit['mu']*1000:.1f} MeV, $\sigma$={fit['sigma']*1000:.1f} MeV")
    else:
        ax.text(0.05, 0.9, 'fit did not converge', color='r', transform=ax.transAxes)

    ax.set_xlabel(hist_def['label'] if hist_def is not None else 'mass [GeV/c$^2$]')
    ax.set_ylabel('entries / bin')
    ax.legend(fontsize=8)
    if title is not None:
        ax.set_title(title)
################################################################################
