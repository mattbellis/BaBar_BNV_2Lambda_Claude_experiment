#!/usr/bin/env python
"""
Generate LaTeX macros and tables for the BAD from results/<channel>.yaml.

Reads the per-channel results files (written by run_stage_N.py scripts) and
writes into BAD_2Lambda/:

    generated_numbers.tex             \\newcommand macros for all key numbers
    generated_table_samples_<ch>.tex  MC/data samples and weights
    generated_table_cutflow_<ch>.tex  diagnostic cutflow
    generated_table_lambdac_modes.tex LambdaC decay-mode fractions
    generated_table_dataskims.tex     per-run luminosity / skim statistics

Numbers are NEVER hand-typed into the LaTeX -- rerun this after any stage
updates a results file:

    python generate_latex_macros.py
"""

import os

from channel_config import CHANNELS, SP_MODE_INFO
import datasets
import results_io

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
BAD_DIR = os.path.normpath(os.path.join(_THIS_DIR, '..', 'BAD_2Lambda'))

# LaTeX macro names cannot contain digits -- translation tables
CHANNEL_PREFIX = {'Lam0Lam0': 'LamLam', 'Lam0LamC': 'LamLamC'}

SPMODE_NAME = {
    '-999': 'Signal',
    '998': 'Uds',
    '1005': 'Ccbar',
    '1235': 'BpBm',
    '1237': 'BzBzbar',
    '3429': 'Tautau',
    '0': 'Data',
}

CUT_NAME = {0: 'CutZero', 1: 'CutOne', 2: 'CutTwo'}

MODE_NAME = {1: 'One', 2: 'Two', 3: 'Three', 4: 'Four'}

################################################################################
def fmt_int(v):
    """Integers with thin-space thousands separators (LaTeX-safe)."""
    return f"{int(v):,}".replace(",", "\\,")

def fmt_float(v, fmt="{:.3g}"):
    return fmt.format(float(v))

################################################################################
def macros_for_channel(channel, res):
    """Return a list of (macro_name, value_string) for one channel."""
    ch = CHANNEL_PREFIX[channel]
    s = res.get('stage01')
    if s is None:
        return []

    m = []

    # Event counts
    m.append((f"{ch}NDataEvents", fmt_int(s['n_events_data'])))
    for spmode, n in s['n_events_mc'].items():
        m.append((f"{ch}NMC{SPMODE_NAME[spmode]}", fmt_int(n)))

    # Scaling weights
    for spmode, w in s['scaling_weights'].items():
        m.append((f"{ch}Weight{SPMODE_NAME[spmode]}", fmt_float(w)))

    # Cutflow
    for spmode, cuts in s['cutflow_mc'].items():
        for cut, n in cuts.items():
            m.append((f"{ch}Cutflow{SPMODE_NAME[spmode]}{CUT_NAME[cut]}", fmt_int(n)))
    for cut, n in s['cutflow_data'].items():
        m.append((f"{ch}CutflowData{CUT_NAME[cut]}", fmt_int(n)))

    # Region definitions [GeV]
    rd = s['region_definitions']
    m.append((f"{ch}SigMesLo", fmt_float(rd['signal MES'][0], "{:.2f}")))
    m.append((f"{ch}SigMesHi", fmt_float(rd['signal MES'][1], "{:.2f}")))
    m.append((f"{ch}SigDeltaELo", fmt_float(rd['signal DeltaE'][0], "{:.2f}")))
    m.append((f"{ch}SigDeltaEHi", fmt_float(rd['signal DeltaE'][1], "{:.2f}")))
    m.append((f"{ch}FitMesLo", fmt_float(rd['fitting MES'][0], "{:.2f}")))
    m.append((f"{ch}FitMesHi", fmt_float(rd['fitting MES'][1], "{:.2f}")))
    m.append((f"{ch}FitDeltaELo", fmt_float(rd['fitting DeltaE'][0], "{:.2f}")))
    m.append((f"{ch}FitDeltaEHi", fmt_float(rd['fitting DeltaE'][1], "{:.2f}")))

    # Data candidate counts by region + blinding check
    rc = s['region_counts_data']
    m.append((f"{ch}NDataCandSignalRegion", fmt_int(rc['signal'])))
    m.append((f"{ch}NDataCandFitRegion", fmt_int(rc['fit'])))
    m.append((f"{ch}NDataCandSidebandOne", fmt_int(rc['sideband1'])))
    m.append((f"{ch}NDataCandSidebandTwo", fmt_int(rc['sideband2'])))

    # Composite mass windows (half-width in MeV)
    for name, comp in s['composites'].items():
        lo, hi = comp['mass_window']
        half_mev = 1000.0 * (hi - lo) / 2.0
        cname = 'LamZero' if name == 'Lambda0' else 'LamC'
        m.append((f"{ch}{cname}MassWindowMeV", fmt_float(half_mev, "{:.1f}")))
        if comp['flight_cut'] is not None:
            m.append((f"{ch}{cname}FlightSigCut", fmt_float(comp['flight_cut'], "{:.0f}")))

    # LambdaC decay-mode fractions (percent, no % sign). Iterate over ALL
    # modes: one with zero surviving candidates is absent from the yaml
    # dict but its macro must still exist (value 0.0).
    if 'lambdac_mode_fractions' in s:
        for tag, mac in [('no_cut', 'NoCut'), ('single_candidate', 'SingleCand')]:
            for mode, name in MODE_NAME.items():
                frac = s['lambdac_mode_fractions'][tag].get(mode, 0.0)
                m.append((f"{ch}ModeFrac{name}{mac}Pct",
                          fmt_float(100 * frac, "{:.1f}")))

    return m

################################################################################
def macros_for_channel_stage02(channel, res):
    """Stage 2 (purity optimization) macros: recommended cut values, LambdaC
    per-mode fit results, and the multi-candidate study (Lam0LamC only)."""
    ch = CHANNEL_PREFIX[channel]
    s = res.get('stage02')
    if s is None:
        return []

    m = []

    # Lambda0: recommended flight cut / mass half-width (both channels)
    m.append((f"{ch}LamZeroRecFlightSigCut", fmt_float(s['lambda0']['recommended_flight_cut'], "{:.0f}")))
    lo, hi = s['lambda0']['recommended_mass_window']
    m.append((f"{ch}LamZeroRecMassWindowMeV", fmt_float(1000 * (hi - lo) / 2.0, "{:.1f}")))

    if 'k0s' in s:
        m.append((f"{ch}KZeroRecFlightSigCut", fmt_float(s['k0s']['recommended_flight_cut'], "{:.0f}")))
        lo, hi = s['k0s']['recommended_mass_window']
        m.append((f"{ch}KZeroRecMassWindowMeV", fmt_float(1000 * (hi - lo) / 2.0, "{:.1f}")))

    if 'lambdac_mode_fits' in s:
        m.append((f"{ch}LamCMassWindowNSigma", fmt_float(s['lambdac_mass_window_nsigma'], "{:.1f}")))
        for mode, name in MODE_NAME.items():
            fit = s['lambdac_mode_fits'].get(mode)
            if fit is None:
                continue
            m.append((f"{ch}LamCMode{name}SigmaMeV", fmt_float(1000 * fit['sigma'], "{:.1f}")))
            m.append((f"{ch}LamCMode{name}MuMeV", fmt_float(1000 * fit['mu'], "{:.1f}")))
            lo, hi = s['lambdac_recommended_windows_per_mode'][mode]
            m.append((f"{ch}LamCMode{name}RecWindowMeV", fmt_float(1000 * (hi - lo) / 2.0, "{:.1f}")))

    if 'multi_candidate_study' in s:
        mc = s['multi_candidate_study']
        m.append((f"{ch}FracZeroGoodBPct", fmt_float(100 * mc['frac_zero_good_B'], "{:.1f}")))
        m.append((f"{ch}FracOneGoodBPct", fmt_float(100 * mc['frac_one_good_B'], "{:.1f}")))
        m.append((f"{ch}FracMultiGoodBPct", fmt_float(100 * mc['frac_multi_good_B'], "{:.1f}")))

    return m

################################################################################
def _pid_ladder_result_macros(prefix, summary):
    """Macros shared by any one Stage 3 ladder-scan result (Lambda0, or one
    LambdaC mode): recommended selector name(s) + FOM/efficiency/background/
    boundary flag. `prefix` is the full macro prefix (channel + target)."""
    m = []
    for particle, selector in summary['recommended_selector'].items():
        pname = {'p': 'P', 'K': 'K', 'pi': 'Pi'}[particle]
        m.append((f"{prefix}{pname}Selector", selector))
    m.append((f"{prefix}Fom", fmt_float(summary['fom_at_recommended'])))
    m.append((f"{prefix}SigEffPct", fmt_float(100 * summary['sig_eff_at_recommended'], "{:.1f}")))
    m.append((f"{prefix}BkgWeighted", fmt_float(summary['bkg_weighted_at_recommended'], "{:.2f}")))
    m.append((f"{prefix}AtBoundary", "yes" if summary['any_at_scan_boundary'] else "no"))
    if 'fom_no_pid_at_all' in summary:
        m.append((f"{prefix}FomNoPidAtAll", fmt_float(summary['fom_no_pid_at_all'])))
        m.append((f"{prefix}BkgNoPidAtAll", fmt_float(summary['bkg_weighted_no_pid_at_all'], "{:.2f}")))
    return m

################################################################################
def macros_for_channel_stage03(channel, res):
    """Stage 3 (PID selector optimization) macros: recommended KM-ladder
    selector(s), Punzi FOM/efficiency/background, and boundary flag, for the
    Lambda0 scan (both channels) and each LambdaC decay mode (Lam0LamC)."""
    ch = CHANNEL_PREFIX[channel]
    s = res.get('stage03')
    if s is None:
        return []

    m = _pid_ladder_result_macros(f"{ch}LamZeroPid", s['lambda0_pid'])

    if 'lambdac_pid_per_mode' in s:
        for mode, name in MODE_NAME.items():
            summary = s['lambdac_pid_per_mode'].get(mode)
            if summary is None:
                continue
            m += _pid_ladder_result_macros(f"{ch}LamCMode{name}Pid", summary)

    if 'multi_candidate_study' in s:
        mc = s['multi_candidate_study']
        m.append((f"{ch}FracZeroGoodBPidPct", fmt_float(100 * mc['frac_zero_good_B'], "{:.1f}")))
        m.append((f"{ch}FracOneGoodBPidPct", fmt_float(100 * mc['frac_one_good_B'], "{:.1f}")))
        m.append((f"{ch}FracMultiGoodBPidPct", fmt_float(100 * mc['frac_multi_good_B'], "{:.1f}")))

    return m

################################################################################
def write_table_pid_selectors(all_results):
    """Stage 3 recommended KM-ladder selector(s) and Punzi FOM/efficiency/
    background per PID target (Lambda0, both channels; LambdaC per mode,
    Lam0LamC), with the boundary-hugging flag surfaced rather than hidden."""
    rows = []  # (label, summary)

    for channel in ('Lam0Lam0', 'Lam0LamC'):
        res = all_results.get(channel, {})
        s = res.get('stage03')
        if s is None:
            continue
        label = CHANNELS[channel]['decay_label']
        rows.append((rf"{label}: $\Lambda^0$", s['lambda0_pid']))
        if 'lambdac_pid_per_mode' in s:
            mode_labels = CHANNELS[channel]['lambdac_modes']
            for mode, summary in s['lambdac_pid_per_mode'].items():
                rows.append((rf"{label}: $\Lambda_c^+$ mode {mode} ({mode_labels[mode]})", summary))

    if not rows:
        return

    lines = [r"% AUTO-GENERATED -- DO NOT EDIT BY HAND",
             r"\begin{table}[h]",
             r"\caption{Stage 3 recommended KM-ladder PID selector(s), the"
             r" Punzi FOM ($a=4$) and signal efficiency/background at that"
             r" operating point, and whether the scan optimum sits at the"
             r" edge of the scanned ladder (see text before applying).}",
             r"\label{tab:pidselectors}",
             r"\centering",
             r"\small",
             r"\resizebox{\textwidth}{!}{%",
             r"\begin{tabular}{lp{5.5cm}rrrc}",
             r"\toprule",
             r"Target & Recommended selector(s) & FOM & $\epsilon_{sig}$ [\%] & "
             r"$B$ (wtd.) & boundary? \\",
             r"\midrule"]

    for label, summary in rows:
        sel_str = ", ".join(f"{p}: {n.replace('KMProtonSelection','').replace('KMKaonMicroSelection','').replace('KMPionMicroSelection','')}"
                            for p, n in summary['recommended_selector'].items())
        lines.append(
            rf"{label} & {sel_str} & {summary['fom_at_recommended']:.3f} & "
            rf"{100*summary['sig_eff_at_recommended']:.1f} & {summary['bkg_weighted_at_recommended']:.2f} & "
            rf"{'yes' if summary['any_at_scan_boundary'] else 'no'} \\")

    lines += [r"\bottomrule", r"\end{tabular}", r"}", r"\end{table}"]

    outname = os.path.join(BAD_DIR, "generated_table_pid_selectors.tex")
    with open(outname, 'w') as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {outname}")

################################################################################
def write_table_lambdac_mode_fits(all_results):
    """LambdaC per-mode mass-resolution fit results and recommended windows."""
    channel = 'Lam0LamC'
    res = all_results.get(channel, {})
    if 'stage02' not in res or 'lambdac_mode_fits' not in res['stage02']:
        return
    s = res['stage02']

    mode_labels = CHANNELS[channel]['lambdac_modes']

    lines = [r"% AUTO-GENERATED -- DO NOT EDIT BY HAND",
             r"\begin{table}[h]",
             r"\caption{Per-mode $\Lambda_c^+$ mass fit (Gaussian + linear"
             r" background, signal MC, $K_S^0$ gate applied to modes 2 and"
             r" 3) and the resulting mass\_window\_nsigma$\times\sigma$"
             r" window.}",
             r"\label{tab:lambdacmodefits}",
             r"\centering",
             r"\begin{tabular}{clrrrr}",
             r"\toprule",
             r"Mode & Decay & $\mu$ [MeV] & $\sigma$ [MeV] & converged & "
             r"window [MeV] \\",
             r"\midrule"]

    for mode, label in mode_labels.items():
        fit = s['lambdac_mode_fits'].get(mode)
        if fit is None:
            continue
        lo, hi = s['lambdac_recommended_windows_per_mode'][mode]
        halfwidth = 1000 * (hi - lo) / 2.0
        lines.append(rf"{mode} & {label} & {1000*fit['mu']:.1f} & {1000*fit['sigma']:.1f} & "
                     rf"{'yes' if fit['converged'] else 'NO'} & $\pm${halfwidth:.1f} \\")

    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]

    outname = os.path.join(BAD_DIR, "generated_table_lambdac_mode_fits.tex")
    with open(outname, 'w') as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {outname}")

################################################################################
def write_macros(all_results):
    lines = ["% AUTO-GENERATED by BNV_2Lambda_analysis/generate_latex_macros.py",
             "% from results/<channel>.yaml -- DO NOT EDIT BY HAND",
             ""]

    # Shared: integrated luminosity (same dataset for both channels)
    for channel, res in all_results.items():
        if 'stage01' in res:
            lumi_invfb = res['stage01']['int_lumi_invpb'] / 1000.0
            lines.append(f"\\newcommand{{\\IntLumiInvFb}}{{{lumi_invfb:.1f}}}")
            break

    for channel, res in all_results.items():
        lines.append(f"% ---- {channel} ----")
        for name, value in macros_for_channel(channel, res):
            lines.append(f"\\newcommand{{\\{name}}}{{{value}}}")
        for name, value in macros_for_channel_stage02(channel, res):
            lines.append(f"\\newcommand{{\\{name}}}{{{value}}}")
        for name, value in macros_for_channel_stage03(channel, res):
            lines.append(f"\\newcommand{{\\{name}}}{{{value}}}")

    outname = os.path.join(BAD_DIR, 'generated_numbers.tex')
    with open(outname, 'w') as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {outname}  ({len(lines)} lines)")

################################################################################
def write_table_samples(channel, res):
    s = res['stage01']
    ch_label = CHANNELS[channel]['decay_label']

    lines = [r"% AUTO-GENERATED -- DO NOT EDIT BY HAND",
             r"\begin{table}[h]",
             rf"\caption{{Monte Carlo and data samples for the {ch_label} channel:"
             r" raw numbers of events in the analysis files and the"
             r" luminosity-scaling weights applied to background MC.}",
             rf"\label{{tab:samples{CHANNEL_PREFIX[channel]}}}",
             r"\centering",
             r"\begin{tabular}{llrr}",
             r"\toprule",
             r"SP mode & Description & \# of events & weight \\",
             r"\midrule"]

    order = ['-999', '998', '1005', '1235', '1237', '3429']
    for spmode in order:
        if spmode not in s['n_events_mc']:
            continue
        label = SP_MODE_INFO[spmode]['label']
        n = fmt_int(s['n_events_mc'][spmode])
        w = fmt_float(s['scaling_weights'][spmode]) if spmode in s['scaling_weights'] \
            else '--'
        lines.append(rf"{spmode} & {label} & {n} & {w} \\")

    lines += [r"\midrule",
              rf"data & collision data (blinded) & {fmt_int(s['n_events_data'])} & -- \\",
              r"\bottomrule",
              r"\end{tabular}",
              r"\end{table}"]

    outname = os.path.join(BAD_DIR, f"generated_table_samples_{channel}.tex")
    with open(outname, 'w') as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {outname}")

################################################################################
def write_table_cutflow(channel, res):
    s = res['stage01']
    ch_label = CHANNELS[channel]['decay_label']

    spmodes = [sp for sp in ['-999', '998', '1005', '1235', '1237', '3429']
               if sp in s['cutflow_mc']]
    cuts = sorted(s['cut_names'].keys())

    colspec = "ll" + "r" * (len(spmodes) + 1)
    header = "Cut & Description & " + \
             " & ".join(SP_MODE_INFO[sp]['label'] for sp in spmodes) + \
             r" & data \\"

    lines = [r"% AUTO-GENERATED -- DO NOT EDIT BY HAND",
             r"\begin{table}[h]",
             rf"\caption{{Numbers of events (raw, unweighted) passing the"
             rf" diagnostic cuts for the {ch_label} channel.}}",
             rf"\label{{tab:cutflow{CHANNEL_PREFIX[channel]}}}",
             r"\centering",
             r"\small",
             rf"\begin{{tabular}}{{{colspec}}}",
             r"\toprule",
             header,
             r"\midrule"]

    for cut in cuts:
        row = [str(cut), s['cut_names'][cut]]
        for sp in spmodes:
            row.append(fmt_int(s['cutflow_mc'][sp][cut]))
        row.append(fmt_int(s['cutflow_data'][cut]))
        lines.append(" & ".join(row) + r" \\")

    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]

    outname = os.path.join(BAD_DIR, f"generated_table_cutflow_{channel}.tex")
    with open(outname, 'w') as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {outname}")

################################################################################
def write_table_lambdac_modes(all_results):
    channel = 'Lam0LamC'
    res = all_results.get(channel, {})
    if 'stage01' not in res or 'lambdac_mode_fractions' not in res['stage01']:
        return
    s = res['stage01']

    mode_labels = CHANNELS[channel]['lambdac_modes']

    lines = [r"% AUTO-GENERATED -- DO NOT EDIT BY HAND",
             r"\begin{table}[h]",
             r"\caption{Fractions of $\Lambda_c^+$ candidates in signal MC"
             r" reconstructed in each decay mode, before any cuts and after"
             r" the single-candidate requirement. The suppression of mode 4"
             r" by the $n_{\Lambda^0} = 1$ requirement is discussed in the text.}",
             r"\label{tab:lambdacmodes}",
             r"\centering",
             r"\begin{tabular}{clrr}",
             r"\toprule",
             r"Mode & Decay & no cut [\%] & single cand. [\%] \\",
             r"\midrule"]

    for mode, label in mode_labels.items():
        f0 = 100 * s['lambdac_mode_fractions']['no_cut'].get(mode, 0.0)
        f1 = 100 * s['lambdac_mode_fractions']['single_candidate'].get(mode, 0.0)
        lines.append(rf"{mode} & {label} & {f0:.1f} & {f1:.1f} \\")

    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]

    outname = os.path.join(BAD_DIR, "generated_table_lambdac_modes.tex")
    with open(outname, 'w') as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {outname}")

################################################################################
def write_table_dataskims():
    """Per-run luminosity and skim statistics (from the dataset CSVs)."""
    df = datasets.read_in_dataset_statistics()
    dfsp = datasets.get_SP_cross_sections_and_labels()

    bbbar_xsec = 0.0
    for mode in (1235, 1237):
        bbbar_xsec += float(dfsp[dfsp['SP Mode'] == mode]['Cross section [nb]'].values[0])

    mask = (df['Data or MC'] == 'Data') & (df['Skim'] == 'LambdaVeryVeryLoose')
    dftmp = df[mask].copy()

    lines = [r"% AUTO-GENERATED -- DO NOT EDIT BY HAND",
             r"\begin{table}[h]",
             r"\caption{Integrated luminosity, numbers of skimmed events"
             r" ({\tt LambdaVeryVeryLoose} skim), and estimated numbers of"
             r" $B\bar{B}$ pairs per run period.}",
             r"\label{tab:dataskims}",
             r"\centering",
             r"\begin{tabular}{lrrr}",
             r"\toprule",
             r"Run & Luminosity [pb$^{-1}$] & \# skimmed events & \# $B\bar{B}$ pairs \\",
             r"\midrule"]

    tot_lumi, tot_n, tot_bb = 0.0, 0, 0.0
    for _, r in dftmp.iterrows():
        lumi = float(r['Luminosity (Data only) 1/pb'])
        n = int(r['# of events (Data or MC)'])
        nbb = lumi * bbbar_xsec * 1000
        tot_lumi += lumi
        tot_n += n
        tot_bb += nbb
        lines.append(rf"{int(r['Run'])} & {lumi:,.1f} & {n:,} & {nbb:,.0f} \\"
                     .replace(",", "\\,"))

    lines += [r"\midrule",
              (rf"Total & {tot_lumi:,.1f} & {tot_n:,} & {tot_bb:,.0f} \\"
               .replace(",", "\\,")),
              r"\bottomrule", r"\end{tabular}", r"\end{table}"]

    outname = os.path.join(BAD_DIR, "generated_table_dataskims.tex")
    with open(outname, 'w') as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {outname}")

################################################################################
if __name__ == '__main__':
    all_results = {}
    for channel in CHANNEL_PREFIX:
        res = results_io.load_results(channel)
        if res:
            all_results[channel] = res
        else:
            print(f"NOTE: no results file for {channel} yet "
                  f"(run run_stage01.py --channel {channel})")

    if not all_results:
        raise SystemExit("No results files found -- nothing to generate.")

    write_macros(all_results)
    for channel, res in all_results.items():
        if 'stage01' in res:
            write_table_samples(channel, res)
            write_table_cutflow(channel, res)
    write_table_lambdac_modes(all_results)
    write_table_lambdac_mode_fits(all_results)
    write_table_pid_selectors(all_results)
    write_table_dataskims()

    print("\nDone. Remember to rerun this script whenever a results file changes.")
