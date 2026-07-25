"""
Channel configuration for the BNV searches

    B0 -> Lambda0 Lambda0   (channel tag 'Lam0Lam0')
    B+ -> LambdaC+ Lambda0  (channel tag 'Lam0LamC')

This module is the single source of truth for everything that differs between
the two channels: file tags, composite-candidate structure, mass windows,
mES/DeltaE region definitions, and histogram definitions.

Every analysis function should take the channel config (or the channel tag)
as an argument rather than hardcoding 'Lambda0' etc.

Usage:
    from channel_config import get_channel_config
    config = get_channel_config('Lam0Lam0')
"""

import copy

# PDG masses [GeV/c^2]
LAMBDA0_MASS_PDG = 1.115683
LAMBDAC_MASS_PDG = 2.28646
K0S_MASS_PDG = 0.497611

################################################################################
# Stage 2 (Lambda/K_S purity) figure-of-merit configuration
#
# The Lambda0 and K_S0 mass-window + flight-significance cuts are optimized
# with S/sqrt(S+B), where S is the (sideband-subtracted) signal-MC yield in
# the mass-peak window and B is the sideband-estimated background under the
# peak from luminosity-weighted background MC. The sideband bands are
# contiguous with the peak window and each have the same width as the peak
# (mirrors the B+ -> p Lambda0 reference purity study:
# BNV_pLambda/Lambda_purity_studies.ipynb). This mirroring is what
# 'FOM_SIDEBAND_WIDTH_MULT = 1.0' encodes: sideband width = 1x the (current,
# possibly-scanned) peak full width, on each side.
################################################################################
FOM_SIDEBAND_WIDTH_MULT = 1.0

def _lambdac_mode_window(mu, sigma, nsigma):
    return [mu - nsigma * sigma, mu + nsigma * sigma]

def _default_purity_scan_ranges():
    """
    Default flight-significance / mass-halfwidth scan grids for the Stage 2
    FOM optimizations. Shared shape for Lambda0 (both channels) and K_S0
    (Lam0LamC); values in GeV for the mass scan, dimensionless significance
    units for the flight scan.
    """
    return {
        'flight_scan': {'lo': 0.0, 'hi': 100.0, 'step': 2.0},
        'mass_halfwidth_scan': {'lo': 0.001, 'hi': 0.010, 'step': 0.0005},
    }

################################################################################
# SP mode information, shared by both channels
#
# 'color' is used for the individual-mode plots, 'label' in legends.
################################################################################
SP_MODE_INFO = {}
SP_MODE_INFO['-999'] = {'label': 'Signal MC',                     'color': 'gold'}
SP_MODE_INFO['998']  = {'label': r'$u\bar{u}/d\bar{d}/s\bar{s}$', 'color': 'tab:blue'}
SP_MODE_INFO['1005'] = {'label': r'$c\bar{c}$',                   'color': 'tab:orange'}
SP_MODE_INFO['1235'] = {'label': r'$B^+B^-$',                     'color': 'tab:green'}
SP_MODE_INFO['1237'] = {'label': r'$B^0\overline{B}^0$',          'color': 'tab:red'}
SP_MODE_INFO['3429'] = {'label': r'$\tau^+\tau^-$',               'color': 'tab:purple'}
SP_MODE_INFO['0']    = {'label': 'Data',                          'color': 'black'}

BACKGROUND_SP_MODES = ['998', '1005', '1235', '1237', '3429']
SIGNAL_SP_MODE = '-999'
DATA_SP_MODE = '0'

################################################################################
# mES / DeltaE region definitions
#
# These mirror the B+ -> p Lambda0 reference analysis for now.
# NOTE: the exact signal windows are to be finalized per channel (CLAUDE.md
# quotes |DeltaE| < 0.05 as approximate; the p-Lambda0 analysis used 0.07).
# Do NOT narrow the assumed signal window below what was used to blind the
# collision-data files -- verify against the blinded region at the Stage 1
# checkpoint (see notebook 01).
################################################################################
def _default_region_definitions():
    region_definitions = {}

    # Low and then high
    region_definitions['signal MES'] = [5.27, 5.3]
    region_definitions['fitting MES'] = [5.2, 5.3]

    region_definitions['signal DeltaE'] = [-0.07, 0.07]
    region_definitions['fitting DeltaE'] = [-0.2, 0.2]

    region_definitions['sideband 1 DeltaE'] = [0.07, 0.14]
    region_definitions['sideband 2 DeltaE'] = [-0.14, -0.07]
    region_definitions['sideband MES'] = [5.27, 5.3]

    # [meslo, meshi, dElo, dEhi] boxes for the Rolke-style inference regions
    region_definitions['inference'] = []
    region_definitions['inference'].append([5.27, 5.3, -0.07, 0.07])
    region_definitions['inference'].append([5.27, 5.3, -0.2, -0.07])
    region_definitions['inference'].append([5.27, 5.3,  0.07, 0.20])
    region_definitions['inference'].append([5.25, 5.27, -0.2, 0.0])
    region_definitions['inference'].append([5.25, 5.27,  0.0, 0.20])
    region_definitions['inference'].append([5.23, 5.25, -0.2, 0.0])
    region_definitions['inference'].append([5.23, 5.25,  0.0, 0.20])
    region_definitions['inference'].append([5.20, 5.23, -0.2, 0.0])
    region_definitions['inference'].append([5.20, 5.23,  0.0, 0.20])

    return region_definitions

################################################################################
# Histogram definitions shared by both channels
################################################################################
COMMON_HIST_DEFS = {}

# B-candidate kinematics
COMMON_HIST_DEFS['BpostFitMes']    = {'nbins': 100, 'lo': 5.2,   'hi': 5.3,   'label': r'$M_{ES}$ [GeV/c$^2$]'}
COMMON_HIST_DEFS['BpostFitDeltaE'] = {'nbins': 100, 'lo': -1.0,  'hi': 1.0,   'label': r'$\Delta E$ [GeV]'}
COMMON_HIST_DEFS['BpostFitMass']   = {'nbins': 100, 'lo': 4.0,   'hi': 6.0,   'label': r'$B$ mass (post-fit) [GeV/c$^2$]'}

# Lambda0 candidates (jagged; in Lam0Lam0 both Lambdas fill the same histogram)
COMMON_HIST_DEFS['Lambda0_unc_Mass']                 = {'nbins': 100, 'lo': 1.105, 'hi': 1.125, 'label': r'$\Lambda^0$ mass [GeV/c$^2$]'}
COMMON_HIST_DEFS['Lambda0FlightLen']                 = {'nbins': 100, 'lo': -1,    'hi': 40,    'label': r'$\Lambda^0$ flight length [cm]'}
COMMON_HIST_DEFS['Lambda0postFitFlight']             = {'nbins': 100, 'lo': -1,    'hi': 40,    'label': r'$\Lambda^0$ post-fit flight length [cm]'}
COMMON_HIST_DEFS['Lambda0postFitFlightSignificance'] = {'nbins': 100, 'lo': -20,   'hi': 300,   'label': r'$\Lambda^0$ flight significance'}
COMMON_HIST_DEFS['Lambda0p3CM']                      = {'nbins': 100, 'lo': 0,     'hi': 3.0,   'label': r'$\Lambda^0$ $p$ (CM) [GeV/c]'}

# Event-shape / continuum-suppression variables (MLP inputs, most likely)
COMMON_HIST_DEFS['BSphr']          = {'nbins': 100, 'lo': 0,     'hi': 0.2,  'label': 'B sphericity'}
COMMON_HIST_DEFS['BCosSphr']       = {'nbins': 100, 'lo': -0.8,  'hi': 1,    'label': 'BCosSphr'}
COMMON_HIST_DEFS['BThrust']        = {'nbins': 100, 'lo': 0.9,   'hi': 1.05, 'label': 'B thrust'}
COMMON_HIST_DEFS['BCosThrust']     = {'nbins': 100, 'lo': 0,     'hi': 1,    'label': 'BCosThrust'}
COMMON_HIST_DEFS['BCosThetaS']     = {'nbins': 150, 'lo': -1.05, 'hi': 1.05, 'label': r'B $\cos\theta_S$'}
COMMON_HIST_DEFS['BCosThetaT']     = {'nbins': 100, 'lo': -1,    'hi': 1,    'label': r'B $\cos\theta_T$'}
COMMON_HIST_DEFS['BLegendreP2']    = {'nbins': 100, 'lo': 0,     'hi': 7,    'label': 'BLegendreP2'}
COMMON_HIST_DEFS['BR2ROE']         = {'nbins': 100, 'lo': 0,     'hi': 1,    'label': 'BR2ROE'}
COMMON_HIST_DEFS['BSphrROE']       = {'nbins': 100, 'lo': 0,     'hi': 1,    'label': 'BSphrROE'}
COMMON_HIST_DEFS['BThrustROE']     = {'nbins': 100, 'lo': 0.5,   'hi': 1,    'label': 'BThrustROE'}
COMMON_HIST_DEFS['R2']             = {'nbins': 100, 'lo': 0,     'hi': 1.05, 'label': 'R2'}
COMMON_HIST_DEFS['R2All']          = {'nbins': 100, 'lo': 0,     'hi': 1,    'label': 'R2All'}
COMMON_HIST_DEFS['thrustMag']      = {'nbins': 100, 'lo': 0.6,   'hi': 1,    'label': 'thrustMag'}
COMMON_HIST_DEFS['thrustMagAll']   = {'nbins': 100, 'lo': 0.6,   'hi': 1,    'label': 'thrustMagAll'}
COMMON_HIST_DEFS['thrustCosTh']    = {'nbins': 100, 'lo': 0,     'hi': 1,    'label': 'thrustCosTh'}
COMMON_HIST_DEFS['thrustCosThAll'] = {'nbins': 100, 'lo': 0,     'hi': 1,    'label': 'thrustCosThAll'}
COMMON_HIST_DEFS['sphericityAll']  = {'nbins': 100, 'lo': 0,     'hi': 0.75, 'label': 'Sphericity (all)'}

# Multiplicities (event-level)
COMMON_HIST_DEFS['nTRK']          = {'nbins': 20, 'lo': 0, 'hi': 20, 'label': '# of charged tracks'}
COMMON_HIST_DEFS['nTracks']       = {'nbins': 20, 'lo': 0, 'hi': 20, 'label': 'nTracks'}
COMMON_HIST_DEFS['nGoodTrkLoose'] = {'nbins': 15, 'lo': 0, 'hi': 15, 'label': '# of good tracks (loose)'}
COMMON_HIST_DEFS['nB']            = {'nbins': 10, 'lo': 0, 'hi': 10, 'label': '# of B candidates'}
COMMON_HIST_DEFS['nLambda0']      = {'nbins': 10, 'lo': 0, 'hi': 10, 'label': r'# of $\Lambda^0$ candidates'}

################################################################################
# Channel definitions
################################################################################
CHANNELS = {}

# ------------------------------------------------------------------------
# B0 -> Lambda0 Lambda0
# ------------------------------------------------------------------------
CHANNELS['Lam0Lam0'] = {
    'name': 'Lam0Lam0',
    'decay_label': r'$B^0 \to \Lambda^0 \Lambda^0$',
    'B_Lund': 511,

    # Both B daughters are Lambda0; Bd1Idx/Bd2Idx point into the Lambda0 arrays
    'B_daughters': [('Lambda0', 'Bd1Idx'), ('Lambda0', 'Bd2Idx')],

    # Composite candidates: how many we require per event (single-candidate
    # selection) and which variables define their purity cuts.
    'composites': {
        'Lambda0': {
            'n_required': 2,
            'mass_var': 'Lambda0_unc_Mass',
            'mass_pdg': LAMBDA0_MASS_PDG,
            'mass_window': [LAMBDA0_MASS_PDG - 0.003, LAMBDA0_MASS_PDG + 0.003],
            'flight_var': 'Lambda0postFitFlightSignificance',
            # Stage 2 S/sqrt(S+B) scan: BOTH the flight-cut and mass-window
            # optima hit the edge of their scanned ranges for this channel
            # (see STATUS.md open issue) -- neither is applied; stays at
            # the p-Lambda0 default pending that fix.
            'flight_cut': 25.0,
            **_default_purity_scan_ranges(),
        },
    },

    'region_definitions': _default_region_definitions(),
    'hist_defs': dict(COMMON_HIST_DEFS),

    # Stage 3 PID selector optimization (Punzi FOM): scans the KM-family
    # ladder (pid_selector.KM_LADDER) for the Lambda0 proton and pion
    # daughters. None until the Stage 3 checkpoint (Matt reviews and
    # applies, same pattern as the Stage 2 cuts above).
    'pid': {
        'lambda0_selector': {'p': None, 'pi': None},
    },
}

# ------------------------------------------------------------------------
# B+ -> LambdaC+ Lambda0
# ------------------------------------------------------------------------
CHANNELS['Lam0LamC'] = {
    'name': 'Lam0LamC',
    'decay_label': r'$B^+ \to \Lambda_c^+ \Lambda^0$',
    'B_Lund': 521,

    'B_daughters': [('LambdaC', 'Bd1Idx'), ('Lambda0', 'Bd2Idx')],

    'composites': {
        'LambdaC': {
            'n_required': 1,
            'mass_var': 'LambdaC_unc_Mass',
            'mass_pdg': LAMBDAC_MASS_PDG,
            # Stage 2 (run_stage02.py): Gaussian+linear-background fit to
            # each mode's signal-MC mass peak (K_S0 gate applied to modes
            # 2/3), mass_window_nsigma * sigma around the fitted mean.
            # Applied 2026-07-24 -- see results/Lam0LamC.yaml (stage02) and
            # the Stage 2 BAD section for the fit plots/justification.
            'mass_windows_per_mode': {
                1: _lambdac_mode_window(2.284779, 0.004209, nsigma=3.0),
                2: _lambdac_mode_window(2.284868, 0.004958, nsigma=3.0),
                3: _lambdac_mode_window(2.284978, 0.004265, nsigma=3.0),
                4: _lambdac_mode_window(2.285011, 0.005228, nsigma=3.0),
            },
            # Width of the per-mode windows, in units of the fitted mode
            # resolution sigma; see the Stage 2 BAD section for the
            # efficiency/purity justification of this choice.
            'mass_window_nsigma': 3.0,
            # Generic-composite fallback (Stage 1 snapshot / any code that
            # still expects a single flat window): envelope of the per-mode
            # windows above. The mode-aware mask in cutflow.py never uses
            # this; it is derived, not authoritative.
            'mass_window': None,  # filled in below, after the dict literal
            # LambdaC flies ~60 um -- no flight-significance cut planned;
            # decided in Stage 2 (note: no LambdaCpostFitFlightSignificance
            # branch exists, only LambdaCFlightLen/LambdaCFlightErr)
            'flight_var': None,
            'flight_cut': None,
        },
        'Lambda0': {
            'n_required': 1,
            'mass_var': 'Lambda0_unc_Mass',
            'mass_pdg': LAMBDA0_MASS_PDG,
            'mass_window': [LAMBDA0_MASS_PDG - 0.003, LAMBDA0_MASS_PDG + 0.003],
            'flight_var': 'Lambda0postFitFlightSignificance',
            # Stage 2 S/sqrt(S+B) scan gave a clean interior optimum at 18
            # (close to the p-Lambda0 default of 25; FOM 321 vs. 273 at 25)
            # -- applied 2026-07-24. The mass-window scan for this composite
            # hit the edge of its scanned range (see STATUS.md open issue),
            # so that value is NOT applied; mass_window stays at the
            # p-Lambda0 default above pending that fix.
            'flight_cut': 18.0,
            **_default_purity_scan_ranges(),
        },
    },

    # K_S0 (daughter of LambdaC in modes 2 and 3 only). Not a B_daughter, so
    # it lives outside 'composites'; the mode-conditional link back to the
    # LambdaC candidate that owns it is built in cutflow.get_lambdac_ks_info.
    'k0s': {
        'mass_var': 'K_SpreFitMass',
        'mass_pdg': K0S_MASS_PDG,
        # Stage 2 S/sqrt(S+B) scan gave a clean interior optimum at +/-26
        # MeV (consistent with the ~7.8 MeV pre-fit resolution noted in
        # Stage 1) -- applied 2026-07-24.
        'mass_window': [K0S_MASS_PDG - 0.026, K0S_MASS_PDG + 0.026],
        'flight_var': 'K_SpostFitFlightSignificance',
        # The flight-significance scan hit the edge of its scanned range
        # (see STATUS.md open issue) -- NOT applied; stays at the Stage 1
        # placeholder pending that fix.
        'flight_cut': 5.0,
        'flight_scan': {'lo': 0.0, 'hi': 100.0, 'step': 2.0},
        # Wider than the shared default: K_SpreFitMass resolution is ~7.8
        # MeV (STATUS.md), so 3 sigma ~ 23 MeV -- the scan needs to cover
        # past that to see the FOM turn over.
        'mass_halfwidth_scan': {'lo': 0.002, 'hi': 0.030, 'step': 0.001},
    },

    'region_definitions': _default_region_definitions(),
    'hist_defs': dict(COMMON_HIST_DEFS),
}

# Derived fallback: LambdaC's flat 'mass_window' is the envelope of its
# per-mode windows (generic consumers only; see comment above).
_lamc_mode_windows = CHANNELS['Lam0LamC']['composites']['LambdaC']['mass_windows_per_mode']
CHANNELS['Lam0LamC']['composites']['LambdaC']['mass_window'] = [
    min(w[0] for w in _lamc_mode_windows.values()),
    max(w[1] for w in _lamc_mode_windows.values()),
]
del _lamc_mode_windows

# LambdaC-specific histograms
CHANNELS['Lam0LamC']['hist_defs']['LambdaC_unc_Mass'] = \
    {'nbins': 100, 'lo': 2.2, 'hi': 2.38, 'label': r'$\Lambda_c^+$ mass [GeV/c$^2$]'}
CHANNELS['Lam0LamC']['hist_defs']['LambdaCFlightLen'] = \
    {'nbins': 100, 'lo': -1, 'hi': 10, 'label': r'$\Lambda_c^+$ flight length [cm]'}
CHANNELS['Lam0LamC']['hist_defs']['nLambdaC'] = \
    {'nbins': 30, 'lo': 0, 'hi': 30, 'label': r'# of $\Lambda_c^+$ candidates'}

# Candidate multiplicities run much higher in this channel (up to ~60 B
# candidates/event seen) -- widen the ranges relative to the common defaults
CHANNELS['Lam0LamC']['hist_defs']['nB'] = \
    {'nbins': 30, 'lo': 0, 'hi': 30, 'label': '# of B candidates'}
CHANNELS['Lam0LamC']['hist_defs']['nLambda0'] = \
    {'nbins': 20, 'lo': 0, 'hi': 20, 'label': r'# of $\Lambda^0$ candidates'}

# K_S0 (daughter of LambdaC in modes 2 and 3)
# NOTE: use the pre-fit mass -- K_SMass is post-fit (mass-constrained) and
# is just a spike at the PDG value
CHANNELS['Lam0LamC']['hist_defs']['K_SpreFitMass'] = \
    {'nbins': 100, 'lo': 0.46, 'hi': 0.54, 'label': r'$K_S^0$ mass (pre-fit) [GeV/c$^2$]'}
CHANNELS['Lam0LamC']['hist_defs']['K_SpostFitFlightSignificance'] = \
    {'nbins': 100, 'lo': -20, 'hi': 300, 'label': r'$K_S^0$ flight significance'}

# Derived fields (added by datasets.add_derived_fields)
CHANNELS['Lam0LamC']['hist_defs']['LambdaCFlightSignificance'] = \
    {'nbins': 100, 'lo': -10, 'hi': 20, 'label': r'$\Lambda_c^+$ flight significance (FlightLen/FlightErr)'}
CHANNELS['Lam0LamC']['hist_defs']['Lambda0FLSig_fromB'] = \
    {'nbins': 100, 'lo': -20, 'hi': 300, 'label': r'$\Lambda^0$ (from $B$) flight significance'}
CHANNELS['Lam0LamC']['hist_defs']['Lambda0FLSig_fromLambdaC'] = \
    {'nbins': 100, 'lo': -20, 'hi': 300, 'label': r'$\Lambda^0$ (from $\Lambda_c^+$) flight significance'}

# LambdaC reconstruction/decay modes (classified from LambdaCnDaus and
# LambdaCd1Lund by cutflow.get_lambdac_decay_mode). Signal MC was generated
# with equal amounts of each.
CHANNELS['Lam0LamC']['lambdac_modes'] = {
    1: r'$p K^- \pi^+$',
    2: r'$p K_S^0$',
    3: r'$p K_S^0 \pi^+ \pi^-$',
    4: r'$\Lambda^0 \pi^+ \pi^+ \pi^-$',
}

# Stage 3 PID selector optimization (Punzi FOM): scans the KM-family ladder
# (pid_selector.KM_LADDER) per LambdaC decay mode -- the daughter slots
# scanned differ by mode (verified empirically against LambdaCdNLund on
# signal MC; see STATUS.md Stage 3 notes):
#   mode 1 (p K- pi+):        p (d1) + K (d2) + pi (d3)
#   mode 2 (p K_S0):           p (d1) only -- K_S0's own daughters are not a
#                              LambdaC daughter slot, no PID cut planned
#   mode 3 (p K_S0 pi+ pi-):   p (d1) + pi (d3, d4 -- same selector, both)
#   mode 4 (Lam0 pi+ pi+ pi-): pi (d2, d3, d4 -- same selector, all three);
#                              the inner Lambda0's own p/pi reuse
#                              'lambda0_selector' below, not scanned again.
# None until the Stage 3 checkpoint (Matt reviews and applies).
CHANNELS['Lam0LamC']['pid'] = {
    'lambda0_selector': {'p': None, 'pi': None},
    'lambdac_selector_per_mode': {
        1: {'p': None, 'K': None, 'pi': None},
        2: {'p': None},
        3: {'p': None, 'pi': None},
        4: {'pi': None},
    },
}


################################################################################
def get_channel_config(channel):
    """
    Return a (deep) copy of the configuration dictionary for a channel.

    channel: 'Lam0Lam0' or 'Lam0LamC'
    """
    if channel not in CHANNELS:
        raise ValueError(f"Unknown channel '{channel}'. Choose from {list(CHANNELS.keys())}")

    return copy.deepcopy(CHANNELS[channel])
################################################################################
