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
            'mass_window': [LAMBDA0_MASS_PDG - 0.003, LAMBDA0_MASS_PDG + 0.003],
            'flight_var': 'Lambda0postFitFlightSignificance',
            'flight_cut': 25.0,  # from p-Lambda0; to be re-optimized in Stage 2
        },
    },

    'region_definitions': _default_region_definitions(),
    'hist_defs': dict(COMMON_HIST_DEFS),
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
            # Placeholder window; to be set from purity studies in Stage 2
            'mass_window': [LAMBDAC_MASS_PDG - 0.010, LAMBDAC_MASS_PDG + 0.010],
            # LambdaC flies ~60 um -- no flight-significance cut foreseen;
            # decide in Stage 2 (note: no LambdaCpostFitFlightSignificance
            # branch exists, only LambdaCFlightLen/LambdaCFlightErr)
            'flight_var': None,
            'flight_cut': None,
        },
        'Lambda0': {
            'n_required': 1,
            'mass_var': 'Lambda0_unc_Mass',
            'mass_window': [LAMBDA0_MASS_PDG - 0.003, LAMBDA0_MASS_PDG + 0.003],
            'flight_var': 'Lambda0postFitFlightSignificance',
            'flight_cut': 25.0,  # from p-Lambda0; to be re-optimized in Stage 2
        },
    },

    'region_definitions': _default_region_definitions(),
    'hist_defs': dict(COMMON_HIST_DEFS),
}

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
