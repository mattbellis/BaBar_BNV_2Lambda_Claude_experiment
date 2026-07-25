"""
PID selector-bitmap decoding for the BNV 2-Lambda analyses.

Ported from BNV_pLambda/myPIDselector.py (reference, do not modify) --
channel-agnostic, so it lives next to datasets.py/cutflow.py rather than
being duplicated per channel.

Each charged-track collection (p/K/pi/mu/e) carries a per-candidate integer
bitmap (<particle>SelectorsMap) where bit i (0 = least significant) records
whether that candidate passes the i-th named PID selector for that particle
hypothesis (see the BaBar PID Selector List wiki page). SELECTORS[particle]
gives the bit -> name mapping, one list per particle, in the same order as
the reference's PIDselector class (validated against it on real signal MC:
same selector, same candidates, identical boolean mask).

<particle>TrkIdx (e.g. 'pTrkIdx') maps each hypothesis-collection candidate
(e.g. the i-th entry in the 'p' collection) to its index in the global
<particle>SelectorsMap array; a composite's daughter index (e.g.
Lambda0d1Idx) then indexes into that hypothesis collection. See
bits_for_hypothesis / passes_selector below.
"""

SELECTORS = {}

SELECTORS['p'] = [
    'VeryLooseLHProtonSelection', 'LooseLHProtonSelection',
    'TightLHProtonSelection', 'VeryTightLHProtonSelection',

    'VeryLooseGLHProtonSelection', 'LooseGLHProtonSelection',
    'TightGLHProtonSelection', 'VeryTightGLHProtonSelection',

    'VeryLooseELHProtonSelection', 'LooseELHProtonSelection',
    'TightELHProtonSelection', 'VeryTightELHProtonSelection',

    'SuperLooseKMProtonSelection', 'VeryLooseKMProtonSelection',
    'LooseKMProtonSelection', 'TightKMProtonSelection',
    'VeryTightKMProtonSelection', 'SuperTightKMProtonSelection',
]

SELECTORS['K'] = [
    'NotPionKaonMicroSelection', 'VeryLooseKaonMicroSelection',
    'LooseKaonMicroSelection', 'TightKaonMicroSelection', 'VeryTightKaonMicroSelection',

    'NotPionNNKaonMicroSelection', 'VeryLooseNNKaonMicroSelection',
    'LooseNNKaonMicroSelection', 'TightNNKaonMicroSelection', 'VeryTightNNKaonMicroSelection',

    'NotPionLHKaonMicroSelection', 'VeryLooseLHKaonMicroSelection',
    'LooseLHKaonMicroSelection', 'TightLHKaonMicroSelection', 'VeryTightLHKaonMicroSelection',

    'VeryLooseGLHKaonMicroSelection', 'LooseGLHKaonMicroSelection',
    'TightGLHKaonMicroSelection', 'VeryTightGLHKaonMicroSelection',

    'NotPionBDTKaonMicroSelection', 'VeryLooseBDTKaonMicroSelection',
    'LooseBDTKaonMicroSelection', 'TightBDTKaonMicroSelection', 'VeryTightBDTKaonMicroSelection',

    'SuperLooseKMKaonMicroSelection', 'VeryLooseKMKaonMicroSelection',
    'LooseKMKaonMicroSelection', 'TightKMKaonMicroSelection',
    'VeryTightKMKaonMicroSelection', 'SuperTightKMKaonMicroSelection',
]

SELECTORS['pi'] = [
    'PidRoyPionSelectionLoose', 'PidRoyPionSelectionNotKaon',

    'VeryLooseLHPionSelection', 'LooseLHPionSelection',
    'TightLHPionSelection', 'VeryTightLHPionSelection',

    'VeryLooseGLHPionSelection', 'LooseGLHPionSelection',
    'TightGLHPionSelection', 'VeryTightGLHPionSelection',

    'SuperLooseKMPionMicroSelection', 'VeryLooseKMPionMicroSelection',
    'LooseKMPionMicroSelection', 'TightKMPionMicroSelection',
    'VeryTightKMPionMicroSelection', 'SuperTightKMPionMicroSelection',
]

SELECTORS['e'] = [
    'NoCalElectronMicroSelection', 'VeryLooseElectronMicroSelection',
    'LooseElectronMicroSelection', 'TightElectronMicroSelection',
    'VeryTightElectronMicroSelection',

    'PidLHElectrons',

    'SuperLooseKMElectronMicroSelection', 'VeryLooseKMElectronMicroSelection',
    'LooseKMElectronMicroSelection', 'TightKMElectronMicroSelection',
    'VeryTightKMElectronMicroSelection', 'SuperTightKMElectronMicroSelection',
]

SELECTORS['mu'] = [
    'MinimumIoniziongMuonMicroSelection', 'VeryLooseMuonMicroSelection',
    'LooseMuonMicroSelection', 'TightMuonMicroSelection', 'VeryTightMuonMicroSelection',

    'NNVeryLooseMuonSelection', 'NNLooseMuonSelection',
    'NNTightMuonSelection', 'NNVeryTightMuonSelection',
    'NNVeryLooseMuonSelectionFakeRate', 'NNLooseMuonSelectionFakeRate',
    'NNTightMuonSelectionFakeRate', 'NNVeryTightMuonSelectionFakeRate',

    'LikeVeryLooseMuonSelection', 'LikeLooseMuonSelection', 'LikeTightMuonSelection',

    'BDTVeryLooseMuonSelection', 'BDTLooseMuonSelection',
    'BDTTightMuonSelection', 'BDTVeryTightMuonSelection',
    'BDTVeryLooseMuonSelectionFakeRate', 'BDTLooseMuonSelectionFakeRate',
    'BDTTightMuonSelectionFakeRate', 'BDTVeryTightMuonSelectionFakeRate',

    'BDTLoPLooseMuonSelection', 'BDTLoPTightMuonSelection',
]

# The KM-family ladder for each particle type, loosest -> tightest -- the
# ladder scanned by Stage 3's Punzi-FOM optimization (see pid_optimization.py
# and the Stage 3 plan in STATUS.md). Kept here, next to SELECTORS, since
# it's a property of the bit -> name mapping above, not of any one channel.
KM_LADDER = {
    'p': [
        'SuperLooseKMProtonSelection', 'VeryLooseKMProtonSelection', 'LooseKMProtonSelection',
        'TightKMProtonSelection', 'VeryTightKMProtonSelection', 'SuperTightKMProtonSelection',
    ],
    'K': [
        'SuperLooseKMKaonMicroSelection', 'VeryLooseKMKaonMicroSelection', 'LooseKMKaonMicroSelection',
        'TightKMKaonMicroSelection', 'VeryTightKMKaonMicroSelection', 'SuperTightKMKaonMicroSelection',
    ],
    'pi': [
        'SuperLooseKMPionMicroSelection', 'VeryLooseKMPionMicroSelection', 'LooseKMPionMicroSelection',
        'TightKMPionMicroSelection', 'VeryTightKMPionMicroSelection', 'SuperTightKMPionMicroSelection',
    ],
}

################################################################################
def bit_index(particle, selector_name):
    """Bit position (0 = least significant) of a named selector for one particle hypothesis."""
    return SELECTORS[particle].index(selector_name)

################################################################################
def bits_for_hypothesis(trkidx, selectors_map):
    """
    Per-candidate integer PID bitmap for one hypothesis collection (e.g. the
    'p' collection): selectors_map[trkidx], aligned with the hypothesis
    collection (same shape as trkidx), not the global track collection
    (selectors_map's own indexing).
    """
    return selectors_map[trkidx]

################################################################################
def passes_selector(bits, selector_name, particle):
    """
    Jagged boolean mask: whether `selector_name` is set in each entry of
    `bits` (as returned by bits_for_hypothesis / indexed further by a
    composite's daughter-index branch), for the given particle type.
    """
    i = bit_index(particle, selector_name)
    return (bits >> i) % 2 == 1
################################################################################
