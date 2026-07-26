"""
Cut/mask-building functions for the BNV 2-Lambda analyses.

Conventions (following the p-Lambda0 reference analysis):

- "candidate" masks are jagged, one boolean per candidate, and can be used to
  slice the candidate-level (jagged) fields.
- "event" masks are flat, one boolean per event, and can be used to slice the
  whole awkward array.
- A set of cuts is collected in a dict-of-dicts, `dcuts`, where each entry has
    dcuts[n]['name']        human-readable name
    dcuts[n]['event']       event-level mask (same length as the full array)
    dcuts[n]['candidates']  candidate-level mask or None
  Cut n=-1 is the AND of everything.

All functions take the channel config dict from channel_config.py so that the
same code serves both B0 -> Lambda0 Lambda0 and B+ -> LambdaC+ Lambda0.
"""

import awkward as ak
import numpy as np

import pid_selector

################################################################################
# Utilities
################################################################################
def munge_mask_shapes(mask_larger, mask_smaller):
    """
    Expand `mask_smaller` (defined on the subset of events where `mask_larger`
    is True) back to the full length of `mask_larger`, filling False elsewhere.
    """
    idx = ak.local_index(mask_larger)

    mask = ak.zeros_like(mask_larger, dtype=bool)
    mask = mask.to_numpy()
    mask[idx[mask_larger]] = mask_smaller

    return mask

################################################################################
def indices_to_booleans(indices, array_to_slice):
    """
    Convert jagged integer indices into a jagged boolean mask over
    `array_to_slice` (True where a candidate's index appears in `indices`).
    """
    whole_set, in_set = ak.unzip(ak.cartesian([
        ak.local_index(array_to_slice), indices], nested=True))

    return ak.any(whole_set == in_set, axis=-1)

################################################################################
# LambdaC decay-mode classification (Lam0LamC channel)
################################################################################
def get_lambdac_decay_mode(data):
    """
    Classify each LambdaC candidate by its reconstruction mode, using the
    number of daughters and the Lund ID of the first daughter:

        1: p K- pi+           (nDaus == 3)
        2: p K_S0             (nDaus == 2)
        3: p K_S0 pi+ pi-     (nDaus == 4, |d1| == 2212)
        4: Lam0 pi+ pi+ pi-   (nDaus == 4, |d1| == 3122)
        0: anything else (should not happen)

    Returns a jagged integer array, one entry per LambdaC candidate.
    """
    ndaus = data['LambdaCnDaus']
    d1 = abs(data['LambdaCd1Lund'])

    mode = ak.zeros_like(ndaus)
    mode = ak.where(ndaus == 3, 1, mode)
    mode = ak.where(ndaus == 2, 2, mode)
    mode = ak.where((ndaus == 4) & (d1 == 2212), 3, mode)
    mode = ak.where((ndaus == 4) & (d1 == 3122), 4, mode)

    return mode

################################################################################
def get_lambdac_decay_mode_per_B(data, config):
    """
    The LambdaC decay mode of each *B* candidate's LambdaC daughter
    (jagged, one entry per B candidate, aligned with e.g. BpostFitMes).
    """
    mode = get_lambdac_decay_mode(data)

    # Which B-daughter index points into the LambdaC collection?
    idxvar = None
    for name, iv in config['B_daughters']:
        if name == 'LambdaC':
            idxvar = iv
    if idxvar is None:
        raise ValueError(f"Channel {config['name']} has no LambdaC among the B daughters")

    return mode[data[idxvar]]

################################################################################
# Event-level candidate counting (single-candidate requirement)
################################################################################
def get_single_candidate_mask(data, config):
    """
    Event mask requiring exactly 1 B candidate and exactly the expected
    number of each composite candidate for this channel:

        Lam0Lam0:  nB == 1 and nLambda0 == 2
        Lam0LamC:  nB == 1 and nLambda0 == 1 and nLambdaC == 1
    """
    mask = ak.num(data['BpostFitMes']) == 1

    for name, comp in config['composites'].items():
        n = ak.num(data[comp['mass_var']])
        mask = mask & (n == comp['n_required'])

    return mask

################################################################################
# mES / DeltaE region masks (candidate-level)
################################################################################
def get_fit_mask(data, config):
    """Candidate mask for the mES/DeltaE fitting region."""
    region_definitions = config['region_definitions']

    mes = data['BpostFitMes']
    de = data['BpostFitDeltaE']

    fit_mask = (mes > region_definitions['fitting MES'][0]) & \
               (mes < region_definitions['fitting MES'][1]) & \
               (de > region_definitions['fitting DeltaE'][0]) & \
               (de < region_definitions['fitting DeltaE'][1])

    return fit_mask

################################################################################
def get_signal_region_mask(data, config):
    """
    Candidate mask for the mES/DeltaE signal region.

    For MC this selects the signal box; for collision data this is the
    BLINDED region and must only ever be used to *verify* that the region is
    empty (or, after human-approved unblinding, in the inference stage).
    """
    region_definitions = config['region_definitions']

    mes = data['BpostFitMes']
    de = data['BpostFitDeltaE']

    signal_mask = (mes > region_definitions['signal MES'][0]) & \
                  (mes < region_definitions['signal MES'][1]) & \
                  (de > region_definitions['signal DeltaE'][0]) & \
                  (de < region_definitions['signal DeltaE'][1])

    return signal_mask

################################################################################
# Composite (Lambda0 / LambdaC / K_S0) purity masks
################################################################################
def get_purity_mask_for_comp(data, comp):
    """
    Candidate mask for one composite-config dict {mass_var, mass_window,
    flight_var, flight_cut}: mass window AND flight-significance cut (where
    configured). Takes the raw comp dict rather than a channel + name so it
    also serves 'k0s', which lives outside config['composites'].
    """
    lo, hi = comp['mass_window']
    m = data[comp['mass_var']]
    mask = (m > lo) & (m < hi)

    if comp['flight_var'] is not None:
        mask = mask & (data[comp['flight_var']] > comp['flight_cut'])

    return mask

################################################################################
def get_composite_purity_mask(data, config, name):
    """
    Candidate mask for one composite species by name ('Lambda0' or
    'LambdaC'). LambdaC is mode-dependent (see get_lambdac_purity_mask) and
    is dispatched there automatically; Lambda0 uses the generic mass+flight
    mask.
    """
    if name == 'LambdaC' and 'mass_windows_per_mode' in config['composites'][name]:
        return get_lambdac_purity_mask(data, config)

    return get_purity_mask_for_comp(data, config['composites'][name])

################################################################################
def get_all_composite_purity_masks(data, config):
    """
    Event and candidate masks for all composite species in this channel.

    Returns (mask_event, candidate_masks) where candidate_masks is a dict
    {name: jagged mask} and mask_event requires the number of candidates
    passing the purity cuts to equal the number required for the channel.
    """
    candidate_masks = {}
    mask_event = None

    for name, comp in config['composites'].items():
        cmask = get_composite_purity_mask(data, config, name)
        candidate_masks[name] = cmask

        npass = ak.num(data[comp['mass_var']][cmask])
        ev = npass == comp['n_required']

        mask_event = ev if mask_event is None else (mask_event & ev)

    return mask_event, candidate_masks

################################################################################
# K_S0 <-> LambdaC linkage and mode-dependent LambdaC purity (Lam0LamC only)
################################################################################
def get_lambdac_ks_info(data):
    """
    For each LambdaC candidate, find its K_S0 daughter (present only in
    modes 2 and 3): checks all four daughter slots for |Lund| == 310 (the
    K_S0 code), mirroring how the mode-4 Lambda0 daughter is picked out via
    LambdaCd1Lund in datasets.add_derived_fields.

    Returns (has_ks, ks_idx): jagged bool / int arrays aligned with the
    LambdaC collection. ks_idx is only meaningful where has_ks is True (it
    defaults to 0, a placeholder index, elsewhere).
    """
    K_S_LUND = 310

    has_ks = ak.zeros_like(data['LambdaCnDaus'], dtype=bool)
    ks_idx = ak.zeros_like(data['LambdaCd1Idx'])

    for i in (1, 2, 3, 4):
        d_lund = abs(data[f'LambdaCd{i}Lund'])
        d_idx = data[f'LambdaCd{i}Idx']

        sel = (d_lund == K_S_LUND)
        has_ks = has_ks | sel
        ks_idx = ak.where(sel, d_idx, ks_idx)

    return has_ks, ks_idx

################################################################################
def get_lambdac_k0s_gate(data, config):
    """
    Jagged mask aligned with the LambdaC collection: True if this LambdaC
    candidate has no K_S0 daughter (modes 1, 4 -- the cut does not apply),
    or has one that passes the K_S0 purity cut (modes 2, 3).
    """
    has_ks, ks_idx = get_lambdac_ks_info(data)

    # ks_idx defaults to placeholder 0 where has_ks is False; events with
    # zero K_S candidates at all would make ks_purity[ks_idx] an out-of-range
    # gather in that case. Pad every event's K_S purity list to length >= 1
    # (with False) so index 0 always exists; its value is never used where
    # has_ks is False (short-circuited by the `~has_ks |` below).
    ks_purity = get_purity_mask_for_comp(data, config['k0s'])
    ks_purity = ak.fill_none(ak.pad_none(ks_purity, 1, axis=1), False)

    return (~has_ks) | (has_ks & ks_purity[ks_idx])

################################################################################
def get_lambdac_mode_mass_window_mask(data, config):
    """
    Jagged mask aligned with the LambdaC collection: mass within the
    per-mode window from config['composites']['LambdaC']['mass_windows_per_mode'].
    """
    mode = get_lambdac_decay_mode(data)
    mass = data[config['composites']['LambdaC']['mass_var']]
    windows = config['composites']['LambdaC']['mass_windows_per_mode']

    mask = ak.zeros_like(mode, dtype=bool)
    for m, (lo, hi) in windows.items():
        mask = mask | ((mode == m) & (mass > lo) & (mass < hi))

    return mask

################################################################################
def get_lambdac_purity_mask(data, config):
    """
    Full Stage 2 LambdaC purity mask (jagged, aligned with the LambdaC
    collection): per-mode mass window AND the K_S0 gate (sequential
    optimization -- the K_S0 cut is fixed first, then the LambdaC windows
    are set with it applied; see run_stage02.py).
    """
    return get_lambdac_mode_mass_window_mask(data, config) & get_lambdac_k0s_gate(data, config)

################################################################################
# Per-B purity mask (all composite daughters pass their purity cuts)
################################################################################
def get_composite_purity_masks_per_B(data, config, candidate_masks=None):
    """
    Jagged mask aligned with the B collection: True if every composite
    daughter of that B candidate passes its purity mask. This is the
    per-candidate building block later stages (PID, antibaryon veto, MLP)
    reuse; get_all_composite_purity_masks's event-level count is a coarser
    summary that does not check the daughters belong to the *same* B.

    candidate_masks: optional pre-computed {name: jagged mask} (e.g. from
    get_all_composite_purity_masks), to avoid recomputing.
    """
    if candidate_masks is None:
        _, candidate_masks = get_all_composite_purity_masks(data, config)

    mask_b = None
    for name, idxvar in config['B_daughters']:
        dmask = candidate_masks[name][data[idxvar]]
        mask_b = dmask if mask_b is None else (mask_b & dmask)

    return mask_b

################################################################################
# Stage 3: PID selector masks
#
# selector_name=None always means "skip this cut" (pass-all) so callers can
# scan one particle at a time, or apply a partial config where some slots
# are still unset (see channel_config.py's 'pid' section).
################################################################################
def _safe_gather(bool_array, idx, valid_mask):
    """
    Gather a jagged boolean array (aligned to some hypothesis/composite
    collection) at per-candidate index `idx`, restricted to rows where
    `valid_mask` is True. Elsewhere `idx` may not even apply to this
    collection (e.g. a different LambdaC decay mode's daughter slot) or may
    be an unused-slot placeholder (0) in an event whose collection is empty
    for this particle type -- both would otherwise raise an out-of-range/
    empty-list error on the gather. Mirrors the ks_idx placeholder + pad_none
    pattern in get_lambdac_k0s_gate; result is only meaningful where
    valid_mask is True.
    """
    safe_idx = ak.where(valid_mask, idx, 0)
    padded = ak.fill_none(ak.pad_none(bool_array, 1, axis=1), False)
    return padded[safe_idx]

################################################################################
def get_lambda0_pid_mask(data, p_selector=None, pi_selector=None):
    """
    Candidate mask for the Lambda0 collection (both channels; also reused
    for Lam0LamC LambdaC mode 4's inner Lambda0, via LambdaCd1Idx): proton
    daughter (Lambda0d1Idx, into the 'p' hypothesis collection) and pion
    daughter (Lambda0d2Idx, into 'pi') PID selector cuts. No empty-collection
    guard needed here: a Lambda0 candidate cannot exist in an event without
    its own p/pi daughters already present in those collections.
    """
    mask = ak.ones_like(data['Lambda0d1Idx'], dtype=bool)

    if p_selector is not None:
        p_bits = pid_selector.bits_for_hypothesis(data['pTrkIdx'], data['pSelectorsMap'])
        p_pass = pid_selector.passes_selector(p_bits, p_selector, 'p')
        mask = mask & p_pass[data['Lambda0d1Idx']]

    if pi_selector is not None:
        pi_bits = pid_selector.bits_for_hypothesis(data['piTrkIdx'], data['piSelectorsMap'])
        pi_pass = pid_selector.passes_selector(pi_bits, pi_selector, 'pi')
        mask = mask & pi_pass[data['Lambda0d2Idx']]

    return mask

################################################################################
def get_lambdac_pid_mask(data, selectors_per_mode, lambda0_pid_mask):
    """
    Candidate mask for the LambdaC collection (Lam0LamC only): mode-dependent
    PID selector cuts on the fixed daughter slots per mode (empirically
    verified against LambdaCdNLund on signal MC -- see STATUS.md Stage 3
    notes):

        mode 1 (p K- pi+):        d1=p, d2=K, d3=pi
        mode 2 (p K_S0):           d1=p
        mode 3 (p K_S0 pi+ pi-):   d1=p, d3=pi, d4=pi (same selector, both)
        mode 4 (Lam0 pi+ pi+ pi-): d2=pi, d3=pi, d4=pi (same selector, all
                                   three); the inner Lambda0 (d1) reuses
                                   lambda0_pid_mask via LambdaCd1Idx rather
                                   than being re-decoded here.

    selectors_per_mode: {mode: {particle: selector_name_or_None}} (see
    channel_config.py's 'pid'/'lambdac_selector_per_mode'); a missing
    particle key, or a None value, skips that cut.
    lambda0_pid_mask: jagged mask aligned with the Lambda0 collection (e.g.
    from get_lambda0_pid_mask), reused for mode 4's inner Lambda0.
    """
    mode = get_lambdac_decay_mode(data)

    p_bits = pid_selector.bits_for_hypothesis(data['pTrkIdx'], data['pSelectorsMap'])
    k_bits = pid_selector.bits_for_hypothesis(data['KTrkIdx'], data['KSelectorsMap'])
    pi_bits = pid_selector.bits_for_hypothesis(data['piTrkIdx'], data['piSelectorsMap'])

    def sel_for(m, particle):
        return selectors_per_mode.get(m, {}).get(particle)

    mask = ak.ones_like(mode, dtype=bool)

    # Mode 1: p (d1), K (d2), pi (d3)
    m1 = (mode == 1)
    s = sel_for(1, 'p')
    if s is not None:
        p_pass = pid_selector.passes_selector(p_bits, s, 'p')
        mask = mask & (~m1 | _safe_gather(p_pass, data['LambdaCd1Idx'], m1))
    s = sel_for(1, 'K')
    if s is not None:
        k_pass = pid_selector.passes_selector(k_bits, s, 'K')
        mask = mask & (~m1 | _safe_gather(k_pass, data['LambdaCd2Idx'], m1))
    s = sel_for(1, 'pi')
    if s is not None:
        pi_pass = pid_selector.passes_selector(pi_bits, s, 'pi')
        mask = mask & (~m1 | _safe_gather(pi_pass, data['LambdaCd3Idx'], m1))

    # Mode 2: p (d1) only
    m2 = (mode == 2)
    s = sel_for(2, 'p')
    if s is not None:
        p_pass = pid_selector.passes_selector(p_bits, s, 'p')
        mask = mask & (~m2 | _safe_gather(p_pass, data['LambdaCd1Idx'], m2))

    # Mode 3: p (d1), pi (d3 and d4, same selector)
    m3 = (mode == 3)
    s = sel_for(3, 'p')
    if s is not None:
        p_pass = pid_selector.passes_selector(p_bits, s, 'p')
        mask = mask & (~m3 | _safe_gather(p_pass, data['LambdaCd1Idx'], m3))
    s = sel_for(3, 'pi')
    if s is not None:
        pi_pass = pid_selector.passes_selector(pi_bits, s, 'pi')
        pi3_ok = (_safe_gather(pi_pass, data['LambdaCd3Idx'], m3) &
                  _safe_gather(pi_pass, data['LambdaCd4Idx'], m3))
        mask = mask & (~m3 | pi3_ok)

    # Mode 4: pi (d2, d3, d4, same selector); inner Lambda0 (d1) via
    # lambda0_pid_mask, not re-decoded
    m4 = (mode == 4)
    s = sel_for(4, 'pi')
    if s is not None:
        pi_pass = pid_selector.passes_selector(pi_bits, s, 'pi')
        pi4_ok = (_safe_gather(pi_pass, data['LambdaCd2Idx'], m4) &
                  _safe_gather(pi_pass, data['LambdaCd3Idx'], m4) &
                  _safe_gather(pi_pass, data['LambdaCd4Idx'], m4))
        mask = mask & (~m4 | pi4_ok)
    inner_lambda0_ok = _safe_gather(lambda0_pid_mask, data['LambdaCd1Idx'], m4)
    mask = mask & (~m4 | inner_lambda0_ok)

    return mask

################################################################################
def get_pid_candidate_masks(data, config):
    """
    The Stage 3 PID candidate masks currently configured for this channel, as
    a {composite_name: jagged mask} dict ready to pass to
    get_composite_purity_masks_per_B. Any selector left None in the config's
    'pid' section is skipped (pass-all), so this is safe to call before
    Stage 3 has been accepted.
    """
    pid_cfg = config.get('pid', {})
    lambda0_selector = pid_cfg.get('lambda0_selector', {}) or {}

    lambda0_pid = get_lambda0_pid_mask(data,
                                       p_selector=lambda0_selector.get('p'),
                                       pi_selector=lambda0_selector.get('pi'))
    masks = {'Lambda0': lambda0_pid}

    if 'LambdaC' in config['composites']:
        masks['LambdaC'] = get_lambdac_pid_mask(
            data, pid_cfg.get('lambdac_selector_per_mode', {}) or {}, lambda0_pid)

    return masks

################################################################################
def get_pid_mask_per_B(data, config):
    """Per-B mask: every composite daughter passes its Stage 3 PID cut."""
    return get_composite_purity_masks_per_B(data, config, get_pid_candidate_masks(data, config))

################################################################################
# Stage 4: antibaryon (antiproton) veto
#
# Physics: the BNV signal final state contains baryons only -- two protons
# (or, for the charge-conjugate B, two antiprotons) and no antibaryon. SM
# decays that produce a baryon almost always produce a compensating
# antibaryon. So an event containing an identified antiproton that is not
# part of the signal B candidate is background-like.
#
# Structure of the ntuple, verified empirically on signal+background MC
# (103658 Lam0Lam0 / 286890 Lam0LamC events; NOT assumed -- see the Stage 4
# BAD section and STATUS.md):
#
#  - '<particle>SelectorsMap' is indexed by TRACK: len == nTRK, one entry per
#    entry of TRKLund. The '<particle>' hypothesis list is a strict SUBSET of
#    the tracks (np is 2-24 smaller than nTRK in every event), and 24% of the
#    tracks passing SuperLooseKMProtonSelection lie OUTSIDE the p list. So the
#    eligible-track pool is a real choice, not a formality -- see `pool`.
#  - sign(TRKLund) is the electric charge (TRKLund takes only +/-211): it
#    agrees with sign(pLund[pTrkIdx]) in 207562 of 207563 proton entries. The
#    single disagreement is an isolated ntuple inconsistency, not a pattern.
#  - The signal baryon sign is consistent within a B candidate: both Lambda0
#    daughters carry the same sign (Lam0Lam0, 0 mismatches), and
#    sign(LambdaCLund) == sign(the Lambda0's proton daughter) (Lam0LamC, 0
#    mismatches, including mode 4's inner Lambda0). So one uniform rule --
#    read the sign off the B's Lambda0 daughter's proton -- serves both
#    channels and every LambdaC mode, and charge conjugation is handled
#    automatically rather than hardcoded.
################################################################################
PROTON_LUND = 2212

def _safe_gather_int(int_array, idx, valid_mask, fill=-1):
    """
    Integer counterpart of _safe_gather: gather `int_array` at per-candidate
    index `idx` only where `valid_mask` is True, filling `fill` elsewhere.
    Invalid entries (a daughter slot that does not apply to this LambdaC
    decay mode, or a placeholder index into an empty collection) would
    otherwise raise an out-of-range error on the gather.
    """
    safe_idx = ak.where(valid_mask, idx, 0)
    padded = ak.fill_none(ak.pad_none(int_array, 1, axis=1), fill)
    return ak.where(valid_mask, padded[safe_idx], fill)

################################################################################
def get_signal_b_track_slots(data, config, scope='all_signal_tracks'):
    """
    Track indices used by each B candidate: the exclusion set for the
    antibaryon veto ("not part of the signal candidate").

    Returns (slot_idx, slot_valid): python lists of jagged arrays aligned
    with the B collection, one pair per final-state daughter slot.
    slot_idx[k] is the TRK-collection index of that slot's track (-1 where
    the slot does not apply); slot_valid[k] is the corresponding boolean.
    Slots are returned rather than a flat per-B x per-track boolean because
    the latter is O(nB x nTRK) per event -- prohibitive at Lam0LamC's ~60 B
    candidates/event.

    scope:
      'all_signal_tracks'      -- every charged track of the B candidate
                                  (recommended: the Lambda0 pions carry
                                  exactly the antiproton charge and would
                                  otherwise fake the veto)
      'proton_daughters_only'  -- only the proton-hypothesis daughters, as in
                                  the p-Lambda0 reference implementation
                                  (BNV_pLambda/babar_analysis_tools.py's
                                  build_antiproton_antimask), kept so the
                                  difference from the reference is a measured
                                  number rather than an assertion

    The LambdaC daughter-slot -> particle assignment per mode is the same
    empirically-verified map used by get_lambdac_pid_mask (see STATUS.md).
    """
    if scope not in ('all_signal_tracks', 'proton_daughters_only'):
        raise ValueError(f"Unknown scope '{scope}'")

    proton_only = (scope == 'proton_daughters_only')

    p_trk = data['pTrkIdx']
    pi_trk = data['piTrkIdx']

    slot_idx, slot_valid = [], []

    def add(idx, valid, is_proton_hyp):
        if proton_only and not is_proton_hyp:
            return
        slot_idx.append(idx)
        slot_valid.append(valid)

    # Track index of each Lambda0 candidate's proton / pion daughter,
    # aligned with the Lambda0 collection.
    lam_p = p_trk[data['Lambda0d1Idx']]
    lam_pi = pi_trk[data['Lambda0d2Idx']]

    # ---- Lambda0 daughters of the B (both channels; twice for Lam0Lam0) ----
    for name, idxvar in config['B_daughters']:
        if name != 'Lambda0':
            continue
        bidx = data[idxvar]
        always = ak.ones_like(bidx, dtype=bool)
        add(lam_p[bidx], always, True)
        add(lam_pi[bidx], always, False)

    if 'LambdaC' not in config['composites']:
        return slot_idx, slot_valid

    # ---- LambdaC daughters (Lam0LamC), dispatched by decay mode ----
    lamc_idxvar = [iv for (n, iv) in config['B_daughters'] if n == 'LambdaC'][0]
    bidx = data[lamc_idxvar]

    mode = get_lambdac_decay_mode_per_B(data, config)
    d1 = data['LambdaCd1Idx'][bidx]
    d2 = data['LambdaCd2Idx'][bidx]
    d3 = data['LambdaCd3Idx'][bidx]
    d4 = data['LambdaCd4Idx'][bidx]

    m1, m2, m3, m4 = (mode == 1), (mode == 2), (mode == 3), (mode == 4)

    # d1: proton for modes 1/2/3, an inner Lambda0 for mode 4
    d1_is_p = m1 | m2 | m3
    add(_safe_gather_int(p_trk, d1, d1_is_p), d1_is_p, True)
    add(_safe_gather_int(lam_p, d1, m4), m4, True)
    add(_safe_gather_int(lam_pi, d1, m4), m4, False)

    # d2: kaon (mode 1), K_S0 -> pi pi (modes 2/3), pion (mode 4)
    add(_safe_gather_int(data['KTrkIdx'], d2, m1), m1, False)
    add(_safe_gather_int(pi_trk, d2, m4), m4, False)

    # NOTE (limitation, modes 2 and 3): the K_S0's two pion tracks CANNOT be
    # resolved from the branches in these parquet files, so they are not
    # excluded here. K_Sd1Idx/K_Sd2Idx do not index any collection present in
    # the file: they exceed npi for 94% of K_S0 candidates, and in TRK space
    # the implied daughter charges are opposite only 55% of the time (random
    # is 50%). Where both indices do happen to fall inside the pi collection
    # (5.6% of candidates) the charges are opposite 100% of the time and d1
    # is always positive, matching K_Sd1Lund = +211 -- i.e. the intended index
    # space IS a pion list, but the one stored here is not it (the K_S0's
    # displaced daughters are evidently a separate/longer list upstream).
    # Kinematic matching against TRK was tried and is not clean either (the
    # daughters are refit at the displaced vertex, so K_SpreFitMass matches a
    # TRK pi+pi- pair to <1 MeV only 45% of the time -- consistent with
    # combinatorics, not identification).
    #
    # Consequence: for modes 2/3 the two K_S0 pions stay in the eligible pool
    # and can fake an antiproton, costing signal efficiency in those modes
    # only. This is measured per mode by run_stage04.py and reported at the
    # checkpoint rather than left implicit -- modes 1 and 4 have no K_S0 and
    # provide the clean comparison.

    # d3: pion for modes 1/3/4;  d4: pion for modes 3/4
    d3_is_pi = m1 | m3 | m4
    add(_safe_gather_int(pi_trk, d3, d3_is_pi), d3_is_pi, False)
    d4_is_pi = m3 | m4
    add(_safe_gather_int(pi_trk, d4, d4_is_pi), d4_is_pi, False)

    return slot_idx, slot_valid

################################################################################
def count_signal_b_tracks(data, config, scope='all_signal_tracks'):
    """
    Number of DISTINCT tracks used by each B candidate (jagged, aligned with
    the B collection). Validation handle for get_signal_b_track_slots: the
    expected RESOLVABLE final-state size is 4 for Lam0Lam0 and, for
    Lam0LamC, 5/3/5/7 for LambdaC modes 1/2/3/4 -- modes 2 and 3 fall 2
    short of their true 5/7 because the K_S0's pion tracks cannot be
    resolved from these files (see the note in get_signal_b_track_slots).
    Any deviation means the daughter walk is wrong, so run_stage04.py
    asserts on it rather than trusting the map.
    """
    slot_idx, slot_valid = get_signal_b_track_slots(data, config, scope)

    n = None
    for j in range(len(slot_idx)):
        # count slot j only if it is valid and its track was not already
        # counted in an earlier slot (see _first_occurrence rationale in
        # get_antibaryon_veto_mask)
        first = slot_valid[j]
        for k in range(j):
            first = first & ~(slot_valid[k] & (slot_idx[k] == slot_idx[j]))
        inc = ak.values_astype(first, np.int64)
        n = inc if n is None else (n + inc)

    return n

################################################################################
def get_antibaryon_veto_mask(data, config, selector=None,
                             scope='all_signal_tracks', pool='all_tracks'):
    """
    Antibaryon (antiproton) veto, as a candidate mask aligned with the B
    collection: True means the candidate SURVIVES the veto (no antibaryon
    found), so it ANDs directly with the Stage 2 purity / Stage 3 PID per-B
    masks.

    A B candidate is vetoed if the event contains a track that
      (a) is not used by THAT B candidate (see `scope`),
      (b) passes `selector` under the proton hypothesis, and
      (c) has electric charge opposite to that candidate's signal baryon
          charge, sign(Lambda0d1Lund) of its Lambda0 daughter -- so B and
          Bbar candidates flip automatically.

    selector=None returns an all-True (pass-all) mask, matching the Stage 3
    convention, so callers can scan the ladder one rung at a time or apply a
    config in which the veto is not yet enabled.

    pool: which tracks are eligible to fire the veto.
      'all_tracks' -- every track, using the per-track pSelectorsMap. This is
                      what the p-Lambda0 reference does, and it matters: 24%
                      of the tracks passing SuperLooseKMProtonSelection are
                      outside the 'p' hypothesis list (see the module header).
      'p_list'     -- only tracks in the 'p' hypothesis list (a weaker veto);
                      provided for the comparison reported at the checkpoint.

    Implementation note: rather than materializing a per-B x per-track
    boolean (O(nB x nTRK), prohibitive at ~60 B candidates/event), this
    counts the event's opposite-charge proton-selector tracks once and
    subtracts the ones belonging to each B candidate. The subtraction is
    exact only if each track is counted once, so slots sharing a track are
    de-duplicated explicitly below.
    """
    nB = ak.num(data['BpostFitMes'])

    if selector is None:
        return ak.ones_like(data['BpostFitMes'], dtype=bool)

    # ---- track-space quantities (per event) ----
    q_trk = np.sign(data['TRKLund'])
    pass_p = pid_selector.passes_selector(data['pSelectorsMap'], selector, 'p')

    if pool == 'p_list':
        in_p_list = indices_to_booleans(data['pTrkIdx'], data['pSelectorsMap'])
        pass_p = pass_p & in_p_list
    elif pool != 'all_tracks':
        raise ValueError(f"Unknown pool '{pool}'")

    # Tracks that could be an antiproton, split by charge, counted per event
    n_pos = ak.sum(pass_p & (q_trk > 0), axis=-1)
    n_neg = ak.sum(pass_p & (q_trk < 0), axis=-1)

    # ---- per-B signal baryon charge ----
    lam_idxvar = [iv for (n, iv) in config['B_daughters'] if n == 'Lambda0'][0]
    q_sig = np.sign(data['Lambda0d1Lund'])[data[lam_idxvar]]

    # Antibaryon charge is the opposite: count the event's candidates of that
    # charge, per B candidate.
    n_opp_event = ak.where(q_sig > 0, n_neg, n_pos)

    # ---- subtract the B candidate's own tracks that would have counted ----
    slot_idx, slot_valid = get_signal_b_track_slots(data, config, scope)

    n_opp_own = None
    for j in range(len(slot_idx)):
        idx_j, valid_j = slot_idx[j], slot_valid[j]

        # De-duplicate: only the first slot referencing a given track counts,
        # otherwise a track shared between two daughter slots is subtracted
        # twice and the candidate is under-vetoed.
        first = valid_j
        for k in range(j):
            first = first & ~(slot_valid[k] & (slot_idx[k] == idx_j))

        counts = (_safe_gather(pass_p, idx_j, valid_j) &
                  (_safe_gather_int(q_trk, idx_j, valid_j, fill=0) == -q_sig))
        inc = ak.values_astype(first & counts, np.int64)
        n_opp_own = inc if n_opp_own is None else (n_opp_own + inc)

    n_opp_other = n_opp_event - n_opp_own

    return n_opp_other == 0

################################################################################
def get_anti_lambda0_veto_mask(data, config):
    """
    Optional add-on evaluated (not adopted) at the Stage 4 checkpoint: veto a
    B candidate if the event contains a Lambda0-collection candidate of the
    OPPOSITE baryon sign to the signal, passing the channel's Lambda0 purity
    cuts, that is not itself a daughter of that B candidate.

    Expected to be largely redundant: the single-candidate requirement
    (nLambda0 == 2 for Lam0Lam0, == 1 for Lam0LamC) already rejects events
    carrying a spare Lambda0bar candidate, and any Lambda0bar contains a
    track the antiproton veto should already catch. Measured so that
    redundancy is a number in results/<channel>.yaml rather than an
    assumption.
    """
    lam_purity = get_composite_purity_mask(data, config, 'Lambda0')
    q_lam = np.sign(data['Lambda0d1Lund'])

    lam_idxvars = [iv for (n, iv) in config['B_daughters'] if n == 'Lambda0']
    q_sig = q_lam[data[lam_idxvars[0]]]

    # The signal B's own Lambda0 daughters always carry the signal sign
    # (verified, see the module header), so they never enter these counts and
    # need no explicit exemption.
    n_anti = ak.where(q_sig > 0,
                      ak.sum(lam_purity & (q_lam < 0), axis=-1),
                      ak.sum(lam_purity & (q_lam > 0), axis=-1))

    return n_anti == 0

################################################################################
# Diagnostic cuts for the early-stage notebooks
################################################################################
def build_diagnostic_cuts(data, config):
    """
    A minimal `dcuts` set for the load-and-look diagnostics (notebook 01):

        0: no cut
        1: single-candidate requirement
        2: 1 + B candidate in the mES/DeltaE fitting region

    The full cutflow (purity, PID, antibaryon veto, ...) is built up in the
    later analysis stages.
    """
    dcuts = {}

    nevents = len(data)

    dcuts[0] = {}
    dcuts[0]['name'] = 'no cut'
    dcuts[0]['event'] = ak.Array(np.ones(nevents, dtype=bool))
    dcuts[0]['candidates'] = None

    # Single-candidate requirement
    mask_single = get_single_candidate_mask(data, config)

    dcuts[1] = {}
    dcuts[1]['name'] = 'single candidate'
    dcuts[1]['event'] = mask_single
    dcuts[1]['candidates'] = None

    # Fitting region
    mask_candidates_fit = get_fit_mask(data, config)
    mask_event_fit = ak.any(mask_candidates_fit, axis=-1)

    dcuts[2] = {}
    dcuts[2]['name'] = 'single cand. + fit region'
    dcuts[2]['event'] = mask_single & mask_event_fit
    dcuts[2]['candidates'] = mask_candidates_fit

    return dcuts

################################################################################
# Cutflow table
################################################################################
def get_numbers_for_cut_flow(data, dcuts, spmodes=None, tag='DEFAULT', verbose=False):
    """
    Count events passing each cut in `dcuts`, per SP mode.

    Returns a DataFrame with columns cut/name/nevents/pct/tag/spmode,
    where pct is relative to the number of events for that SP mode
    before any cuts.
    """
    import pandas as pd

    if spmodes is None:
        spmodes = np.unique(np.array(data['spmode'].to_list()))
        if verbose:
            print(f"Using all spmodes in the file: {spmodes}")

    df_dict = {"cut": [], "name": [], "nevents": [], "pct": [], "tag": [], "spmode": []}

    for spmode in spmodes:
        mask_sp = data['spmode'] == spmode

        norg = ak.sum(mask_sp)

        # "org" row (skipped if dcuts already starts with a no-cut entry 0)
        if 0 not in dcuts:
            df_dict["cut"].append(0)
            df_dict["name"].append("org")
            df_dict["nevents"].append(norg)
            df_dict["pct"].append(100.0)
            df_dict["tag"].append(tag)
            df_dict["spmode"].append(spmode)

        for key in dcuts.keys():
            mask_event = dcuts[key]['event']
            n = ak.sum(mask_sp & mask_event)

            df_dict["cut"].append(key)
            df_dict["name"].append(dcuts[key]["name"])
            df_dict["nevents"].append(n)
            df_dict["pct"].append(100 * n / norg if norg > 0 else 0.0)
            df_dict["tag"].append(tag)
            df_dict["spmode"].append(spmode)

    return pd.DataFrame.from_dict(df_dict)
################################################################################
