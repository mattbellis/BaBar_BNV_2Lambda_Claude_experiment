"""
Dataset loading and MC-scaling utilities for the BNV 2-Lambda analyses.

The parquet files live in <repo>/data/ and follow the naming convention

    Background_and_signal_SP_modes_All_runs_<channel>.parquet   (MC)
    Data_All_runs_<channel>_BLINDED.parquet                     (collision data)

where <channel> is 'Lam0Lam0' or 'Lam0LamC'.

BLINDING: load_datasets() only ever opens the _BLINDED collision files.
Unblinding requires explicit human approval (see CLAUDE.md); the UNBLINDED
flag raises until that approval is recorded here.
"""

import os
import time

import awkward as ak
import numpy as np
import pandas as pd

# Directory of this file, so the defaults below work regardless of
# the notebook/script working directory
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))

# Default location of the parquet files: <repo>/data
DEFAULT_DATA_DIR = os.path.normpath(os.path.join(_THIS_DIR, '..', 'data'))

################################################################################
def load_datasets(channel, topdir=None, sp_file_tag='Background_and_signal_SP_modes',
                  collision_file_tag='Data', subset='all', sp_or_data=None,
                  UNBLINDED=False):
    """
    Load the MC ("SP") and/or collision-data awkward arrays for a channel.

    channel:      'Lam0Lam0' or 'Lam0LamC'
    topdir:       directory holding the parquet files (default: <repo>/data)
    subset:       'all' -> 'All_runs' files ('run1' -> 'Only_Run_1', if produced)
    sp_or_data:   'sp' for MC only, 'col' for collision only, None for both
    UNBLINDED:    must remain False; raises otherwise (blinding policy)

    Returns (data_sp, data_collision); entries are None if not requested.
    """

    if UNBLINDED is True:
        raise RuntimeError(
            "Unblinding requires explicit human approval (see CLAUDE.md). "
            "This function only loads the _BLINDED collision files."
        )

    if topdir is None:
        topdir = DEFAULT_DATA_DIR

    subset_tag = ""
    if subset.lower() == 'run1':
        subset_tag = 'Only_Run_1'
    elif subset.lower() == 'all':
        subset_tag = 'All_runs'
    else:
        raise ValueError(f"Unknown subset '{subset}' (use 'all' or 'run1')")

    data_sp, data_collision = None, None

    if sp_or_data in ('sp', 'SP', None):
        start = time.time()
        filename = f"{topdir}/{sp_file_tag}_{subset_tag}_{channel}.parquet"
        print(f"Opening {filename}...")
        data_sp = ak.from_parquet(filename)
        print(f"Took {time.time()-start:.3f} seconds\n")

    if sp_or_data in ('col', 'collision', None):
        start = time.time()
        filename = f"{topdir}/{collision_file_tag}_{subset_tag}_{channel}_BLINDED.parquet"
        print(f"Opening {filename}...")
        data_collision = ak.from_parquet(filename)
        print(f"Took {time.time()-start:.3f} seconds\n")

    return data_sp, data_collision
################################################################################

################################################################################
def read_in_dataset_statistics(infilename=None):
    """Numbers of events/luminosity per run and skim (copied from p-Lambda0)."""
    if infilename is None:
        infilename = os.path.join(_THIS_DIR, 'dataset_statistics.csv')
    return pd.read_csv(infilename)
################################################################################

################################################################################
def get_SP_cross_sections_and_labels(infilename=None):
    """Cross sections and labels for the SP background modes."""
    if infilename is None:
        infilename = os.path.join(_THIS_DIR, 'SP_cross_sections_and_labels.csv')
    return pd.read_csv(infilename)
################################################################################

################################################################################
def scaling_value(spmode, dataset_information=None, cs_data=None, verbose=False):
    """
    Weight to scale one generated MC event of a given SP mode to the
    luminosity of the full dataset:

        w = (cross section * integrated luminosity) / (# generated events)

    spmode can be a string or int. Data ('0' / -1) and signal MC ('-999')
    return 1 -- signal normalization is handled separately.
    """
    mode = int(spmode)

    if mode in (0, -1, -999):
        return 1.0

    if dataset_information is None:
        dataset_information = read_in_dataset_statistics()
    if cs_data is None:
        cs_data = get_SP_cross_sections_and_labels()

    # Number of generated (unskimmed) MC events for this SP mode
    mc_mask = (dataset_information["SP mode"] == mode) & \
              (dataset_information["Skim"] != "LambdaVeryVeryLoose")
    nevents_mc = dataset_information[mc_mask]["# of events (Data or MC)"].sum()

    if nevents_mc == 0:
        print(f"WARNING: no generated-event counts found for SP mode {mode}; returning weight 1")
        return 1.0

    # Cross section [nb]
    cs_mask = cs_data["SP Mode"] == mode
    if cs_mask.sum() == 0:
        print(f"WARNING: no cross section found for SP mode {mode}; returning weight 1")
        return 1.0
    cs = cs_data[cs_mask]["Cross section [nb]"].values[0]

    # Integrated luminosity of the data [1/pb]
    mask = (dataset_information['Data or MC'] == 'Data') & \
           (dataset_information['Skim'] != 'LambdaVeryVeryLoose')
    int_lumi = dataset_information[mask]['Luminosity (Data only) 1/pb'].sum()

    # Factor of 1000 converts nb * 1/pb
    n_exp_in_data = cs * int_lumi * 1000

    scaling = n_exp_in_data / nevents_mc

    if verbose:
        print(f"- Cross section for SP mode {mode}:      {cs} nb")
        print(f"- # of events generated for SP-{mode}: {nevents_mc:13d}")
        print(f"- Number expected in data:             {n_exp_in_data:.1f}")
        print(f"- Integrated luminosity:               {int_lumi:.1f} 1/pb")
        print(f"- Scaling value:                       {scaling:.4f}")

    return scaling
################################################################################

################################################################################
def get_scaling_weights(spmodes, verbose=False):
    """
    Convenience: return a dict {spmode(str): weight} for a list of SP modes.
    """
    dataset_information = read_in_dataset_statistics()
    cs_data = get_SP_cross_sections_and_labels()

    weights = {}
    for spmode in spmodes:
        weights[str(spmode)] = scaling_value(spmode,
                                             dataset_information=dataset_information,
                                             cs_data=cs_data,
                                             verbose=verbose)
    return weights
################################################################################

################################################################################
def add_derived_fields(data, config):
    """
    Add channel-dependent derived fields to an awkward array (in place):

    All channels:
      - Lambda0FLSig_fromB: flight significance of the Lambda0(s) that are
        direct B daughters (one entry per B-daughter Lambda0)

    Lam0LamC only:
      - LambdaCFlightSignificance: LambdaCFlightLen / LambdaCFlightErr
        (no post-fit flight-significance branch exists for the LambdaC)
      - Lambda0FLSig_fromLambdaC: flight significance of the Lambda0 used
        inside a LambdaC candidate (mode 4, LambdaC -> Lam0 pi pi pi;
        empty list for events with no such LambdaC)

    These can then be histogrammed with the standard machinery
    (see the corresponding entries in the channel's hist_defs).
    """
    flsig = data['Lambda0postFitFlightSignificance']

    # Lambda0(s) that are direct B daughters
    idxvars = [iv for (name, iv) in config['B_daughters'] if name == 'Lambda0']
    parts = [flsig[data[iv]] for iv in idxvars]
    data['Lambda0FLSig_fromB'] = parts[0] if len(parts) == 1 else ak.concatenate(parts, axis=1)

    if 'LambdaC' in config['composites']:
        data['LambdaCFlightSignificance'] = data['LambdaCFlightLen'] / data['LambdaCFlightErr']

        # Lambda0 daughters of LambdaC candidates (LambdaCd1Idx indexes the
        # Lambda0 collection only when d1 *is* a Lambda0)
        sel = abs(data['LambdaCd1Lund']) == 3122
        data['Lambda0FLSig_fromLambdaC'] = flsig[data['LambdaCd1Idx'][sel]]

    print(f"Added derived fields for channel {config['name']}")

################################################################################
def event_counts_by_spmode(data):
    """
    Return a DataFrame with the raw number of events for each SP mode
    in an awkward array (uses the event-level 'spmode' field).
    """
    spmodes = np.array(data['spmode'].to_list())
    vals, counts = np.unique(spmodes, return_counts=True)

    df = pd.DataFrame({'spmode': vals, 'nevents': counts})
    return df.sort_values('nevents', ascending=False).reset_index(drop=True)
################################################################################
