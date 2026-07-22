#!/usr/bin/env python
"""
Stage 0/1: compute the load-and-diagnostics numbers and write them to
results/<channel>.yaml (the single source of truth for the BAD/paper).

This computes the same quantities shown in notebooks/01_load_and_diagnostics.ipynb:
event counts per SP mode, MC scaling weights, the diagnostic cutflow,
the blinding verification count, and (Lam0LamC) the LambdaC decay-mode
fractions.

Usage:
    python run_stage01.py --channel Lam0Lam0
    python run_stage01.py --channel Lam0LamC
    python run_stage01.py --channel all
"""

import argparse

import awkward as ak
import numpy as np

from channel_config import get_channel_config, BACKGROUND_SP_MODES, SIGNAL_SP_MODE
import cutflow
import datasets
import results_io

################################################################################
def region_candidate_counts(data, config):
    """Candidate counts (no cuts) in the signal/fit/sideband mES-DeltaE regions."""
    rd = config['region_definitions']
    mes = ak.flatten(data['BpostFitMes'])
    de = ak.flatten(data['BpostFitDeltaE'])

    def box(meslo, meshi, delo, dehi):
        return int(ak.sum((mes > meslo) & (mes < meshi) & (de > delo) & (de < dehi)))

    return {
        'signal': box(*rd['signal MES'], *rd['signal DeltaE']),
        'fit': box(*rd['fitting MES'], *rd['fitting DeltaE']),
        'sideband1': box(*rd['sideband MES'], *rd['sideband 1 DeltaE']),
        'sideband2': box(*rd['sideband MES'], *rd['sideband 2 DeltaE']),
    }

################################################################################
def run_stage01(channel):
    print(f"\n========== Stage 0/1 for {channel} ==========")
    config = get_channel_config(channel)

    data_sp, data_collision = datasets.load_datasets(channel)
    datasets.add_derived_fields(data_sp, config)
    datasets.add_derived_fields(data_collision, config)

    out = {}

    # ------------------------------------------------------------------
    # Event counts per SP mode
    # ------------------------------------------------------------------
    df_sp_counts = datasets.event_counts_by_spmode(data_sp)
    out['n_events_mc'] = {str(r['spmode']): int(r['nevents']) for _, r in df_sp_counts.iterrows()}
    out['n_events_data'] = int(len(data_collision))

    # ------------------------------------------------------------------
    # Integrated luminosity and MC scaling weights
    # ------------------------------------------------------------------
    dataset_information = datasets.read_in_dataset_statistics()
    mask = (dataset_information['Data or MC'] == 'Data') & \
           (dataset_information['Skim'] != 'LambdaVeryVeryLoose')
    out['int_lumi_invpb'] = float(dataset_information[mask]['Luminosity (Data only) 1/pb'].sum())

    bkg_present = [s for s in BACKGROUND_SP_MODES if s in out['n_events_mc']]
    weights = datasets.get_scaling_weights(bkg_present)
    out['scaling_weights'] = {k: float(v) for k, v in weights.items()}

    # ------------------------------------------------------------------
    # Diagnostic cutflow (raw event counts per SP mode per cut)
    # ------------------------------------------------------------------
    dcuts_sp = cutflow.build_diagnostic_cuts(data_sp, config)
    dcuts_col = cutflow.build_diagnostic_cuts(data_collision, config)

    out['cut_names'] = {int(k): v['name'] for k, v in dcuts_sp.items()}

    df_cf_sp = cutflow.get_numbers_for_cut_flow(data_sp, dcuts_sp, tag=channel)
    out['cutflow_mc'] = {}
    for spmode, grp in df_cf_sp.groupby('spmode'):
        out['cutflow_mc'][str(spmode)] = {int(r['cut']): int(r['nevents']) for _, r in grp.iterrows()}

    df_cf_col = cutflow.get_numbers_for_cut_flow(data_collision, dcuts_col, tag=channel)
    out['cutflow_data'] = {int(r['cut']): int(r['nevents']) for _, r in df_cf_col.iterrows()}

    # ------------------------------------------------------------------
    # Region definitions (recorded so the BAD quotes exactly what ran)
    # and candidate counts per region
    # ------------------------------------------------------------------
    rd = config['region_definitions']
    out['region_definitions'] = {k: v for k, v in rd.items() if k != 'inference'}

    out['region_counts_data'] = region_candidate_counts(data_collision, config)

    mask_sig = data_sp['spmode'] == SIGNAL_SP_MODE
    out['region_counts_signal_mc'] = region_candidate_counts(data_sp[mask_sig], config)

    # ------------------------------------------------------------------
    # Blinding verification: data candidates in the signal region MUST be 0
    # ------------------------------------------------------------------
    nsig_data = out['region_counts_data']['signal']
    out['blinding_check'] = {'n_data_candidates_in_signal_region': nsig_data,
                             'passed': bool(nsig_data == 0)}
    if nsig_data != 0:
        print(f"*** WARNING: BLINDING CHECK FAILED: {nsig_data} data candidates "
              f"in the assumed signal region ***")

    # ------------------------------------------------------------------
    # Composite-candidate configuration snapshot (mass windows etc.)
    # ------------------------------------------------------------------
    out['composites'] = {
        name: {'n_required': comp['n_required'],
               'mass_window': [float(comp['mass_window'][0]), float(comp['mass_window'][1])],
               'flight_var': comp['flight_var'],
               'flight_cut': comp['flight_cut']}
        for name, comp in config['composites'].items()
    }

    # ------------------------------------------------------------------
    # LambdaC decay-mode fractions (Lam0LamC only), signal MC
    # ------------------------------------------------------------------
    if 'LambdaC' in config['composites']:
        lamc_mode = cutflow.get_lambdac_decay_mode(data_sp)

        def fractions(mode_flat):
            vals, counts = np.unique(ak.to_numpy(mode_flat), return_counts=True)
            return {int(v): float(c / counts.sum()) for v, c in zip(vals, counts)}

        out['lambdac_mode_fractions'] = {
            'no_cut': fractions(ak.flatten(lamc_mode[mask_sig])),
            'single_candidate': fractions(
                ak.flatten(lamc_mode[dcuts_sp[1]['event'] & mask_sig])),
        }

    results_io.update_results(channel, 'stage01', out)

    return out

################################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--channel', required=True,
                        choices=['Lam0Lam0', 'Lam0LamC', 'all'])
    args = parser.parse_args()

    channels = ['Lam0Lam0', 'Lam0LamC'] if args.channel == 'all' else [args.channel]
    for channel in channels:
        run_stage01(channel)
