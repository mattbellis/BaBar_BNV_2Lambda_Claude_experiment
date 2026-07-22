"""
Read/write the per-channel results files: results/<channel>.yaml

These files are the SINGLE SOURCE OF TRUTH for analysis numbers
(event counts, efficiencies, yields, systematics, limits, ...).
Each analysis stage writes its key results here, and
generate_latex_macros.py turns them into LaTeX macros/tables for the BAD.
Numbers are NEVER hand-typed into LaTeX (see CLAUDE.md).

Layout: one top-level section per analysis stage, e.g.

    channel: Lam0Lam0
    stage01:
      last_updated: '2026-07-22'
      n_events: {...}
      cutflow: {...}
"""

import datetime
import os

import yaml

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))

# <repo>/results
RESULTS_DIR = os.path.normpath(os.path.join(_THIS_DIR, '..', 'results'))

################################################################################
def results_filename(channel):
    return os.path.join(RESULTS_DIR, f"{channel}.yaml")

################################################################################
def load_results(channel):
    """Return the results dict for a channel ({} if no file yet)."""
    filename = results_filename(channel)
    if not os.path.exists(filename):
        return {}
    with open(filename) as f:
        return yaml.safe_load(f) or {}

################################################################################
def update_results(channel, section, values, verbose=True):
    """
    Merge `values` (a dict) into the `section` (e.g. 'stage01') of a
    channel's results file, stamping the date. Existing keys within the
    section are replaced wholesale; other sections are untouched.
    """
    os.makedirs(RESULTS_DIR, exist_ok=True)

    results = load_results(channel)
    results['channel'] = channel

    results[section] = dict(values)
    results[section]['last_updated'] = datetime.date.today().isoformat()

    filename = results_filename(channel)
    with open(filename, 'w') as f:
        yaml.safe_dump(results, f, sort_keys=False, default_flow_style=False)

    if verbose:
        print(f"Wrote section '{section}' to {filename}")

    return results
################################################################################
