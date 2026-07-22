#!/bin/bash
# Copy the diagnostic plots produced by the notebooks into the BAD figures
# directory. Run from BNV_2Lambda_analysis/ after (re)running the notebooks:
#
#     ./copy_plots_to_BAD.sh
#
# The copies in BAD_2Lambda/figures/ ARE committed to git (the notebook
# output directory notebooks/plots/ is not), so the BAD always builds.

set -e
cd "$(dirname "$0")"

SRC=notebooks/plots
DEST=../BAD_2Lambda/figures

mkdir -p "$DEST"
cp -rv "$SRC"/* "$DEST"/
