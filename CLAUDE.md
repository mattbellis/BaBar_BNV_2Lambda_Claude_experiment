# CLAUDE.md - BaBar BNV Search: B -> Lambda0 Lambda0 and B -> LambdaC Lambda0
 
## Project overview
 
This repository contains two searches for baryon-number-violating (Delta-B = 2)
decays of B mesons, using BaBar data (424.3 fb^-1 at the Upsilon(4S)):
 
- **B -> Lambda0 Lambda0**
- **B -> LambdaC Lambda0**
These follow the completed **B+ -> p Lambda0** analysis (upper limit
6.1e-8 at 90% CL), whose code lives in this repo **for reference only** -- do not
modify it. The new analyses reuse its strategy: blind search in the mES / Delta-E
plane, Lambda0 purity cuts (flight-length significance + mass window), PID
selectors optimized with the Punzi figure of merit, an antibaryon veto, an MLP
classifier on event-shape variables, and a Rolke-style upper-limit calculation
with sideband-based background estimation.
 
Deliverables: analysis code, a BAD (BaBar Analysis Document, internal LaTeX
note), and eventually a paper draft.
 
## Data
 
- ROOT files have already been converted to **parquet + awkward array** format
  (signal MC, background MC, and collision data). Work from the parquet files.
- Do not attempt to reprocess or re-download data.
## Collaboration rules (ground rules -- read carefully)
 
1. **The human runs the real analysis.** Claude writes and edits code; the human
   executes full analysis runs. Claude may run quick tests on small samples with
   permission, but the pipeline must always be runnable by the human alone,
   end-to-end, without Claude.
2. **Staged development with checkpoints.** The analysis proceeds in stages (see
   below). **Stop at each stage boundary** and wait for the human to review
   results, plots, and code before proceeding. Do not run ahead.
3. **Blinding is sacrosanct.** The signal region in collision data (approx.
   5.27 < mES < 5.30 GeV/c^2 and |Delta-E| < 0.05 GeV; exact windows to be set
   per channel) stays masked until the human explicitly approves unblinding.
   Never write or run code that looks at collision data in the signal region
   before that approval. If a request seems to require it, flag it and stop.
4. **Plan before implementing.** For any nontrivial design decision (background
   treatment, reconstruction choices, new cuts, classifier changes), discuss and
   agree on a plan before writing code.
5. **Keep WORKFLOW-README.md current.** Any change to how the pipeline is run --
   new script, changed arguments, new dependency, reordered steps -- must be
   reflected in WORKFLOW-README.md in the same commit. That document must always
   let the human run everything independently.
## Analysis stages
 
Each stage ends with a validation checkpoint (plots/numbers the human reviews)
and a git commit.
 
1. Preselection / skims -- loose cuts, single-candidate requirements
2. Lambda0 (and LambdaC) purity -- flight-length significance, mass windows
3. PID selector optimization (Punzi FOM)
4. Antibaryon veto
5. MLP classifier -- training, feature validation, output-cut optimization
   (Punzi FOM with Delta-E-sideband background estimate)
6. Efficiency, sideband studies, systematics
7. Unblinding (human approval required) and Rolke upper limit
## Code conventions
 
- **Parametrize by channel from day one.** One codebase; a channel configuration
  selects Lambda0-Lambda0 vs. LambdaC-Lambda0 (and reconstruction mode). No
  copied-and-diverging per-channel directories.
- Python + awkward/parquet, matching the style of the p-Lambda0 reference code.
- **Single source of truth for numbers:** analysis stages write key results
  (efficiencies, yields, tau, systematics, limits) to results/<channel>.yaml.
  A small script generates LaTeX macros from these for the BAD/paper. Never
  hand-type analysis numbers into LaTeX.
- Prefer scripts/notebooks invocable per stage (e.g., a Makefile or
  run_stage_N.py --channel ...) so the human can drive the pipeline.
## BAD (internal note)
 
- LaTeX, in this repo; build with latexmk; fix compile errors from the log.
- **Written incrementally**: draft the corresponding BAD section as part of each
  stage checkpoint, while decisions and plots are fresh. Use the p-Lambda0 BAD
  as the structural template.
- Prose must describe what the code actually does -- read the code, don't guess.
- All numbers come from the generated LaTeX macros (see above).
## Git
 
- Commit at each validated stage checkpoint, with clear messages describing the
  physics/analysis content of the change.
- Ask before committing; **never push without explicit approval.**
- Never commit large data files (parquet, ROOT); respect .gitignore.
