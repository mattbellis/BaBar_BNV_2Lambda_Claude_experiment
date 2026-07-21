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
 
