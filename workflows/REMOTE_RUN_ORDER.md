# Remote run order

Run these commands on the assigned analysis server after setting the project
root. Do not run them in the local manuscript workspace.

```bash
export FAERS_PROJECT_ROOT=/root/autodl-tmp/faers-antithyroid-safety-signals
export PYTHONPATH="$FAERS_PROJECT_ROOT/code:$PYTHONPATH"
```

The source files in this release may be copied into the server project's code
directory. The raw data and all report-level intermediate files remain on the
server.

## Primary analysis

```bash
python src/analysis/run_faers_full_window.py "$FAERS_PROJECT_ROOT"
python src/analysis/run_gate1_sensitivity.py "$FAERS_PROJECT_ROOT"
python src/analysis/remote_faers_supplement.py "$FAERS_PROJECT_ROOT" build
python src/analysis/remote_faers_supplement.py "$FAERS_PROJECT_ROOT" ror
python src/analysis/remote_faers_supplement.py "$FAERS_PROJECT_ROOT" ic
python src/analysis/remote_faers_supplement.py "$FAERS_PROJECT_ROOT" temporal
python src/analysis/remote_faers_supplement.py "$FAERS_PROJECT_ROOT" interaction
```

## JADER and cross-database summary

```bash
python src/analysis/run_jader_pilot.py "$FAERS_PROJECT_ROOT"
python src/analysis/make_jader_concordance.py
```

The JADER script uses `FAERS_PROJECT_ROOT` for its input and output paths.
Keep JADER estimates separate from FAERS estimates and do not pool them.

## Figures

```bash
python src/figures/make_figure1_nature.py
python src/figures/make_figure2_calendar_nature.py
python src/figures/make_figure3_cross_database_nature.py
python src/figures/make_figure4_ic_nature.py
```

## Published-case evidence layer

The case-evidence scripts require the server-side literature evidence
directory and network access. Execute them in the order recorded by the
existing project protocol, then run the verification and QC scripts. The
resulting publication-level evidence is contextual and is not an independent
patient count or a systematic review.

## Quality control

```bash
python src/qc/remote_submission_qc.py
perl src/qc/remote_submission_qc.pl
```

Before public release, inspect the generated audit JSON and TSV files on the
server. Do not commit them when they contain report-level or restricted
publication data.
