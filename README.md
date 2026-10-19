# FAERS and JADER antithyroid pharmacovigilance analysis

This repository contains the analysis code for the study of hepatic and
hematologic disproportional-reporting signals for methimazole and
propylthiouracil in FAERS and JADER.

Repository: <https://github.com/Hajimi-Sudo/faers-jader-antithyroid-safety-signals>

The code is designed for execution on an analysis server. Raw FAERS/JADER
files, report-level derived data, publication full text, credentials, and
manuscript build artifacts are intentionally excluded from this release.
The analysis is a pharmacovigilance signal study. Its estimates are not
incidence, relative-risk, causal-effect, treatment-effect, or confirmed
drug-drug-interaction estimates.

## Repository layout

```text
src/
  analysis/       FAERS, JADER, sensitivity, and aggregate-table scripts
  figures/        deterministic figure-generation scripts
  case_evidence/ published-case evidence search and review scripts
  qc/             submission and evidence-quality checks
  manuscript/     remote aggregate manuscript-generation script
workflows/        server-side execution order and output contract
docs/             method and data-boundary notes
```

## Runtime boundary

Set `FAERS_PROJECT_ROOT` to the project directory on the analysis server. The
expected default is `/root/autodl-tmp/faers-antithyroid-safety-signals`.
Do not place raw or report-level FAERS/JADER data in a public repository.
The scripts read and write only beneath `FAERS_PROJECT_ROOT` when that
variable is set.

## Dependencies

- Python 3.10 or newer
- Python standard library
- `matplotlib` and `numpy` for figure generation
- Perl 5 for the optional Perl implementations and verification checks
- Network access from the server only for the PubMed/Europe PMC evidence
  acquisition scripts

The exact server environment and frozen data-window identifiers should be
recorded in the project audit directory rather than committed with the code.

## Reproducibility notes

The primary FAERS window and outcome dictionaries are frozen in the scripts
and manuscript methods. Serious outcomes use `DE`, `LT`, `HO`, `DS`, `CA`, and
`OT`. `RI` is excluded because it denotes intervention to prevent permanent
impairment. The report unit is the latest case version per case ID, with
database-specific rules documented in the manuscript.

See [`workflows/REMOTE_RUN_ORDER.md`](workflows/REMOTE_RUN_ORDER.md) for the
server-only execution order and expected aggregate outputs.

## License

Add the repository license selected by the authors before public release. The
source code is provided for reproducibility of the accompanying study and
does not grant access to restricted or redistributed adverse-event records.
