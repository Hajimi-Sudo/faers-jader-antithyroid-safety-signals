#!/usr/bin/env python3
"""Print a compact remote-only digest for manual case-evidence review."""
import csv
import os
from pathlib import Path

root = Path(os.environ.get("FAERS_PROJECT_ROOT", "/root/autodl-tmp/faers-antithyroid-safety-signals"))
path = root / "literature_case_evidence" / "manual_case_review_sheet.tsv"
fields = [
    "age_sex_snippets", "dose_snippets", "latency_snippets",
    "dechallenge_snippets", "rechallenge_snippets",
    "alternative_cause_snippets", "outcome_snippets",
]
with path.open(encoding="utf-8", newline="") as handle:
    rows = csv.DictReader(handle, delimiter="\t")
    for i, row in enumerate(rows, 1):
        print(f"\n[{i}] PMID {row['pmid']} | {row['year']} | {row['drug_class']} | {row['outcome_class']}")
        print(row["title"])
        for field in fields:
            value = " ".join(row[field].split())
            if value:
                print(f"{field}: {value[:750]}")
