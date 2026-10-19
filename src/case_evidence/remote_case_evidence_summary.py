#!/usr/bin/env python3
"""Create a conservative article-level case-evidence summary on the server.

This produces a review aid, not an automated causality assessment.
"""
from __future__ import annotations

import csv
import os
import re
from pathlib import Path

ROOT = Path(os.environ.get("FAERS_PROJECT_ROOT", "/root/autodl-tmp/faers-antithyroid-safety-signals"))
BASE = ROOT / "literature_case_evidence"


def yes_no(pattern: str, text: str) -> str:
    return "reported" if re.search(pattern, text, re.I) else "not_identified_in_screening"


def first(pattern: str, text: str, n: int = 900) -> str:
    m = re.search(pattern, text, re.I)
    return " ".join(text[max(0, m.start() - 160):m.end() + n].split()) if m else ""


def main() -> None:
    rows = []
    with (BASE / "manual_case_review_sheet.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            blob = " ".join(row.get(k, "") for k in (
                "age_sex_snippets", "dose_snippets", "latency_snippets",
                "dechallenge_snippets", "rechallenge_snippets",
                "alternative_cause_snippets", "outcome_snippets"
            ))
            rows.append({
                "pmid": row["pmid"],
                "pmcid": row["pmcid"],
                "doi": row["doi"],
                "year": row["year"],
                "title": row["title"],
                "drug_class": row["drug_class"],
                "outcome_class": row["outcome_class"],
                "age_sex_status": "reported" if row["age_sex_snippets"] else "unclear",
                "dose_status": "reported" if row["dose_snippets"] else "unclear",
                "latency_status": "reported" if row["latency_snippets"] else "unclear",
                "dechallenge_status": yes_no(r"discontinu|withdraw|stop|improv|recover|resolv|normaliz", row["dechallenge_snippets"]),
                "rechallenge_status": yes_no(r"rechallenge|re-exposure|reexposure|reintroduced|readminister|second exposure", row["rechallenge_snippets"]),
                "alternative_cause_status": yes_no(r"thyroid storm|thyrotoxicosis|viral hepatitis|autoimmune hepatitis|hepatitis B|sepsis|ischemia|shock liver|Graves", row["alternative_cause_snippets"]),
                "outcome_status": "reported" if row["outcome_snippets"] else "unclear",
                "case_count_status": "unclear",
                "causality_status": "manual_adjudication_required",
                "include_in_final_case_summary": "pending",
                "evidence_note": "Screening-derived fields only; verify against full text before manuscript use.",
            })
    out = BASE / "case_evidence_summary.tsv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} conservative case-evidence rows")


if __name__ == "__main__":
    main()
