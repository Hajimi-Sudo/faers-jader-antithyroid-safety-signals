#!/usr/bin/env python3
"""Build a remote article-level case review sheet from screened full texts."""
from __future__ import annotations

import csv
import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(os.environ.get("FAERS_PROJECT_ROOT", "/root/autodl-tmp/faers-antithyroid-safety-signals"))
BASE = ROOT / "literature_case_evidence"

PATTERNS = {
    "drug": re.compile(r"\b(?:methimazole|thiamazole|propylthiouracil|PTU)\b", re.I),
    "age_sex": re.compile(r"\b(?:\d{1,3}[- ]year[- ]old|aged?\s+\d{1,3}|\d{1,3}\s*(?:yo|y/o)|male|female|woman|man|girl|boy)\b", re.I),
    "dose": re.compile(r"\b(?:\d+(?:\.\d+)?\s*(?:mg|g|mcg|μg)(?:\s*(?:daily|per day|/day|once|twice|three times))?|dose|dosage)\b", re.I),
    "latency": re.compile(r"\b(?:after\s+\d+\s+(?:day|week|month|year)s?|\d+\s+(?:day|week|month|year)s?\s+(?:after|of|following)|within\s+\d+\s+(?:day|week|month|year)s?)\b", re.I),
    "dechallenge": re.compile(r"\b(?:discontinu|withdraw|stop|cessat|withdrawal|after stopping|after discontinuing|upon discontinuation|improv|recover|resolv|normaliz)\w*\b", re.I),
    "rechallenge": re.compile(r"\b(?:rechallenge|re-exposure|reexposure|reintroduced|readminister|second exposure)\b", re.I),
    "alternative": re.compile(r"\b(?:thyroid storm|thyrotoxicosis|viral hepatitis|autoimmune hepatitis|hepatitis B|sepsis|ischemia|shock liver|Graves' disease)\b", re.I),
    "outcome": re.compile(r"\b(?:death|died|fatal|transplant|liver transplantation|hospitali[sz]|ICU|intensive care|sepsis|discharged|surviv|resolved|recovered|mortality)\w*\b", re.I),
}


def text_from_xml(path: Path) -> str:
    try:
        root = ET.parse(path).getroot()
        return " ".join(" ".join(root.itertext()).split())
    except Exception:
        return ""


def snippets(blob: str, pattern: re.Pattern[str], limit: int = 4) -> str:
    out = []
    for match in list(pattern.finditer(blob))[:limit]:
        out.append(blob[max(0, match.start() - 220):min(len(blob), match.end() + 420)])
    return " || ".join(out)


def main() -> None:
    rows = []
    with (BASE / "screening_audit.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("evidence_tier") != "direct_case_signal":
                continue
            path = BASE / "fulltext_xml" / f"{row.get('pmcid', '')}.xml"
            blob = text_from_xml(path)
            rows.append({
                "pmid": row.get("pmid", ""),
                "pmcid": row.get("pmcid", ""),
                "doi": row.get("doi", ""),
                "year": row.get("year", ""),
                "title": row.get("title", ""),
                "drug_class": row.get("drug_class", ""),
                "outcome_class": row.get("outcome_class", ""),
                "age_sex_snippets": snippets(blob, PATTERNS["age_sex"]),
                "dose_snippets": snippets(blob, PATTERNS["dose"]),
                "latency_snippets": snippets(blob, PATTERNS["latency"]),
                "dechallenge_snippets": snippets(blob, PATTERNS["dechallenge"]),
                "rechallenge_snippets": snippets(blob, PATTERNS["rechallenge"]),
                "alternative_cause_snippets": snippets(blob, PATTERNS["alternative"]),
                "outcome_snippets": snippets(blob, PATTERNS["outcome"]),
                "manual_verification_status": "manual_verification_pending",
                "inclusion_decision": "pending",
                "reviewer_notes": "",
            })
    fields = list(rows[0]) if rows else []
    with (BASE / "manual_case_review_sheet.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} review candidates")


if __name__ == "__main__":
    main()
