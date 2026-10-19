#!/usr/bin/env python3
"""Screen remote PubMed candidates and extract article-level case-evidence fields.

The output is an auditable screening aid. Heuristic matches are not treated as
clinical adjudication; uncertain fields remain explicitly marked.
"""
from __future__ import annotations

import csv
import json
import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(os.environ.get("FAERS_PROJECT_ROOT", "/root/autodl-tmp/faers-antithyroid-safety-signals"))
BASE = ROOT / "literature_case_evidence"
XML = BASE / "fulltext_xml"

DRUG_PATTERNS = {
    "methimazole": re.compile(r"\b(?:methimazole|thiamazole)\b", re.I),
    "propylthiouracil": re.compile(r"\b(?:propylthiouracil|ptu)\b", re.I),
}
OUTCOME_PATTERNS = {
    "hepatic": re.compile(r"\b(?:hepatotoxicity|hepatitis|hepatic failure|liver injury|liver failure|jaundice|cholestasis)\b", re.I),
    "hematologic": re.compile(r"\b(?:agranulocytosis|neutropenia|pancytopenia|bone marrow failure|myelosuppression)\b", re.I),
}
DIRECT_PATTERNS = {
    "drug_attribution": re.compile(
        r"(?:methimazole|thiamazole|propylthiouracil|\bPTU\b)\s*[- ]?(?:induced|associated|related|caused|toxicity|hepatotoxicity|agranulocytosis|neutropenia|liver injury|liver failure)",
        re.I,
    ),
    "dechallenge": re.compile(
        r"(?:discontinu|withdraw|stop|cessat|withdrawal|after stopping|after discontinuing|upon discontinuation|dechallenge).{0,220}(?:improv|recover|normal|resolv|decreas)",
        re.I,
    ),
    "rechallenge": re.compile(r"rechallenge|re-exposure|reexposure|reintroduced|reintroduced", re.I),
    "alternative_cause": re.compile(r"thyroid storm|thyrotoxicosis|viral hepatitis|autoimmune hepatitis|sepsis|ischemia|shock liver", re.I),
}
EXCLUSION_PATTERNS = {
    "review": re.compile(r"\breview\b|systematic review|meta-analysis", re.I),
    "database_only": re.compile(r"FAERS|JADER|disproportionality|spontaneous reporting", re.I),
    "nonhuman": re.compile(r"\b(?:rat|mice|mouse|animal model|in vitro|cell line)\b", re.I),
}


def clean(value: str) -> str:
    return " ".join((value or "").split())


def xml_text(path: Path) -> str:
    try:
        root = ET.parse(path).getroot()
        return clean(" ".join(root.itertext()))
    except Exception:
        return ""


def extract_case_sections(path: Path) -> str:
    """Return a bounded text excerpt around case-report sections for review."""
    try:
        root = ET.parse(path).getroot()
    except Exception:
        return ""
    chunks = []
    for sec in root.findall(".//sec"):
        title = clean(" ".join(sec.findtext("title", default="").split()))
        if re.search(r"case|patient|clinical|outcome|treatment|laboratory", title, re.I):
            body = clean(" ".join(sec.itertext()))
            if body:
                chunks.append(f"[{title}] {body}")
    return " ".join(chunks)[:12000]


def first_match(pattern: re.Pattern[str], text: str) -> str:
    match = pattern.search(text)
    return match.group(0) if match else ""


def main() -> None:
    source = BASE / "pubmed_candidates.tsv"
    screening_rows = []
    with source.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            title_abs = clean(" ".join([row.get("title", ""), row.get("abstract", "")]))
            full_path = Path(row.get("fulltext_xml", "")) if row.get("fulltext_xml") else XML / f"{row.get('pmcid', '')}.xml"
            has_fulltext = full_path.exists() and full_path.stat().st_size > 1000
            blob = title_abs + " " + (xml_text(full_path)[:300000] if has_fulltext else "")
            drugs = [name for name, pattern in DRUG_PATTERNS.items() if pattern.search(blob)]
            outcomes = [name for name, pattern in OUTCOME_PATTERNS.items() if pattern.search(blob)]
            pubtypes = row.get("pubtypes", "")
            reasons = []
            if not drugs:
                reasons.append("no_target_drug")
            if not outcomes:
                reasons.append("no_prespecified_outcome")
            if not has_fulltext:
                reasons.append("no_accessible_fulltext")
            if EXCLUSION_PATTERNS["review"].search(row.get("title", "") + " " + row.get("pubtypes", "")):
                reasons.append("review_or_synthesis")
            if EXCLUSION_PATTERNS["database_only"].search(row.get("title", "") + " " + row.get("abstract", "")):
                reasons.append("database_or_disproportionality_study")
            if EXCLUSION_PATTERNS["nonhuman"].search(row.get("title", "") + " " + row.get("abstract", "")):
                reasons.append("nonhuman_or_in_vitro")
            likely_case = "Case Reports" in pubtypes or bool(re.search(r"case report|case series|patient", row.get("title", "") + " " + row.get("abstract", ""), re.I))
            if not likely_case:
                reasons.append("not_case_report_like")
            include = not reasons and likely_case
            case_excerpt = extract_case_sections(full_path) if include else ""
            excerpt_blob = title_abs + " " + case_excerpt
            direct_attribution = bool(DIRECT_PATTERNS["drug_attribution"].search(row.get("title", "") + " " + excerpt_blob))
            dechallenge = bool(DIRECT_PATTERNS["dechallenge"].search(excerpt_blob))
            rechallenge = bool(DIRECT_PATTERNS["rechallenge"].search(excerpt_blob))
            alternative_cause = bool(DIRECT_PATTERNS["alternative_cause"].search(excerpt_blob))
            carbimazole_only = bool(re.search(r"\bcarbimazole\b", excerpt_blob, re.I)) and not bool(re.search(r"\bmethimazole\b", excerpt_blob, re.I))
            if carbimazole_only:
                reasons.append("carbimazole_only_or_not_methimazole")
                include = False
            evidence_tier = "direct_case_signal" if include and direct_attribution else "case_context_only" if include else "excluded"
            n_cases = "reported_or_unclear"
            if include:
                counts = re.findall(r"\b(?:one|two|three|four|five|six|seven|eight|nine|ten|\d+)\s+(?:patients?|cases?)\b", excerpt_blob, re.I)
                if counts:
                    n_cases = "; ".join(counts[:5])
            screening_rows.append({
                "pmid": row.get("pmid", ""),
                "pmcid": row.get("pmcid", ""),
                "doi": row.get("doi", ""),
                "year": row.get("year", ""),
                "journal": row.get("journal", ""),
                "title": row.get("title", ""),
                "drug_class": ";".join(drugs),
                "outcome_class": ";".join(outcomes),
                "has_accessible_fulltext": str(has_fulltext).lower(),
                "likely_case_report": str(likely_case).lower(),
                "include_for_manual_case_review": str(include).lower(),
                "exclusion_reasons": ";".join(reasons),
                "case_count_text": n_cases,
                "direct_drug_attribution_signal": str(direct_attribution).lower(),
                "dechallenge_signal": str(dechallenge).lower(),
                "rechallenge_signal": str(rechallenge).lower(),
                "alternative_cause_signal": str(alternative_cause).lower(),
                "evidence_tier": evidence_tier,
                "case_excerpt": case_excerpt,
            })

    fields = list(screening_rows[0]) if screening_rows else []
    out = BASE / "screening_audit.tsv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(screening_rows)
    included = [x for x in screening_rows if x["include_for_manual_case_review"] == "true"]
    summary = {
        "candidate_articles": len(screening_rows),
        "accessible_fulltext": sum(x["has_accessible_fulltext"] == "true" for x in screening_rows),
        "likely_case_reports": sum(x["likely_case_report"] == "true" for x in screening_rows),
        "manual_case_review_candidates": len(included),
        "direct_case_signal_candidates": sum(x["evidence_tier"] == "direct_case_signal" for x in screening_rows),
        "dechallenge_signal_candidates": sum(x["dechallenge_signal"] == "true" for x in screening_rows),
        "rechallenge_signal_candidates": sum(x["rechallenge_signal"] == "true" for x in screening_rows),
        "alternative_cause_signal_candidates": sum(x["alternative_cause_signal"] == "true" for x in screening_rows),
        "included_by_drug_outcome": {},
        "claim_boundary": "screening aid only; no automated case-level causality or clinical adjudication",
    }
    for row in included:
        key = f"{row['drug_class']}|{row['outcome_class']}"
        summary["included_by_drug_outcome"][key] = summary["included_by_drug_outcome"].get(key, 0) + 1
    (BASE / "screening_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
