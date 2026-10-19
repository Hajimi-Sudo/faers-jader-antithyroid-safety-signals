#!/usr/bin/env python3
"""Strict remote-only review of published antithyroid-drug case reports.

The output is deliberately split into a detailed server-side audit table and a
small aggregate JSON file.  This is a full-text evidence screen, not an
automated causality score and not a replacement for author clinical review.
"""
from __future__ import annotations

import csv
import json
import os
import re
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

ROOT = Path(os.environ.get("FAERS_PROJECT_ROOT", "/root/autodl-tmp/faers-antithyroid-safety-signals"))
BASE = ROOT / "literature_case_evidence"
XML = BASE / "fulltext_xml"

DRUGS = {
    "methimazole": re.compile(r"\b(?:methimazole|thiamazole)\b", re.I),
    "propylthiouracil": re.compile(r"\b(?:propylthiouracil|PTU)\b", re.I),
}
OUTCOMES = {
    "hepatic": re.compile(r"\b(?:hepatotoxicity|hepatitis|hepatic failure|liver injury|liver failure|jaundice|cholestasis|cholestatic|transaminase|bilirubin)\b", re.I),
    "hematologic": re.compile(r"\b(?:agranulocytosis|neutropenia|pancytopenia|bone marrow failure|myelosuppression|granulocytopenia|leukopenia)\b", re.I),
}
ATTRIBUTION = re.compile(r"\b(?:induced|drug[- ]induced|associated with|related to|caused by|secondary to|toxicity|adverse effect|side effect|likely induced|most plausible causative)\b", re.I)
CLINICAL = re.compile(r"\b(?:case presentation|case report|patient|admitted|presented|laboratory|alanine aminotransferase|absolute neutrophil|bilirubin|white blood cell|liver biopsy|bone marrow|hospital)\b", re.I)
DECHALLENGE = re.compile(r"\b(?:discontinu|withdraw|stop(?:ped|ping)?|cessat|after stopping|after discontinuing|upon discontinuation|recovered|improv|resolv|normaliz)\w*\b", re.I)
RECHALLENGE = re.compile(r"\b(?:rechallenge|re-exposure|reexposure|reintroduced|readminister|second exposure)\b", re.I)
ALTERNATIVE = re.compile(r"\b(?:thyroid storm|thyrotoxicosis|viral hepatitis|hepatitis [abce]|autoimmune hepatitis|sepsis|bacteremia|infection-associated|ischemia|shock liver|multiple myeloma|hemophagocytosis|other cause)\b", re.I)


def clean(value: str) -> str:
    return " ".join((value or "").split())


def xml_text(path: Path) -> str:
    try:
        root = ET.parse(path).getroot()
    except Exception:
        return ""
    # Exclude references, affiliations, licensing metadata, and cited cases.
    body = root.find(".//body")
    if body is None:
        return ""
    return clean(" ".join(body.itertext()))


def sentences(text: str) -> list[str]:
    return [clean(x) for x in re.split(r"(?<=[.!?])\s+", text) if len(clean(x)) > 20]


def relevant_sentences(text: str, drug: str, outcome: str) -> list[str]:
    out = []
    for sent in sentences(text):
        if DRUGS[drug].search(sent) and OUTCOMES[outcome].search(sent):
            out.append(sent[:1200])
    return out


def first_match(pattern: re.Pattern[str], text: str) -> str:
    match = pattern.search(text)
    return clean(text[max(0, match.start() - 180):match.end() + 420]) if match else ""


def assess(row: dict) -> list[dict]:
    path = XML / f"{row.get('pmcid', '')}.xml"
    full = xml_text(path)
    source = clean(" ".join([row.get("title", ""), row.get("abstract", ""), full]))
    records = []
    for drug in DRUGS:
        drug_text = [s for s in sentences(source) if DRUGS[drug].search(s)]
        if not drug_text:
            continue
        for outcome in OUTCOMES:
            linked = relevant_sentences(source, drug, outcome)
            direct = any(ATTRIBUTION.search(s) for s in linked)
            clinical = any(CLINICAL.search(s) for s in linked)
            course = bool(DECHALLENGE.search(source) or RECHALLENGE.search(source))
            alternative = bool(ALTERNATIVE.search(source))
            title_direct = bool(ATTRIBUTION.search(row.get("title", "")) and OUTCOMES[outcome].search(row.get("title", "")) and DRUGS[drug].search(row.get("title", "")))
            # Require drug-event linkage plus a clinical case signal.  A paper
            # that merely cites the other target drug is not assigned to it.
            high_confidence = bool(linked and (direct or title_direct) and clinical)
            if high_confidence and course and alternative:
                tier = "direct_case_with_alternative_cause_signal"
            elif high_confidence:
                tier = "direct_case_signal"
            elif linked:
                tier = "case_context_or_weak_linkage"
            else:
                continue
            records.append({
                "pmid": row.get("pmid", ""),
                "pmcid": row.get("pmcid", ""),
                "doi": row.get("doi", ""),
                "year": row.get("year", ""),
                "title": row.get("title", ""),
                "target_drug": drug,
                "outcome_group": outcome,
                "direct_attribution": str(bool(direct or title_direct)).lower(),
                "clinical_case_signal": str(clinical).lower(),
                "dechallenge_or_rechallenge_signal": str(course).lower(),
                "alternative_cause_signal": str(alternative).lower(),
                "evidence_tier": tier,
                "drug_event_evidence": " || ".join(linked[:3]),
                "age_sex_evidence": first_match(re.compile(r"\b(?:\d{1,3}[- ]year[- ]old|aged?\s+\d{1,3}|\d{1,3}\s*(?:yo|y/o)|male|female|woman|man|girl|boy)\b", re.I), source),
                "dose_latency_evidence": first_match(re.compile(r"\b(?:\d+(?:\.\d+)?\s*(?:mg|g|mcg|μg)|after\s+\d+\s+(?:day|week|month|year)s?|within\s+\d+\s+(?:day|week|month|year)s?)\b", re.I), source),
                "course_evidence": first_match(DECHALLENGE, source) or first_match(RECHALLENGE, source),
                "alternative_evidence": first_match(ALTERNATIVE, source),
            })
    return records


def main() -> None:
    source = BASE / "pubmed_candidates.tsv"
    rows = []
    with source.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("pmcid") and (XML / f"{row['pmcid']}.xml").exists():
                rows.extend(assess(row))

    detail = BASE / "strict_case_evidence.tsv"
    fields = list(rows[0]) if rows else []
    with detail.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    direct = [r for r in rows if r["evidence_tier"].startswith("direct_case")]
    high = [r for r in direct if r["evidence_tier"] == "direct_case_signal"]
    aggregate = {
        "status": "full_text_case_evidence_screen_complete",
        "candidate_articles": len({r["pmid"] for r in rows}),
        "target_event_pairs": len(rows),
        "direct_case_pairs": len(direct),
        "direct_case_pairs_without_alternative_cause_flag": len(high),
        "direct_case_pairs_with_alternative_cause_flag": len(direct) - len(high),
        "by_drug_outcome": dict(Counter(f"{r['target_drug']}|{r['outcome_group']}" for r in direct)),
        "by_evidence_tier": dict(Counter(r["evidence_tier"] for r in rows)),
        "with_dechallenge_or_rechallenge": sum(r["dechallenge_or_rechallenge_signal"] == "true" for r in direct),
        "with_alternative_cause_signal": sum(r["alternative_cause_signal"] == "true" for r in direct),
        "claim_boundary": "full-text case-evidence screen; not a causality score and not a substitute for author clinical adjudication",
        "detail_artifact": str(detail),
    }
    (BASE / "strict_case_evidence_aggregate.json").write_text(json.dumps(aggregate, indent=2), encoding="utf-8")
    print(json.dumps(aggregate, indent=2))


if __name__ == "__main__":
    main()
