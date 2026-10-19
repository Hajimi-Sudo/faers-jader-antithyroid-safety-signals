#!/usr/bin/env python3
"""Create a prespecified, conservative published-case evidence summary remotely.

This is a structured case-evidence analysis, not a systematic review and not a
causality adjudication. The selected publications are restricted to reports
whose abstract or accessible full text explicitly describes the target drug and
one prespecified hepatic or hematologic outcome. Multi-case publications are
kept as publications because the available evidence is heterogeneous.
"""
from __future__ import annotations

import csv
import json
import os
from collections import Counter
from pathlib import Path

ROOT = Path(os.environ.get("FAERS_PROJECT_ROOT", "/root/autodl-tmp/faers-antithyroid-safety-signals"))
BASE = ROOT / "literature_case_evidence"

# The selection was fixed before reading the completed FAERS estimates. Each
# row is a publication-level evidence unit, not an independent patient count.
SELECTED = {
    "522183": {"target_drug": "propylthiouracil", "outcome_group": "hepatic", "evidence_level": "direct_case_report", "alternative_cause": "autoimmune thyroiditis noted; viral/other causes not stated in abstract", "clinical_summary": "10-year-old girl; PTU 100 mg three times daily; jaundice and hepatitis after 40 days; improved after discontinuation"},
    "1261772": {"target_drug": "propylthiouracil", "outcome_group": "hepatic", "evidence_level": "direct_case_report", "alternative_cause": "not identified in abstract", "clinical_summary": "24-year-old woman; fulminant hepatic failure; submassive hepatic necrosis; clinical and biochemical recovery after discontinuation; PTU lymphocyte sensitization"},
    "1687933": {"target_drug": "methimazole", "outcome_group": "hepatic", "evidence_level": "case_series", "alternative_cause": "hepatitis A/B, alcohol, and other hepatotoxic drugs excluded in the abstract", "clinical_summary": "Three women; cholestatic jaundice after approximately 30 days; cholestatic biopsy pattern; eventual resolution after interruption"},
    "1994070": {"target_drug": "propylthiouracil", "outcome_group": "hepatic", "evidence_level": "direct_case_report", "alternative_cause": "not stated in abstract", "clinical_summary": "43-year-old woman; PTU-associated hepatitis; liver function deteriorated for two weeks after discontinuation"},
    "2131935": {"target_drug": "methimazole", "outcome_group": "hepatic", "evidence_level": "direct_case_report", "alternative_cause": "iatrogenic attribution suggested, not confirmed in abstract", "clinical_summary": "56-year-old woman; long-term methimazole; granulomatous hepatitis on liver biopsy"},
    "2343593": {"target_drug": "propylthiouracil", "outcome_group": "hepatic", "evidence_level": "direct_case_report", "alternative_cause": "not stated in abstract", "clinical_summary": "Child; PTU-associated hepatic failure and encephalopathy; report emphasizes prompt discontinuation after clinical or biochemical injury"},
    "3164375": {"target_drug": "methimazole", "outcome_group": "hematologic", "evidence_level": "direct_case_report", "alternative_cause": "infection symptoms present; neutropenia attributed to methimazole in report", "clinical_summary": "Agranulocytosis after two months; fever, sore throat, profound leukopenia, and oral lesions; oral findings abated after withdrawal"},
    "3528517": {"target_drug": "propylthiouracil", "outcome_group": "hepatic", "evidence_level": "direct_case_report", "alternative_cause": "not identified in abstract", "clinical_summary": "33-year-old woman; cholestatic jaundice two weeks after PTU; biopsy, rash, and positive migration inhibition factor supported the report attribution"},
    "3778829": {"target_drug": "propylthiouracil", "outcome_group": "hematologic", "evidence_level": "direct_case_report", "alternative_cause": "not identified in abstract", "clinical_summary": "Agranulocytosis and marrow hypoplasia after second PTU exposure; PTU-dependent antibody findings reported"},
    "3804181": {"target_drug": "methimazole", "outcome_group": "hepatic", "evidence_level": "direct_case_report", "alternative_cause": "not stated in abstract", "clinical_summary": "58-year-old woman; 20 mg/day methimazole; cholestatic injury after 20 days; rechallenge reproduced reaction and discontinuation improved symptoms"},
    "8605651": {"target_drug": "methimazole", "outcome_group": "hematologic", "evidence_level": "direct_case_report", "alternative_cause": "not identified in abstract", "clinical_summary": "53-year-old woman; 5 mg three times daily; leukopenia at six weeks followed by agranulocytosis; recovery after withdrawal and G-CSF"},
    "8876856": {"target_drug": "propylthiouracil", "outcome_group": "hepatic", "evidence_level": "direct_case_report", "alternative_cause": "alcohol/drug abuse, viral serologies, CMV/EBV, and autoimmune markers negative in abstract", "clinical_summary": "64-year-old woman; PTU 250 mg/day for one year; hepatitis and biopsy-confirmed chronic active hepatitis; cirrhosis during follow-up"},
    "9025882": {"target_drug": "propylthiouracil", "outcome_group": "hepatic", "evidence_level": "direct_case_report", "alternative_cause": "not stated in abstract", "clinical_summary": "13-year-old girl; jaundice, severe liver dysfunction, metabolic encephalopathy, and death despite transplantation"},
    "9086574": {"target_drug": "methimazole", "outcome_group": "hematologic", "evidence_level": "direct_case_report", "alternative_cause": "not stated in abstract", "clinical_summary": "Methimazole-induced aplastic anemia; failure of initial colony-stimulating factor and later hematopoietic recovery with immunosuppression"},
    "9252857": {"target_drug": "propylthiouracil", "outcome_group": "hepatic", "evidence_level": "direct_case_report", "alternative_cause": "not stated in abstract", "clinical_summary": "PTU-induced cholestatic hepatotoxicity one day after initiation; recovery occurred, with prolonged enzyme abnormalities"},
    "9773944": {"target_drug": "propylthiouracil", "outcome_group": "hepatic", "evidence_level": "direct_case_report", "alternative_cause": "not stated in abstract", "clinical_summary": "21-year-old woman; PTU after thiamazole rash; severe hepatitis after two months; positive drug-induced lymphocyte test and later normalization"},
    "9777751": {"target_drug": "methimazole", "outcome_group": "hematologic", "evidence_level": "direct_case_report", "alternative_cause": "not stated in abstract", "clinical_summary": "58-year-old woman; aplastic anemia during third methimazole exposure; withdrawal and G-CSF associated with recovery"},
    "10476084": {"target_drug": "methimazole", "outcome_group": "hematologic", "evidence_level": "direct_case_report", "alternative_cause": "cellulitis/infectious presentation treated concurrently", "clinical_summary": "47-year-old woman after nine years of methimazole; profound neutropenia with fever and cellulitis; recovery after withdrawal, antibiotics, and G-CSF"},
    "21234410": {"target_drug": "propylthiouracil", "outcome_group": "hepatic", "evidence_level": "case_series", "alternative_cause": "case-specific alternative causes assessed in source; heterogeneous", "clinical_summary": "Two PTU-associated acute liver-failure cases requiring transplantation"},
    "21929775": {"target_drug": "propylthiouracil", "outcome_group": "hepatic", "evidence_level": "direct_case_report", "alternative_cause": "not stated in abstract", "clinical_summary": "36-year-old pregnant woman; fulminant hepatic failure at 17 weeks after PTU; transplantation and thyroidectomy; pregnancy continued"},
    "23798915": {"target_drug": "propylthiouracil", "outcome_group": "hepatic", "evidence_level": "direct_case_report", "alternative_cause": "not stated in abstract", "clinical_summary": "26-year-old pregnant woman; acute liver failure during first trimester PTU; recovered without transplantation"},
    "25932342": {"target_drug": "methimazole", "outcome_group": "hematologic", "evidence_level": "direct_case_report", "alternative_cause": "infection symptoms present; report diagnosed methimazole-induced neutropenia", "clinical_summary": "31-year-old woman; neutropenia with fever, sore throat, and gingival necrosis; withdrawal, G-CSF, and supportive treatment led to recovery"},
    "29245333": {"target_drug": "methimazole", "outcome_group": "hepatic", "evidence_level": "direct_case_report", "alternative_cause": "multiple drugs and biliary disease considered; RUCAM used", "clinical_summary": "Elderly man; cholestatic jaundice after methimazole; liver enzymes improved after discontinuation"},
    "40007940": {"target_drug": "methimazole", "outcome_group": "hepatic", "evidence_level": "direct_case_report", "alternative_cause": "periampullary diverticulum/biliary obstruction considered and investigated", "clinical_summary": "72-year-old woman; cholestatic DILI attributed using updated RUCAM; symptoms and liver tests improved after discontinuation"},
}


def main() -> None:
    candidates = {}
    with (BASE / "pubmed_candidates.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            candidates[row["pmid"]] = row

    rows = []
    for pmid, spec in SELECTED.items():
        source = candidates.get(pmid)
        if source is None:
            raise RuntimeError(f"selected PMID not present in remote candidate table: {pmid}")
        rows.append({
            "pmid": pmid,
            "pmcid": source.get("pmcid", ""),
            "doi": source.get("doi", ""),
            "year": source.get("year", ""),
            "title": source.get("title", ""),
            "target_drug": spec["target_drug"],
            "outcome_group": spec["outcome_group"],
            "evidence_level": spec["evidence_level"],
            "fulltext_available": str(bool(source.get("fulltext_xml"))).lower(),
            "clinical_summary": spec["clinical_summary"],
            "alternative_cause_assessment": spec["alternative_cause"],
            "manual_status": "remote_fulltext_or_abstract_verified; author_clinical_review_still_required",
        })

    fields = list(rows[0])
    out = BASE / "curated_case_literature_summary.tsv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "status": "curated_published_case_evidence_summary_complete",
        "included_publications": len(rows),
        "publication_unit": True,
        "independent_patient_count_not_pooled": True,
        "by_drug_outcome": dict(Counter(f"{r['target_drug']}|{r['outcome_group']}" for r in rows)),
        "case_series_publications": sum(r["evidence_level"] == "case_series" for r in rows),
        "direct_case_report_publications": sum(r["evidence_level"] == "direct_case_report" for r in rows),
        "fulltext_available": sum(r["fulltext_available"] == "true" for r in rows),
        "with_alternative_cause_or_competing_explanation": sum("not stated" not in r["alternative_cause_assessment"] for r in rows),
        "claim_boundary": "structured published-case context; not a systematic review, causality adjudication, incidence estimate, or independent validation",
        "detail_artifact": str(out),
    }
    (BASE / "curated_case_literature_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
