#!/usr/bin/env python3
"""Build remote-only submission QC artifacts for the FAERS/JADER manuscript.

This script reads only aggregate QC and article-level screening artifacts. It
does not read FAERS/JADER raw files or report-level analysis matrices.
"""
from __future__ import annotations

import csv
import json
import os
from collections import Counter
from pathlib import Path

ROOT = Path(os.environ.get("FAERS_PROJECT_ROOT", "/root/autodl-tmp/faers-antithyroid-safety-signals"))
CASE = ROOT / "literature_case_evidence"
RESULTS = ROOT / "results"
OUT = ROOT / "submission_qc"
OUT.mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    screen = read_tsv(CASE / "screening_audit.tsv")
    curated = read_tsv(CASE / "curated_case_literature_summary.tsv")
    manual = read_tsv(CASE / "manual_case_review_sheet.tsv")
    final_qc = json.loads((RESULTS / "final_qc_report.json").read_text(encoding="utf-8"))

    included = [r for r in screen if r.get("include_for_manual_case_review") == "true"]
    direct = [r for r in included if r.get("evidence_tier") == "direct_case_signal"]
    pending = [r for r in manual if r.get("manual_verification_status") != "verified"]
    flow = [
        ("PubMed/Europe PMC candidate articles", len(screen), "search and metadata retrieval"),
        ("Manual case-review candidates", len(included), "heuristic screen; not final inclusion"),
        ("Direct-attribution screen candidates", len(direct), "automated attribution signal"),
        ("Curated publication-level evidence units", len(curated), "contextual summary"),
        ("Author-verified evidence units", len(manual) - len(pending), "required before submission"),
    ]
    with (OUT / "case_evidence_flow.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["stage", "n", "interpretation"])
        writer.writerows(flow)

    mapping = """# READUS-PV Submission Mapping\n\nStatus: BLOCKED pending author-level clinical verification and completion of the author declarations.\n\n| READUS-PV domain | Evidence artifact | Status | Remaining action |\n|---|---|---|---|\n| Background and objective | manuscript.tex, Introduction | present | retain reporting-signal wording |\n| Study design and data source | manuscript.tex, Methods 2.1 | present | add final access dates/release identifiers |\n| Study population and report unit | final_qc_report.json; full_window_summary.txt | present | retain database-specific duplicate rules |\n| Exposure definition | analysis scripts; manuscript.tex, Methods 2.3 | present | place full alias dictionary in supplement |\n| Outcome and serious-outcome definitions | manuscript.tex, Methods 2.4; OT audit outputs | present after correction | keep DE/LT/HO/DS/CA/OT and explicitly exclude RI |\n| Statistical analysis | aggregate ROR/PRR/IC outputs; Methods 2.5 | present | add exact software/version details |\n| Results flow and exclusions | case_evidence_flow.tsv; population table | partial | reconcile final author-verified case flow |\n| Disproportionality results | supplement_ror_prr.tsv; supplement_ic.tsv | present | retain denominators and uncertainty intervals |\n| Case-by-case clinical evidence | manual_case_review_sheet.tsv; curated summary | partial, critical | verify each included publication and mark inclusion/exclusion |\n| Interpretation and limitations | manuscript.tex, Discussion | present | add publication-selection bias after verification |\n| Data/code/declarations | manuscript.tex, Declarations | partial | complete funding, COI, CRediT, protocol and code statements |\n\nThe case-evidence screen is a review aid. Its counts are publication-level\nevidence units and are not independent patient counts, a systematic review, a\ncausality adjudication, or validation of FAERS/JADER.\n"""
    (OUT / "READUS_PV_MAPPING.md").write_text(mapping, encoding="utf-8")

    qc = {
        "status": "blocked_pending_author_case_verification",
        "raw_data_read": False,
        "report_level_matrix_read": False,
        "faers_final_qc_present": bool(final_qc),
        "case_screen": {
            "candidate_articles": len(screen),
            "manual_case_review_candidates": len(included),
            "direct_attribution_screen_candidates": len(direct),
            "curated_publication_units": len(curated),
            "manual_review_rows": len(manual),
            "manual_verification_pending": len(pending),
        },
        "blocking_items": [
            "author-level verification of each included publication-level evidence unit",
            "completion of READUS-PV checklist sign-off",
            "completion of author funding, conflict-of-interest, CRediT, protocol, and code statements",
            "final reference identity and DOI verification",
        ],
        "claim_boundary": "reporting signals and contextual published-case evidence only; no incidence, causality, pooled risk, or confirmed interaction",
    }
    (OUT / "submission_qc.json").write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
