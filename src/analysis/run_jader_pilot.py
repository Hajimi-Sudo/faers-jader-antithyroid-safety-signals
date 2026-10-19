#!/usr/bin/env python3
"""Remote-only JADER pilot for the antithyroid signal study.

The script is intended to run on the assigned server. It keeps only aggregate
outputs and uses the official Japanese generic-name field for exposure.
"""

import csv
import glob
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict


HEPATIC = {
    "\u6025\u6027\u809d\u4e0d\u5168",
    "\u809d\u4e0d\u5168",
    "\u5287\u75c7\u809d\u708e",
    "\u4e2d\u6bd2\u6027\u809d\u708e",
    "\u85ac\u7269\u6027\u809d\u969c\u5bb3",
    "\u809d\u969c\u5bb3",
    "\u809d\u6bd2\u6027",
    "\u809d\u6027\u8133\u75c7",
    "\u80c6\u6c41\u3046\u3063\u6ede\u6027\u809d\u708e",
    "\u9ec4\u75b8",
}
HEMATOLOGIC = {
    "\u7121\u9846\u7c92\u7403\u75c7",
    "\u767a\u71b1\u6027\u597d\u4e2d\u7403\u6e1b\u5c11\u75c7",
    "\u597d\u4e2d\u7403\u6e1b\u5c11\u75c7",
    "\u6c4e\u8840\u7403\u6e1b\u5c11\u75c7",
    "\u518d\u751f\u4e0d\u826f\u6027\u8ca7\u8840",
    "\u9aa8\u9ac4\u4e0d\u5168",
    "\u9aa8\u9ac4\u6291\u5236",
}

METHIMAZOLE_NAMES = {"\u30c1\u30a2\u30de\u30be\u30fc\u30eb", "THIAMAZOLE"}
PTU_NAMES = {"\u30d7\u30ed\u30d4\u30eb\u30c1\u30aa\u30a6\u30e9\u30b7\u30eb", "PROPYLTHIOURACIL"}
JAPANESE_QUARTERS = {"\u7b2c\u4e00": 1, "\u7b2c\u4e8c": 2, "\u7b2c\u4e09": 3, "\u7b2c\u56db": 4}


def csv_files(root, prefix):
    return sorted(glob.glob(os.path.join(root, "raw", "jader", "work", "part_*", f"{prefix}*.csv")))


def read_rows(path):
    with open(path, encoding="cp932", newline="") as handle:
        reader = csv.reader(handle)
        next(reader, None)
        yield from reader


def report_key(row):
    return row[0], row[1]


def drug_kind(name):
    value = name.strip()
    if value in METHIMAZOLE_NAMES:
        return 1
    if value in PTU_NAMES:
        return 2
    return 0


def latest_keys(keys):
    latest = {}
    for identifier, report_count in keys:
        try:
            count = int(report_count)
        except ValueError:
            count = -1
        old = latest.get(identifier)
        if old is None or count > old[0]:
            latest[identifier] = (count, report_count)
    return {(identifier, report_count) for identifier, (_, report_count) in latest.items()}


def period_index(value):
    match = re.fullmatch(r"(\d{4})\u30fb(\u7b2c\u4e00|\u7b2c\u4e8c|\u7b2c\u4e09|\u7b2c\u56db)", value.strip())
    if not match:
        return None
    return int(match.group(1)) * 4 + JAPANESE_QUARTERS[match.group(2)]


def ror(a, b, c, d):
    aa, bb, cc, dd = (x + 0.5 for x in (a, b, c, d))
    value = aa * dd / (bb * cc)
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return value, math.exp(math.log(value) - 1.96 * se), math.exp(math.log(value) + 1.96 * se)


def main(root):
    drug_paths = csv_files(root, "drug")
    demo_paths = csv_files(root, "demo")
    reac_paths = csv_files(root, "reac")
    if len(drug_paths) != 36 or len(demo_paths) != 36 or len(reac_paths) != 36:
        raise RuntimeError(f"expected 36 files per table, got drug={len(drug_paths)} demo={len(demo_paths)} reac={len(reac_paths)}")
    start_period = int(os.environ.get("JADER_PERIOD_START", "0"))
    end_period = int(os.environ.get("JADER_PERIOD_END", "999999"))

    # First freeze the report unit inside the prespecified time window. This
    # prevents a later report revision from leaking into an earlier window.
    latest_all = {}
    for path in demo_paths:
        for row in read_rows(path):
            if len(row) < 11:
                continue
            key = report_key(row)
            index = period_index(row[6])
            if index is None or not (start_period <= index <= end_period):
                continue
            try:
                count = int(row[1])
            except ValueError:
                count = -1
            old = latest_all.get(row[0])
            if old is not None and (count, key) < (old[0], old[1]):
                continue
            latest_all[row[0]] = (count, key, row)

    reports = {}
    for _, key, row in latest_all.values():
        reports[key] = {
                "identifier": row[0],
                "report_count": row[1],
                "sex": row[2],
                "age": row[3],
                "period": row[6],
                "status": row[7],
                "report_type": row[8],
                "e2b": row[10],
                "group": 0,
                "hepatic": 0,
                "hematologic": 0,
            }

    # Mark target exposure only after the windowed report unit is frozen.
    target_rows = 0
    target_keys = set()
    target_ps = defaultdict(set)
    target_presence = defaultdict(set)
    target_roles = Counter()
    generic_names = Counter()
    drug_rows = 0
    for path in drug_paths:
        for row in read_rows(path):
            drug_rows += 1
            if len(row) < 5:
                continue
            key = report_key(row)
            if key not in reports:
                continue
            kind = drug_kind(row[4])
            if not kind:
                continue
            target_rows += 1
            target_keys.add(key)
            target_presence[key].add(kind)
            target_roles[(kind, row[3])] += 1
            generic_names[(kind, row[4])] += 1
            if row[3] == "\u88ab\u7591\u85ac":
                target_ps[key].add(kind)

    for key, kinds in target_presence.items():
        ps_kinds = target_ps.get(key, set())
        reports[key]["group"] = 3 if len(kinds) > 1 else next(iter(ps_kinds)) if ps_kinds else 4

    latest_target_keys = target_keys
    duplicate_identifier_versions = Counter(identifier for identifier, _ in target_keys)

    group_counts = Counter(report["group"] for report in reports.values())
    missingness = Counter()
    for report in reports.values():
        for field in ("sex", "age", "period", "report_type", "e2b"):
            if not report[field]:
                missingness[field] += 1

    reaction_rows = 0
    target_reaction_rows = 0
    pt_counts = Counter()
    outcome_counts = Counter()
    for path in reac_paths:
        for row in read_rows(path):
            reaction_rows += 1
            key = report_key(row) if len(row) >= 2 else None
            if key not in reports or len(row) < 5:
                continue
            target_reaction_rows += 1
            pt = row[3].strip()
            pt_counts[pt] += 1
            outcome_counts[row[4].strip()] += 1
            if pt in HEPATIC:
                reports[key]["hepatic"] = 1
            if pt in HEMATOLOGIC:
                reports[key]["hematologic"] = 1

    output_dir = os.path.join(root, "results")
    os.makedirs(output_dir, exist_ok=True)
    output = os.path.join(output_dir, "jader_pilot_ror.tsv")
    with open(output + ".tmp", "w", encoding="utf-8", newline="") as handle:
        handle.write("target\toutcome\ta\tb\tc\td\tROR\tCI95_low\tCI95_high\n")
        for target, label in ((1, "METHIMAZOLE"), (2, "PROPYLTHIOURACIL")):
            for outcome in ("hepatic", "hematologic", "any_target_outcome"):
                a = b = c = d = 0
                for report in reports.values():
                    if report["group"] not in (0, target):
                        continue
                    event = bool(report["hepatic"] or report["hematologic"]) if outcome == "any_target_outcome" else bool(report[outcome])
                    if report["group"] == target:
                        if event:
                            a += 1
                        else:
                            b += 1
                    elif event:
                        c += 1
                    else:
                        d += 1
                value, low, high = ror(a, b, c, d)
                handle.write("\t".join((label, outcome, str(a), str(b), str(c), str(d), f"{value:.6f}", f"{low:.6f}", f"{high:.6f}")) + "\n")
    os.replace(output + ".tmp", output)

    summary = {
        "table_files": {"drug": len(drug_paths), "demo": len(demo_paths), "reac": len(reac_paths)},
        "drug_rows": drug_rows,
        "target_drug_rows_exact_generic_name": target_rows,
        "target_report_versions": len(target_keys),
        "latest_target_report_versions": len(latest_target_keys),
        "latest_target_demo_rows": len(reports),
        "target_identifier_count": len(duplicate_identifier_versions),
        "identifiers_with_multiple_report_counts": sum(1 for value in duplicate_identifier_versions.values() if value > 1),
        "roles": {f"{kind}:{role}": count for (kind, role), count in target_roles.items()},
        "generic_names": {f"{kind}:{name}": count for (kind, name), count in generic_names.items()},
        "reaction_rows": reaction_rows,
        "latest_report_rows_in_window": len(reports),
        "latest_report_reaction_rows": target_reaction_rows,
        "latest_group_counts": dict(group_counts),
        "missingness_in_latest_reports": dict(missingness),
        "jader_period_start": start_period,
        "jader_period_end": end_period,
        "target_pt_top50": pt_counts.most_common(50),
        "target_outcomes": outcome_counts.most_common(),
        "outcome_dictionary": {"hepatic": sorted(HEPATIC), "hematologic": sorted(HEMATOLOGIC)},
        "report_unit": "latest report_count per identifier; report-count duplicates remain an audited limitation",
        "background": "all latest JADER report versions with no methimazole/thiamazole or propylthiouracil record in any drug role; reports with target drugs in non-PS roles are excluded",
        "claim_boundary": "JADER disproportional reporting signal only; no incidence or causality",
    }
    with open(os.path.join(output_dir, "jader_pilot_summary.json.tmp"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    os.replace(os.path.join(output_dir, "jader_pilot_summary.json.tmp"), os.path.join(output_dir, "jader_pilot_summary.json"))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/root/autodl-tmp/faers-antithyroid-safety-signals")
