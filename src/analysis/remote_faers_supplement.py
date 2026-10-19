#!/usr/bin/env python3
"""Server-only FAERS supplement driver.

Builds a compact report matrix from the already downloaded quarterly files and
runs independent supplementary analyses on that matrix. This file is never
run against local data.
"""

import argparse
import csv
import gzip
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)
import run_faers_full_window as fw


OUTCOMES = ("hepatic", "hematologic", "any_target_outcome", "serious")
LABELS = {1: "METHIMAZOLE", 2: "PROPYLTHIOURACIL"}
GROUPS = {0: "NEITHER", 1: "METHIMAZOLE", 2: "PROPYLTHIOURACIL", 3: "MIXED_TARGETS"}
TARGET_RE = re.compile(
    r"(?:METHIMAZOLE|THIAMAZOLE|TAPAZOLE|THYROZOL|THIAMAZOL|"
    r"PROPYLTHIOURACIL|PROP[- ]?THIOURACIL|PROPHYLTHIOURACIL|PTU|PROPYCIL|TIURAGYL)"
)


def cont(r):
    a, b, c, d = (x + 0.5 for x in r)
    value = a * d / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return value, math.exp(math.log(value) - 1.96 * se), math.exp(math.log(value) + 1.96 * se)


def prr(r):
    a, b, c, d = (x + 0.5 for x in r)
    value = (a / (a + b)) / (c / (c + d))
    se = math.sqrt(1 / a - 1 / (a + b) + 1 / c - 1 / (c + d))
    return value, math.exp(math.log(value) - 1.96 * se), math.exp(math.log(value) + 1.96 * se)


def files(root, quarter, prefix):
    return fw.files_for(root, quarter, prefix)


def build(root):
    quarters = fw.DEFAULT_QUARTERS
    results = os.path.join(root, "results")
    os.makedirs(results, exist_ok=True)
    latest = {}
    for quarter in quarters:
        for row in fw.rows(files(root, quarter, "DEMO")):
            if len(row) < 17 or not row[1]:
                continue
            version = int(row[2]) if row[2].isdigit() else 0
            old = latest.get(row[1])
            if old is None or (version, row[0]) > (old[0], old[1]):
                latest[row[1]] = [version, row[0], quarter, row[4], row[15], row[16]]
    reports = {
        item[1]: {
            "quarter": item[2], "event_dt": item[3], "age": item[4], "sex": item[5],
            "ps": 0, "psss": 0, "hepatic": 0, "hematologic": 0, "serious": 0,
        }
        for item in latest.values()
    }
    del latest

    flags = {}
    target_presence = {}
    for quarter in quarters:
        for row in fw.rows(files(root, quarter, "DRUG")):
            if len(row) < 6 or row[0] not in reports:
                continue
            kind = fw.target_kind(row)
            if kind:
                target_presence.setdefault(row[0], set()).update((1, 2) if kind == 3 else (kind,))
            if not kind or row[3] not in ("PS", "SS"):
                continue
            item = flags.setdefault(row[0], [0, 0, 0, 0])
            if kind in (1, 3):
                item[0] = 1
            if kind in (2, 3):
                item[1] = 1
            if row[3] == "PS" and kind in (1, 3):
                item[2] = 1
            if row[3] == "PS" and kind in (2, 3):
                item[3] = 1
    target_pids = set(flags)
    for pid, item in flags.items():
        any_m, any_p, ps_m, ps_p = item
        present = target_presence.get(pid, set())
        reports[pid]["ps"] = 3 if len(present) > 1 else 1 if ps_m and not ps_p else 2 if ps_p and not ps_m else 4 if present else 0
        reports[pid]["psss"] = 3 if len(present) > 1 else 1 if any_m and not any_p else 2 if any_p and not any_m else 4 if present else 0
    for pid, present in target_presence.items():
        if len(present) > 1:
            reports[pid]["ps"] = reports[pid]["psss"] = 3
        elif reports[pid]["ps"] == 0:
            reports[pid]["ps"] = reports[pid]["psss"] = 4
    del flags

    for quarter in quarters:
        for row in fw.rows(files(root, quarter, "REAC")):
            if len(row) < 3 or row[0] not in reports:
                continue
            outcome = fw.PRIMARY_PT.get(row[2].strip().upper())
            if outcome:
                reports[row[0]][outcome] = 1
        for row in fw.rows(files(root, quarter, "OUTC")):
            if len(row) >= 3 and row[0] in reports and row[2].strip().upper() in fw.SERIOUS_CODES:
                reports[row[0]]["serious"] = 1

    matrix = os.path.join(results, "analysis_matrix.tsv.gz")
    tmp = matrix + ".tmp"
    with gzip.open(tmp, "wt", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["primaryid", "quarter", "event_dt", "sex", "age", "ps", "psss", *OUTCOMES])
        for pid, item in reports.items():
            outcome_values = [
                item["hepatic"],
                item["hematologic"],
                int(item["hepatic"] or item["hematologic"]),
                item["serious"],
            ]
            writer.writerow([pid, item["quarter"], item["event_dt"], item["sex"], item["age"], item["ps"], item["psss"], *outcome_values])
    os.replace(tmp, matrix)

    co_meds = defaultdict(set)
    for quarter in quarters:
        for row in fw.rows(files(root, quarter, "DRUG")):
            if len(row) < 6 or row[0] not in target_pids or row[3] in ("PS", "SS"):
                continue
            name = row[4].strip().upper()
            if name and not TARGET_RE.search(name) and name not in {"UNKNOWN", "UNSPECIFIED"}:
                co_meds[row[0]].add(name)
    med_path = os.path.join(results, "target_concomitant.tsv.gz")
    tmp = med_path + ".tmp"
    with gzip.open(tmp, "wt", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["primaryid", "drug"])
        for pid, names in co_meds.items():
            for name in sorted(names):
                writer.writerow([pid, name])
    os.replace(tmp, med_path)
    summary = {
        "reports": len(reports),
        "target_ps_reports": sum(1 for x in reports.values() if x["ps"] in (1, 2)),
        "target_psss_reports": sum(1 for x in reports.values() if x["psss"] in (1, 2)),
        "matrix": matrix,
        "concomitant_rows": sum(len(x) for x in co_meds.values()),
        "quarters": len(quarters),
    }
    with open(os.path.join(results, "supplement_build_summary.json.tmp"), "w", encoding="utf-8") as out:
        json.dump(summary, out, indent=2)
    os.replace(os.path.join(results, "supplement_build_summary.json.tmp"), os.path.join(results, "supplement_build_summary.json"))


def read_matrix(root):
    path = os.path.join(root, "results", "analysis_matrix.tsv.gz")
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def strata(row):
    yield "ALL"
    if row["sex"] in ("F", "M"):
        yield "SEX_" + row["sex"]
    if row["age"] and row["age"] in "AENPT":
        yield "AGE_" + row["age"]


def ror_task(root):
    counts = defaultdict(lambda: [0, 0, 0, 0])
    for row in read_matrix(root):
        for exposure in ("ps", "psss"):
            group = int(row[exposure])
            if group not in (0, 1, 2):
                continue
            for stratum in strata(row):
                for outcome in OUTCOMES:
                    event = int(row[outcome])
                    for target in (1, 2):
                        if group not in (0, target):
                            continue
                        cell = counts[(exposure, stratum, outcome, target)]
                        if group == target:
                            cell[0 if event else 1] += 1
                        elif event:
                            cell[2] += 1
                        else:
                            cell[3] += 1
    path = os.path.join(root, "results", "supplement_ror_prr.tsv")
    with open(path + ".tmp", "w", encoding="utf-8", newline="") as out:
        fields = ["exposure", "stratum", "outcome", "target", "a", "b", "c", "d", "ROR", "ROR95_low", "ROR95_high", "PRR", "PRR95_low", "PRR95_high"]
        out.write("\t".join(fields) + "\n")
        for (exposure, stratum, outcome, group), cell in sorted(counts.items()):
            rv = cont(cell); pv = prr(cell)
            out.write("\t".join([exposure, stratum, outcome, LABELS[group], *map(str, cell), *[f"{x:.6f}" for x in (*rv, *pv)]]) + "\n")
    os.replace(path + ".tmp", path)


def ic_task(root):
    counts = defaultdict(lambda: [0, 0, 0, 0])
    for row in read_matrix(root):
        group = int(row["ps"])
        if group not in (0, 1, 2):
            continue
        for outcome in OUTCOMES:
            event = int(row[outcome])
            for target in (1, 2):
                if group not in (0, target):
                    continue
                cell = counts[(outcome, target)]
                if group == target:
                    cell[0 if event else 1] += 1
                elif event:
                    cell[2] += 1
                else:
                    cell[3] += 1
    path = os.path.join(root, "results", "supplement_ic.tsv")
    with open(path + ".tmp", "w", encoding="utf-8") as out:
        out.write("outcome\ttarget\ta\tb\tc\td\tIC\tIC025_approx\n")
        for (outcome, group), cell in sorted(counts.items()):
            a, b, c, d = cell
            n = a + b + c + d
            ic = math.log2(((a + 0.5) * (n + 1.0)) / ((a + b + 0.5) * (a + c + 0.5)))
            se = math.sqrt(1 / (a + 0.5) + 1 / (b + 0.5) + 1 / (c + 0.5) + 1 / (d + 0.5)) / math.log(2)
            out.write("\t".join([outcome, LABELS[group], *map(str, cell), f"{ic:.6f}", f"{ic - 1.96 * se:.6f}"]) + "\n")
    os.replace(path + ".tmp", path)


def temporal_task(root):
    counts = defaultdict(lambda: [0, 0, 0, 0])
    for row in read_matrix(root):
        group = int(row["ps"])
        if group not in (0, 1, 2):
            continue
        quarter = row["quarter"]
        for outcome in OUTCOMES:
            event = int(row[outcome])
            for target in (1, 2):
                if group not in (0, target):
                    continue
                cell = counts[(quarter, outcome, target)]
                if group == target:
                    cell[0 if event else 1] += 1
                elif event:
                    cell[2] += 1
                else:
                    cell[3] += 1
    path = os.path.join(root, "results", "supplement_quarterly_ror.tsv")
    with open(path + ".tmp", "w", encoding="utf-8") as out:
        out.write("quarter\toutcome\ttarget\ta\tb\tc\td\tROR\tCI95_low\tCI95_high\n")
        for (quarter, outcome, group), cell in sorted(counts.items()):
            value = cont(cell)
            out.write("\t".join([quarter, outcome, LABELS[group], *map(str, cell), *[f"{x:.6f}" for x in value]]) + "\n")
    os.replace(path + ".tmp", path)


def interaction_task(root):
    outcomes = {pid: {x: 0 for x in OUTCOMES} for pid in set()}
    group = {}
    for row in read_matrix(root):
        if int(row["ps"]) in (1, 2):
            group[row["primaryid"]] = int(row["ps"])
            outcomes[row["primaryid"]] = {x: int(row[x]) for x in OUTCOMES}
    meds = defaultdict(set)
    path = os.path.join(root, "results", "target_concomitant.tsv.gz")
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["primaryid"] in group:
                meds[row["drug"]].add(row["primaryid"])
    path = os.path.join(root, "results", "supplement_concomitant_ror.tsv")
    with open(path + ".tmp", "w", encoding="utf-8") as out:
        out.write("target\tconcomitant_drug\tn_reports\toutcome\ta\tb\tc\td\tROR\tCI95_low\tCI95_high\n")
        for drug, pids in sorted(meds.items(), key=lambda x: (-len(x[1]), x[0])):
            if len(pids) < 20:
                continue
            for target in (1, 2):
                selected = {pid for pid in pids if group.get(pid) == target}
                all_target = {pid for pid, value in group.items() if value == target}
                if len(selected) < 10 or len(all_target - selected) < 10:
                    continue
                for outcome in OUTCOMES:
                    a = sum(outcomes[pid][outcome] for pid in selected)
                    b = len(selected) - a
                    c = sum(outcomes[pid][outcome] for pid in all_target - selected)
                    d = len(all_target - selected) - c
                    value = cont((a, b, c, d))
                    out.write("\t".join([LABELS[target], drug, str(len(selected)), outcome, str(a), str(b), str(c), str(d), *[f"{x:.6f}" for x in value]]) + "\n")
    os.replace(path + ".tmp", path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root")
    parser.add_argument("task", choices=("build", "ror", "ic", "temporal", "interaction"))
    args = parser.parse_args()
    {"build": build, "ror": ror_task, "ic": ic_task, "temporal": temporal_task, "interaction": interaction_task}[args.task](args.root)


if __name__ == "__main__":
    main()
