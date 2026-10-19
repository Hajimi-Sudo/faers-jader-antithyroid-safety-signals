#!/usr/bin/env python3
"""Remote-only Gate 1 sensitivity and calendar-stability analysis."""

import csv
import datetime as dt
import glob
import math
import os
import re
import sys
from collections import defaultdict

import run_faers_full_window as fw


PERIODS = (
    ("P1_2014Q4_2017Q4", 2014, 2017),
    ("P2_2018Q1_2021Q4", 2018, 2021),
    ("P3_2022Q1_2025Q1", 2022, 2025),
)


def prr_line(row):
    a, b, c, d = (float(row[key]) + 0.5 for key in ("a", "b", "c", "d"))
    prr = (a / (a + b)) / (c / (c + d))
    se = math.sqrt(1 / a - 1 / (a + b) + 1 / c - 1 / (c + d))
    low = math.exp(math.log(prr) - 1.96 * se)
    high = math.exp(math.log(prr) + 1.96 * se)
    return prr, low, high


def run_prr(root):
    source = os.path.join(root, "results", "full_window_ror.tsv")
    output = os.path.join(root, "results", "gate1_prr.tsv")
    with open(source, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        with open(output + ".tmp", "w", encoding="utf-8", newline="") as out:
            fields = reader.fieldnames + ["PRR", "PRR95_low", "PRR95_high"]
            out.write("\t".join(fields) + "\n")
            for row in reader:
                prr, low, high = prr_line(row)
                row.update(PRR=f"{prr:.6f}", PRR95_low=f"{low:.6f}", PRR95_high=f"{high:.6f}")
                out.write("\t".join(row[field] for field in fields) + "\n")
    os.replace(output + ".tmp", output)


def period_for_date(value):
    if not value or not re.fullmatch(r"\d{8}", value):
        return None
    year = int(value[:4])
    for label, start, end in PERIODS:
        if start <= year <= end:
            return label
    return None


def load_reports(root):
    quarters = fw.DEFAULT_QUARTERS
    latest = {}
    stats = defaultdict(int)
    for quarter in quarters:
        for row in fw.rows(fw.files_for(root, quarter, "DEMO")):
            if len(row) < 17 or not row[1]:
                continue
            stats["demo_rows"] += 1
            version = int(row[2]) if row[2].isdigit() else 0
            old = latest.get(row[1])
            if old is None or (version, row[0]) > (old[0], old[1]):
                latest[row[1]] = [version, row[0], row[4], row[15], row[16]]
    reports = {
        item[1]: {
            "period": period_for_date(item[2]),
            "age_group": item[3],
            "sex": item[4],
            "group": 0,
            "hepatic": 0,
            "hematologic": 0,
            "serious": 0,
        }
        for item in latest.values()
    }
    del latest

    target_flags = {}
    target_presence = {}
    target_seq = {}
    for quarter in quarters:
        for row in fw.rows(fw.files_for(root, quarter, "DRUG")):
            if len(row) < 6 or row[0] not in reports:
                continue
            stats["drug_rows"] += 1
            kind = fw.target_kind(row)
            if kind:
                target_presence.setdefault(row[0], set()).update((1, 2) if kind == 3 else (kind,))
            if not kind or row[3] not in ("PS", "SS"):
                continue
            flags = target_flags.setdefault(row[0], [0, 0, 0, 0])
            if kind in (1, 3):
                flags[0] = 1
            if kind in (2, 3):
                flags[1] = 1
            if row[3] == "PS":
                if kind in (1, 3):
                    flags[2] = 1
                if kind in (2, 3):
                    flags[3] = 1
            if row[2]:
                target_seq.setdefault(row[0], set()).add(row[2])
    for primaryid, flags in target_flags.items():
        any_m, any_p, ps_m, ps_p = flags
        present = target_presence.get(primaryid, set())
        reports[primaryid]["group"] = 3 if len(present) > 1 else 1 if ps_m and not ps_p else 2 if ps_p and not ps_m else 4 if present else 0
    for primaryid, present in target_presence.items():
        if len(present) > 1:
            reports[primaryid]["group"] = 3
        elif reports[primaryid]["group"] == 0:
            reports[primaryid]["group"] = 4

    for quarter in quarters:
        for row in fw.rows(fw.files_for(root, quarter, "REAC")):
            if len(row) < 3 or row[0] not in reports:
                continue
            outcome = fw.PRIMARY_PT.get(row[2].strip().upper())
            if outcome:
                reports[row[0]][outcome] = 1
        for row in fw.rows(fw.files_for(root, quarter, "OUTC")):
            if len(row) >= 3 and row[0] in reports and row[2].strip().upper() in fw.SERIOUS_CODES:
                reports[row[0]]["serious"] = 1
    return reports


def period_ror(reports, period, outcome, target):
    a = b = c = d = 0
    for report in reports.values():
        if report["period"] != period or report["group"] not in (0, target):
            continue
        event = bool(report[outcome]) if outcome != "any_target_outcome" else bool(report["hepatic"] or report["hematologic"])
        if report["group"] == target:
            if event:
                a += 1
            else:
                b += 1
        elif event:
            c += 1
        else:
            d += 1
    aa, bb, cc, dd = (x + 0.5 for x in (a, b, c, d))
    ror = aa * dd / (bb * cc)
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return a, b, c, d, ror, math.exp(math.log(ror) - 1.96 * se), math.exp(math.log(ror) + 1.96 * se)


def run_calendar(root):
    reports = load_reports(root)
    output = os.path.join(root, "results", "gate1_calendar_ror.tsv")
    with open(output + ".tmp", "w", encoding="utf-8") as out:
        out.write("period\toutcome\ttarget\ta\tb\tc\td\tROR\tCI95_low\tCI95_high\n")
        for period, _, _ in PERIODS:
            for outcome in ("hepatic", "hematologic", "any_target_outcome", "serious"):
                for target, label in ((1, "METHIMAZOLE"), (2, "PROPYLTHIOURACIL")):
                    values = period_ror(reports, period, outcome, target)
                    a, b, c, d, ror, low, high = values
                    out.write("\t".join((period, outcome, label, str(a), str(b), str(c), str(d), f"{ror:.6f}", f"{low:.6f}", f"{high:.6f}")) + "\n")
    os.replace(output + ".tmp", output)
    with open(os.path.join(root, "results", "gate1_calendar_summary.txt"), "w", encoding="utf-8") as out:
        out.write("Calendar stability design: event-year periods fixed before inspection\n")
        out.write("P1=2014-2017; P2=2018-2021; P3=2022-2025\n")
        out.write(f"reports={len(reports)}\n")
        out.write("Sparse strata with total target/background event cells <3 are descriptive only.\n")


def main(root):
    run_prr(root)
    run_calendar(root)


if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    main(sys.argv[1] if len(sys.argv) > 1 else "/root/autodl-tmp/faers-antithyroid-safety-signals")
