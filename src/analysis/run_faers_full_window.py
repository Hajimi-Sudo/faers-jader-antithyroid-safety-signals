#!/usr/bin/env python3
"""Memory-bounded FAERS full-window analysis for the remote server only."""

import csv
import datetime as dt
import glob
import math
import os
import re
import sys


DEFAULT_QUARTERS = """
2014q4 2015q1 2015q2 2015q3 2015q4 2016q1 2016q2 2016q3 2016q4 2017q1
2017q2 2017q3 2017q4 2018q1 2018q2 2018q3 2018q4 2019q1 2019q2 2019q3
2019q4 2020q1 2020q2 2020q3 2020q4 2021q1 2021q2 2021q3 2021q4 2022q1
2022q2 2022q3 2022q4 2023q1 2023q2 2023q3 2023q4 2024q1 2024q2 2024q3
2024q4 2025q1
""".split()

PRIMARY_PT = {
    "ACUTE HEPATIC FAILURE": "hepatic",
    "HEPATIC FAILURE": "hepatic",
    "FULMINANT HEPATITIS": "hepatic",
    "TOXIC HEPATITIS": "hepatic",
    "DRUG-INDUCED LIVER INJURY": "hepatic",
    "LIVER INJURY": "hepatic",
    "HEPATOTOXICITY": "hepatic",
    "HEPATIC ENCEPHALOPATHY": "hepatic",
    "CHOLESTATIC HEPATITIS": "hepatic",
    "JAUNDICE": "hepatic",
    "AGRANULOCYTOSIS": "hematologic",
    "FEBRILE NEUTROPENIA": "hematologic",
    "NEUTROPENIA": "hematologic",
    "PANCYTOPENIA": "hematologic",
    "APLASTIC ANAEMIA": "hematologic",
    "APLASTIC ANEMIA": "hematologic",
    "BONE MARROW FAILURE": "hematologic",
    "MYELOSUPPRESSION": "hematologic",
}
# FDA FAERS OUTC codes for death, life-threatening, hospitalization, disability,
# congenital anomaly, and other serious/important medical event.
SERIOUS_CODES = {"DE", "LT", "HO", "DS", "CA", "OT"}
TARGET_LABELS = {1: "METHIMAZOLE", 2: "PROPYLTHIOURACIL"}
METHIMAZOLE_RE = re.compile(
    r"(?:^|[^A-Z])(?:METHIMAZOLE|THIAMAZOLE|TAPAZOLE|THYROZOL|THIAMAZOL|"
    r"THACAPZOL|THYCAPZOL|TIRODRIL|STRUMAZOL|UNIMAZOLE|METHIZOL|METIMAZOL|"
    r"METHEMAZOL|METIBASOL|THIMAZOLE)(?:[^A-Z]|$)"
)
PTU_RE = re.compile(
    r"(?:^|[^A-Z])(?:PROPYLTHIOURACIL|PROP[- ]?THIOURACIL|"
    r"PROPHYLTHIOURACIL|PTU|PROPYCIL|TIURAGYL)(?:[^A-Z]|$)"
)


def files_for(root, quarter, prefix):
    pattern = os.path.join(root, "raw", "work", quarter, "**", "*.txt")
    matches = []
    for path in glob.glob(pattern, recursive=True):
        name = os.path.basename(path)
        if re.fullmatch(prefix + r"\d+Q[1-4](?:_new)?\.txt", name, re.I):
            matches.append(path)
    if len(matches) != 1:
        raise RuntimeError(f"expected one {prefix} file for {quarter}, found {matches}")
    return matches[0]


def rows(path):
    with open(path, "r", encoding="latin1", newline="") as handle:
        reader = csv.reader(handle, delimiter="$", quoting=csv.QUOTE_NONE)
        for row in reader:
            if row and row[0].lower() != "primaryid":
                yield row


def target_kind(row):
    text = (row[4] + " " + row[5]).upper()
    methimazole = bool(METHIMAZOLE_RE.search(text))
    ptu = bool(PTU_RE.search(text))
    if methimazole and ptu:
        return 3
    if methimazole:
        return 1
    if ptu:
        return 2
    return 0


def date_days(start, event):
    if not start or not event or not re.fullmatch(r"\d{8}", start) or not re.fullmatch(r"\d{8}", event):
        return None
    try:
        first = dt.datetime.strptime(start, "%Y%m%d").date()
        second = dt.datetime.strptime(event, "%Y%m%d").date()
        return (second - first).days
    except ValueError:
        return None


def ror_line(reports, outcome_name, stratum, filter_index, filter_value, target):
    a = b = c = d = 0
    for report in reports.values():
        if filter_index is not None and report[filter_index] != filter_value:
            continue
        group = report[2]
        if group not in (target, 0):
            continue
        if outcome_name == "hepatic":
            event = bool(report[3])
        elif outcome_name == "hematologic":
            event = bool(report[4])
        elif outcome_name == "any_target_outcome":
            event = bool(report[3] or report[4])
        else:
            event = bool(report[5])
        if group == target:
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
    low = math.exp(math.log(ror) - 1.96 * se)
    high = math.exp(math.log(ror) + 1.96 * se)
    return "\t".join(
        [stratum, outcome_name, "ps_group", TARGET_LABELS[target], "NEITHER",
         str(a), str(b), str(c), str(d), f"{ror:.6f}", f"{low:.6f}", f"{high:.6f}"]
    )


def main(root):
    quarters = os.environ.get("FAERS_QUARTERS", " ".join(DEFAULT_QUARTERS)).replace(",", " ").split()
    progress_path = os.path.join(root, "results", "full_progress.log")
    os.makedirs(os.path.join(root, "results"), exist_ok=True)

    latest = {}
    stats = {"demo_rows": 0, "drug_rows": 0, "outc_rows": 0, "reac_rows": 0, "ther_rows": 0}
    repeated_caseids = 0
    extra_rows_removed = 0
    versioned_caseids = 0
    for quarter in quarters:
        with open(progress_path, "a", encoding="utf-8") as progress:
            progress.write(f"DEMO {quarter}\n")
        for row in rows(files_for(root, quarter, "DEMO")):
            if len(row) < 17 or not row[1]:
                continue
            stats["demo_rows"] += 1
            primaryid, caseid = row[0], row[1]
            version = int(row[2]) if row[2].isdigit() else 0
            old = latest.get(caseid)
            if old is None:
                latest[caseid] = [version, primaryid, row[4], row[13], row[15], row[16], False, False]
                continue
            extra_rows_removed += 1
            if not old[6]:
                repeated_caseids += 1
                old[6] = True
            if version != old[0] and not old[7]:
                versioned_caseids += 1
                old[7] = True
            if (version, primaryid) > (old[0], old[1]):
                changed = old[7]
                repeated = old[6]
                latest[caseid] = [version, primaryid, row[4], row[13], row[15], row[16], repeated, changed]

    reports = {}
    for version, primaryid, event_dt, age, age_group, sex, _, _ in latest.values():
        reports[primaryid] = [event_dt, age_group, 0, 0, 0, 0, sex, None, 0]
    del latest

    target_flags = {}
    target_presence = {}
    target_seq = {}
    for quarter in quarters:
        with open(progress_path, "a", encoding="utf-8") as progress:
            progress.write(f"DRUG {quarter}\n")
        for row in rows(files_for(root, quarter, "DRUG")):
            if len(row) < 6 or row[0] not in reports:
                continue
            stats["drug_rows"] += 1
            kind = target_kind(row)
            if kind:
                target_presence.setdefault(row[0], set()).update((1, 2) if kind == 3 else (kind,))
            role = row[3]
            if not kind or role not in ("PS", "SS"):
                continue
            flags = target_flags.setdefault(row[0], [0, 0, 0, 0])
            if kind in (1, 3):
                flags[0] = 1
            if kind in (2, 3):
                flags[1] = 1
            if role == "PS":
                if kind in (1, 3):
                    flags[2] = 1
                if kind in (2, 3):
                    flags[3] = 1
            if row[2]:
                target_seq.setdefault(row[0], {})[row[2]] = kind

    for primaryid, flags in target_flags.items():
        report = reports[primaryid]
        any_m, any_p, ps_m, ps_p = flags
        present = target_presence.get(primaryid, set())
        report[2] = 3 if len(present) > 1 else 1 if ps_m and not ps_p else 2 if ps_p and not ps_m else 4 if present else 0
        report[8] = 3 if len(present) > 1 else 1 if any_m and not any_p else 2 if any_p and not any_m else 4 if present else 0
    for primaryid, present in target_presence.items():
        if len(present) > 1:
            reports[primaryid][2] = reports[primaryid][8] = 3
        elif reports[primaryid][2] == 0:
            reports[primaryid][2] = reports[primaryid][8] = 4
    for quarter in quarters:
        with open(progress_path, "a", encoding="utf-8") as progress:
            progress.write(f"REAC_OUTC_THER {quarter}\n")
        for row in rows(files_for(root, quarter, "REAC")):
            if len(row) < 3 or row[0] not in reports:
                continue
            stats["reac_rows"] += 1
            outcome = PRIMARY_PT.get(row[2].strip().upper())
            if outcome == "hepatic":
                reports[row[0]][3] = 1
            elif outcome == "hematologic":
                reports[row[0]][4] = 1
        for row in rows(files_for(root, quarter, "OUTC")):
            if len(row) < 3 or row[0] not in reports:
                continue
            stats["outc_rows"] += 1
            if row[2].strip().upper() in SERIOUS_CODES:
                reports[row[0]][5] = 1
        for row in rows(files_for(root, quarter, "THER")):
            if len(row) < 5 or row[0] not in target_seq or row[2] not in target_seq[row[0]]:
                continue
            stats["ther_rows"] += 1
            start = row[3]
            old_start = reports[row[0]][7]
            if start and (old_start is None or start < old_start):
                reports[row[0]][7] = start

    out_tmp = os.path.join(root, "results", "full_window_ror.tsv.tmp")
    sum_tmp = os.path.join(root, "results", "full_window_summary.txt.tmp")
    with open(out_tmp, "w", encoding="utf-8") as out:
        out.write("stratum\toutcome\texposure_definition\ttarget\tbackground\ta\tb\tc\td\tROR\tCI95_low\tCI95_high\n")
        strata = [("ALL", None, None), ("SEX_F", 6, "F"), ("SEX_M", 6, "M")]
        strata += [(f"AGE_{age}", 1, age) for age in "AENPT"]
        for name, index, value in strata:
            for outcome_name in ("hepatic", "hematologic", "any_target_outcome", "serious"):
                for target in (1, 2):
                    out.write(ror_line(reports, outcome_name, name, index, value, target) + "\n")
    os.replace(out_tmp, os.path.join(root, "results", "full_window_ror.tsv"))

    group_counts = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}
    any_counts = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}
    valid_onset = 0
    for report in reports.values():
        group_counts[report[2]] += 1
        any_counts[report[8]] += 1
        if date_days(report[7], report[0]) is not None and date_days(report[7], report[0]) >= 0:
            valid_onset += 1
    labels = {0: "NEITHER_TARGET", 1: "METHIMAZOLE", 2: "PROPYLTHIOURACIL", 3: "MIXED_TARGETS", 4: "TARGET_NON_PRIMARY"}
    with open(sum_tmp, "w", encoding="utf-8") as out:
        out.write("FAERS antithyroid analysis\n")
        out.write(f"window={','.join(quarters)}\nquarters={len(quarters)}\n")
        out.write("report_unit=latest caseversion per caseid; primaryid tiebreak\nprimary_exposure=role_cod PS\n")
        out.write("sensitivity_exposure=role_cod PS or SS\nbackground=NEITHER_TARGET: no methimazole/thiamazole or propylthiouracil record in any drug role\nexcluded_from_comparison=TARGET_NON_PRIMARY and MIXED_TARGETS\n")
        out.write("outcome=explicit frozen PT list\nror=0.5 continuity correction; Wald 95% CI\n\nraw_rows\n")
        out.write("\n".join(f"{key}={stats[key]}" for key in sorted(stats)))
        out.write("\n\ncase_deduplication\n")
        out.write(f"latest_caseids={len(reports)}\nrepeated_caseids={repeated_caseids}\nextra_rows_removed={extra_rows_removed}\ncaseids_with_version_change={versioned_caseids}\n\nps_groups\n")
        out.write("\n".join(f"{labels[key]}={group_counts[key]}" for key in range(5)))
        out.write("\n\nps_or_ss_groups\n")
        out.write("\n".join(f"{labels[key]}={any_counts[key]}" for key in range(5)))
        out.write(f"\n\ntarget_alias_rows\nMETHIMAZOLE_rows={sum(1 for flags in target_flags.values() if flags[0])}\nPROPYLTHIOURACIL_rows={sum(1 for flags in target_flags.values() if flags[1])}\n")
        out.write(f"\nvalid_onset_days={valid_onset}\nror_file=results/full_window_ror.tsv\n")
    os.replace(sum_tmp, os.path.join(root, "results", "full_window_summary.txt"))
    with open(progress_path, "a", encoding="utf-8") as progress:
        progress.write("DONE\n")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/root/autodl-tmp/faers-antithyroid-safety-signals")
