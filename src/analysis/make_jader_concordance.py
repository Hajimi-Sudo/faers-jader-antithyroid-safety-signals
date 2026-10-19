#!/usr/bin/env python3
"""Create aggregate JADER PRR and FAERS-JADER pilot concordance tables."""

import csv
import math
import os


ROOT = os.environ.get("FAERS_PROJECT_ROOT", "/root/autodl-tmp/faers-antithyroid-safety-signals")
RESULTS = os.path.join(ROOT, "results")


def main():
    source = os.path.join(RESULTS, "jader_window_pilot_ror.tsv")
    prr_path = os.path.join(RESULTS, "jader_window_pilot_prr.tsv")
    rows = []
    with open(source, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames + ["PRR", "PRR95_low", "PRR95_high"]
        for row in reader:
            a, b, c, d = (float(row[key]) + 0.5 for key in ("a", "b", "c", "d"))
            value = (a / (a + b)) / (c / (c + d))
            se = math.sqrt(1 / a - 1 / (a + b) + 1 / c - 1 / (c + d))
            row.update(
                PRR=f"{value:.6f}",
                PRR95_low=f"{math.exp(math.log(value) - 1.96 * se):.6f}",
                PRR95_high=f"{math.exp(math.log(value) + 1.96 * se):.6f}",
            )
            rows.append(row)
    with open(prr_path + ".tmp", "w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(fields) + "\n")
        for row in rows:
            handle.write("\t".join(row[field] for field in fields) + "\n")
    os.replace(prr_path + ".tmp", prr_path)

    faers = {}
    with open(os.path.join(RESULTS, "full_window_ror.tsv"), encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["stratum"] == "ALL" and row["outcome"] in ("hepatic", "hematologic", "any_target_outcome"):
                faers[(row["target"], row["outcome"])] = row
    jader = {(row["target"], row["outcome"]): row for row in rows}
    output = os.path.join(RESULTS, "cross_database_pilot_concordance.tsv")
    fields = [
        "target", "outcome",
        "FAERS_a", "FAERS_b", "FAERS_c", "FAERS_d", "FAERS_ROR", "FAERS_CI95_low", "FAERS_CI95_high",
        "JADER_a", "JADER_b", "JADER_c", "JADER_d", "JADER_ROR", "JADER_CI95_low", "JADER_CI95_high",
        "JADER_PRR", "JADER_PRR95_low", "JADER_PRR95_high",
    ]
    with open(output + ".tmp", "w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(fields) + "\n")
        for key in sorted(jader):
            frow = faers[key]
            jrow = jader[key]
            values = [
                key[0], key[1],
                frow["a"], frow["b"], frow["c"], frow["d"], frow["ROR"], frow["CI95_low"], frow["CI95_high"],
                jrow["a"], jrow["b"], jrow["c"], jrow["d"], jrow["ROR"], jrow["CI95_low"], jrow["CI95_high"],
                jrow["PRR"], jrow["PRR95_low"], jrow["PRR95_high"],
            ]
            handle.write("\t".join(values) + "\n")
    os.replace(output + ".tmp", output)
    print(f"WROTE {prr_path}")
    print(f"WROTE {output}")


if __name__ == "__main__":
    main()
