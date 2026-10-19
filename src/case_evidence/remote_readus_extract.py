#!/usr/bin/env python3
"""Extract readable READUS-PV statement text on the remote server."""
from pathlib import Path
import os
import re
import xml.etree.ElementTree as ET

root = Path(os.environ.get("FAERS_PROJECT_ROOT", "/root/autodl-tmp/faers-antithyroid-safety-signals"))
path = root / "literature_case_evidence" / "readus_statement.xml"
out = root / "literature_case_evidence" / "readus_statement.txt"
tree = ET.parse(path)
text = " ".join(" ".join(tree.getroot().itertext()).split())
text = re.sub(r"\s+", " ", text)
out.write_text(text, encoding="utf-8")
print(f"wrote {len(text)} characters to {out}")
for key in ("checklist", "reporting", "disproportionality", "case reports", "limitations", "data source"):
    print(f"\n### {key}")
    for match in list(re.finditer(key, text, re.I))[:5]:
        print(text[max(0, match.start() - 300):match.end() + 700])
