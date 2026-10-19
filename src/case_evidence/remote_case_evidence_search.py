#!/usr/bin/env python3
"""Remote-only PubMed/PMC evidence search for antithyroid-drug case evidence.

This script is intended to run on the assigned server. It writes only article-level
metadata and downloaded full-text XML into the remote project root.
"""
from __future__ import annotations

import csv
import json
import os
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


ROOT = Path(os.environ.get("FAERS_PROJECT_ROOT", "/root/autodl-tmp/faers-antithyroid-safety-signals"))
OUT = ROOT / "literature_case_evidence"
FULLTEXT = OUT / "fulltext_xml"
OUT.mkdir(parents=True, exist_ok=True)
FULLTEXT.mkdir(parents=True, exist_ok=True)

QUERIES = {
    "hepatic_methimazole": (
        '(methimazole[Title/Abstract] OR thiamazole[Title/Abstract]) AND '
        '(hepatotoxicity[Title/Abstract] OR "liver injury"[Title/Abstract] OR '
        '"hepatic failure"[Title/Abstract] OR hepatitis[Title/Abstract]) AND '
        '(case reports[Publication Type] OR "case report"[Title] OR "case series"[Title])'
    ),
    "hepatic_ptu": (
        '(propylthiouracil[Title/Abstract] OR PTU[Title/Abstract]) AND '
        '(hepatotoxicity[Title/Abstract] OR "liver injury"[Title/Abstract] OR '
        '"hepatic failure"[Title/Abstract] OR hepatitis[Title/Abstract]) AND '
        '(case reports[Publication Type] OR "case report"[Title] OR "case series"[Title])'
    ),
    "hematologic_methimazole": (
        '(methimazole[Title/Abstract] OR thiamazole[Title/Abstract]) AND '
        '(agranulocytosis[Title/Abstract] OR neutropenia[Title/Abstract] OR '
        '"bone marrow"[Title/Abstract]) AND '
        '(case reports[Publication Type] OR "case report"[Title] OR "case series"[Title])'
    ),
    "hematologic_ptu": (
        '(propylthiouracil[Title/Abstract] OR PTU[Title/Abstract]) AND '
        '(agranulocytosis[Title/Abstract] OR neutropenia[Title/Abstract] OR '
        '"bone marrow"[Title/Abstract]) AND '
        '(case reports[Publication Type] OR "case report"[Title] OR "case series"[Title])'
    ),
    "broad_antithyroid_cases": (
        '(("antithyroid agents"[MeSH Terms] OR antithyroid[Title/Abstract] OR '
        'methimazole[Title/Abstract] OR propylthiouracil[Title/Abstract]) AND '
        '(hepatic[Title/Abstract] OR liver[Title/Abstract] OR agranulocytosis[Title/Abstract] '
        'OR neutropenia[Title/Abstract])) AND '
        '(case reports[Publication Type] OR "case report"[Title] OR "case series"[Title])'
    ),
}


def fetch(url: str, timeout: int = 30) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "FAERS-case-evidence-search/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def text(node: ET.Element | None) -> str:
    if node is None:
        return ""
    return " ".join("".join(node.itertext()).split())


def esearch(term: str) -> list[str]:
    query = urllib.parse.urlencode({"db": "pubmed", "term": term, "retmax": "1000", "retmode": "xml"})
    root = ET.fromstring(fetch("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + query))
    return [x.text for x in root.findall(".//Id") if x.text]


def parse_pubmed(payload: bytes) -> dict[str, dict]:
    root = ET.fromstring(payload)
    records: dict[str, dict] = {}
    for article in root.findall(".//PubmedArticle"):
        pmid = text(article.find(".//PMID"))
        if not pmid:
            continue
        title = text(article.find(".//ArticleTitle"))
        abstract = " ".join(text(x) for x in article.findall(".//AbstractText"))
        journal = text(article.find(".//Journal/Title"))
        year = text(article.find(".//PubDate/Year")) or text(article.find(".//PubDate/MedlineDate"))
        pubtypes = [text(x) for x in article.findall(".//PublicationType")]
        authors = []
        for author in article.findall(".//AuthorList/Author"):
            name = " ".join(filter(None, [text(author.find(".//LastName")), text(author.find(".//ForeName"))]))
            if name:
                authors.append(name)
        # Restrict identifiers to this article's PubmedData block. A broad
        # `.//ArticleId` search also sees identifiers in cited references.
        ids = {
            x.attrib.get("IdType", ""): text(x)
            for x in article.findall("./PubmedData/ArticleIdList/ArticleId")
        }
        mesh = [text(x.find("DescriptorName")) for x in article.findall(".//MeshHeading")]
        records[pmid] = {
            "pmid": pmid,
            "title": title,
            "abstract": abstract,
            "journal": journal,
            "year": year,
            "pubtypes": "; ".join(pubtypes),
            "authors": "; ".join(authors),
            "doi": ids.get("doi", ""),
            "pmcid": ids.get("pmc", ""),
            "mesh": "; ".join(mesh),
        }
    return records


def classify(record: dict) -> tuple[str, str]:
    blob = (record["title"] + " " + record["abstract"] + " " + record["mesh"]).lower()
    drugs = []
    if re.search(r"methimazole|thiamazole", blob):
        drugs.append("methimazole")
    if re.search(r"propylthiouracil|\bptu\b", blob):
        drugs.append("propylthiouracil")
    outcomes = []
    if re.search(r"hepat|liver|hepatic", blob):
        outcomes.append("hepatic")
    if re.search(r"agranulocyt|neutropeni|bone marrow", blob):
        outcomes.append("hematologic")
    return ";".join(drugs), ";".join(outcomes)


def download_fulltext(record: dict) -> str:
    pmcid = record.get("pmcid", "")
    if not pmcid:
        return ""
    path = FULLTEXT / f"{pmcid}.xml"
    if path.exists() and path.stat().st_size > 1000:
        return str(path)
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
    for attempt in range(3):
        try:
            payload = fetch(url, timeout=15)
            if b"<" not in payload[:200] or len(payload) < 1000:
                return ""
            tmp = path.with_suffix(".xml.tmp")
            tmp.write_bytes(payload)
            tmp.replace(path)
            return str(path)
        except Exception as exc:
            if attempt == 2:
                (OUT / "search_errors.log").open("a", encoding="utf-8").write(
                    f"fulltext={pmcid}\terror={type(exc).__name__}: {exc}\n"
                )
            time.sleep(0.2 * (attempt + 1))
    return ""


def main() -> None:
    protocol = {
        "database": "PubMed via NCBI E-utilities; full text via Europe PMC",
        "search_date_utc": time.strftime("%Y-%m-%d", time.gmtime()),
        "queries": QUERIES,
        "scope": "Published case reports and case series describing methimazole/thiamazole or propylthiouracil-associated hepatic or hematologic adverse events",
        "article_unit": "publication, with case-level extraction only when explicitly reported in accessible full text",
        "planned_exclusion": ["FAERS/JADER-only analyses", "animal or in-vitro studies", "reviews without extractable index cases", "articles without a target drug or prespecified outcome"],
        "claim_boundary": "clinical case evidence complements, but does not validate or establish causality for, spontaneous-report disproportionality signals",
    }
    (OUT / "search_protocol.json").write_text(json.dumps(protocol, indent=2), encoding="utf-8")

    counts = []
    all_pmids: set[str] = set()
    for label, term in QUERIES.items():
        try:
            ids = esearch(term)
        except Exception as exc:
            ids = []
            (OUT / "search_errors.log").open("a", encoding="utf-8").write(
                f"query={label}\terror={type(exc).__name__}: {exc}\n"
            )
        counts.append({"query_id": label, "count": len(ids), "pmids": ";".join(ids)})
        all_pmids.update(ids)
        with (OUT / "query_counts.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["query_id", "count", "pmids"], delimiter="\t")
            writer.writeheader()
            writer.writerows(counts)
        time.sleep(0.35)

    ids = sorted(all_pmids, key=int)
    records: dict[str, dict] = {}
    for offset in range(0, len(ids), 100):
        chunk = ids[offset:offset + 100]
        query = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(chunk), "retmode": "xml"})
        try:
            records.update(parse_pubmed(fetch("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + query)))
        except Exception as exc:
            (OUT / "search_errors.log").open("a", encoding="utf-8").write(
                f"efetch_offset={offset}\terror={type(exc).__name__}: {exc}\n"
            )
        time.sleep(0.35)

    for record in records.values():
        drug, outcome = classify(record)
        record["drug_class"] = drug
        record["outcome_class"] = outcome
    fields = ["pmid", "pmcid", "doi", "year", "journal", "title", "authors", "pubtypes", "drug_class", "outcome_class", "fulltext_xml", "abstract", "mesh"]
    with (OUT / "pubmed_candidates.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for record in records.values():
            record["fulltext_xml"] = ""
        writer.writerows(records[pmid] for pmid in ids if pmid in records)

    with ThreadPoolExecutor(max_workers=12) as pool:
        futures = {pool.submit(download_fulltext, record): record for record in records.values()}
        for future in as_completed(futures):
            futures[future]["fulltext_xml"] = future.result()

    with (OUT / "pubmed_candidates.tsv.tmp").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records[pmid] for pmid in ids if pmid in records)
    (OUT / "pubmed_candidates.tsv.tmp").replace(OUT / "pubmed_candidates.tsv")
    print(json.dumps({"status": "ok", "unique_pmids": len(records), "fulltext_xml": sum(bool(x["fulltext_xml"]) for x in records.values()), "out": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
