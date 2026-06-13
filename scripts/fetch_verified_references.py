"""Fetch a reproducible, topic-balanced PubMed reference library."""

from __future__ import annotations

import csv
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


OUT = Path("submission_package/ijms_submission_2026-06-13/references_verified.csv")
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "TIMP1-IJMS-manuscript/1.0"})
SESSION.mount(
    "https://",
    HTTPAdapter(
        max_retries=Retry(
            total=6,
            connect=6,
            read=6,
            backoff_factor=1.0,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=("GET",),
        )
    ),
)

QUERIES = [
    (
        "AKI_CKD_transition",
        '("acute kidney injury"[Title/Abstract] AND "chronic kidney disease"[Title/Abstract]) '
        'AND (maladaptive repair[Title/Abstract] OR transition[Title/Abstract]) AND 2015:2026[dp]',
        8,
    ),
    (
        "tubular_maladaptive_repair",
        '("kidney"[Title/Abstract] OR renal[Title/Abstract]) AND '
        '("maladaptive repair"[Title/Abstract] OR "failed repair"[Title/Abstract]) AND 2019:2026[dp]',
        8,
    ),
    (
        "kidney_single_cell_spatial",
        '("kidney"[Title/Abstract] OR renal[Title/Abstract]) AND '
        '("single-cell"[Title/Abstract] OR "single nucleus"[Title/Abstract] OR spatial[Title/Abstract]) '
        'AND (injury[Title/Abstract] OR fibrosis[Title/Abstract]) AND 2021:2026[dp]',
        8,
    ),
    (
        "kidney_fibrosis_ECM_TGFb",
        '("kidney fibrosis"[Title/Abstract] OR "renal fibrosis"[Title/Abstract]) AND '
        '("extracellular matrix"[Title/Abstract] OR TGF-beta[Title/Abstract] OR '
        '"TGF-β"[Title/Abstract]) AND 2020:2026[dp]',
        8,
    ),
    (
        "kidney_senescence",
        '("kidney"[Title/Abstract] OR renal[Title/Abstract]) AND '
        '(senescence[Title/Abstract] OR "cell cycle arrest"[Title/Abstract]) '
        'AND (injury[Title/Abstract] OR fibrosis[Title/Abstract]) AND 2019:2026[dp]',
        7,
    ),
    (
        "TIMP1_fibrosis_injury",
        '(TIMP1[Title/Abstract] OR "tissue inhibitor of metalloproteinases 1"[Title/Abstract]) '
        'AND (kidney[Title/Abstract] OR renal[Title/Abstract] OR fibrosis[Title/Abstract] '
        'OR injury[Title/Abstract]) AND 2015:2026[dp]',
        10,
    ),
    (
        "transcriptomic_reproducibility",
        '(transcriptomic[Title/Abstract] OR "gene expression"[Title/Abstract]) AND '
        '(kidney[Title/Abstract] OR renal[Title/Abstract]) AND '
        '(cohort[Title/Abstract] OR reproducib*[Title/Abstract] OR atlas[Title/Abstract]) '
        'AND 2020:2026[dp]',
        7,
    ),
]

def text(node: ET.Element | None) -> str:
    if node is None:
        return ""
    return "".join(node.itertext()).strip()


def search(term: str, retmax: int) -> list[str]:
    response = SESSION.get(
        f"{BASE}/esearch.fcgi",
        params={
            "db": "pubmed",
            "term": term,
            "retmode": "json",
            "retmax": retmax,
            "sort": "relevance",
        },
        timeout=60,
    )
    response.raise_for_status()
    return response.json()["esearchresult"]["idlist"]


def fetch(pmids: list[str]) -> list[ET.Element]:
    response = SESSION.get(
        f"{BASE}/efetch.fcgi",
        params={"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"},
        timeout=60,
    )
    response.raise_for_status()
    return ET.fromstring(response.content).findall(".//PubmedArticle")


def parse(article: ET.Element) -> dict[str, str]:
    citation = article.find("./MedlineCitation")
    article_node = citation.find("./Article")
    journal = article_node.find("./Journal")
    issue = journal.find("./JournalIssue")
    pub_date = issue.find("./PubDate")
    year = text(pub_date.find("./Year")) or text(pub_date.find("./MedlineDate"))[:4]
    authors = []
    for author in article_node.findall("./AuthorList/Author"):
        collective = text(author.find("./CollectiveName"))
        if collective:
            authors.append(collective)
            continue
        last = text(author.find("./LastName"))
        initials = text(author.find("./Initials"))
        if last:
            authors.append(f"{last} {initials}".strip())
    ids = {
        item.attrib.get("IdType", ""): text(item)
        for item in article.findall("./PubmedData/ArticleIdList/ArticleId")
    }
    return {
        "pmid": text(citation.find("./PMID")),
        "doi": ids.get("doi", ""),
        "year": year,
        "authors": "; ".join(authors),
        "title": text(article_node.find("./ArticleTitle")),
        "journal": text(journal.find("./Title")),
        "journal_abbrev": text(journal.find("./ISOAbbreviation")),
        "volume": text(issue.find("./Volume")),
        "issue": text(issue.find("./Issue")),
        "pages_or_elocation": text(article_node.find("./Pagination/MedlinePgn"))
        or text(article_node.find("./ELocationID")),
        "publication_types": "; ".join(
            text(x) for x in article_node.findall("./PublicationTypeList/PublicationType")
        ),
    }


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    topic_by_pmid: dict[str, str] = {}
    ordered_pmids: list[str] = []
    for topic, query, count in QUERIES:
        for pmid in search(query, count):
            if pmid not in topic_by_pmid:
                topic_by_pmid[pmid] = topic
                ordered_pmids.append(pmid)
        time.sleep(0.35)

    records = []
    for start in range(0, len(ordered_pmids), 50):
        for node in fetch(ordered_pmids[start : start + 50]):
            record = parse(node)
            record["topic"] = topic_by_pmid.get(record["pmid"], "")
            record["verified_source"] = (
                f"https://pubmed.ncbi.nlm.nih.gov/{record['pmid']}/"
            )
            records.append(record)
        time.sleep(0.35)

    by_pmid = {row["pmid"]: row for row in records}
    ordered = [by_pmid[pmid] for pmid in ordered_pmids if pmid in by_pmid]
    fields = [
        "topic",
        "pmid",
        "doi",
        "year",
        "authors",
        "title",
        "journal",
        "journal_abbrev",
        "volume",
        "issue",
        "pages_or_elocation",
        "publication_types",
        "verified_source",
    ]
    with OUT.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(ordered)
    print(f"Wrote {len(ordered)} verified references to {OUT}")


if __name__ == "__main__":
    main()
