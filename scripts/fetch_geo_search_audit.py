"""Record the reproducible NCBI GEO search and all returned accessions."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import requests


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "timp1_validation" / "geo_search_hits_2026-06-14.csv"
API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
QUERY = (
    '(human[Organism]) AND (kidney OR renal) AND '
    '(AKI OR "acute kidney injury" OR CKD OR fibrosis) AND '
    '("expression profiling by array"[DataSet Type] OR '
    '"expression profiling by high throughput sequencing"[DataSet Type])'
)


def main() -> None:
    search = requests.get(
        f"{API}/esearch.fcgi",
        params={
            "db": "gds",
            "term": QUERY,
            "retmax": 500,
            "retmode": "json",
        },
        timeout=120,
    )
    search.raise_for_status()
    ids = search.json()["esearchresult"]["idlist"]
    rows: list[dict[str, str]] = []
    for start in range(0, len(ids), 100):
        batch = ids[start : start + 100]
        response = requests.get(
            f"{API}/esummary.fcgi",
            params={
                "db": "gds",
                "id": ",".join(batch),
                "retmode": "json",
            },
            timeout=120,
        )
        response.raise_for_status()
        result = response.json()["result"]
        for uid in batch:
            item = result.get(uid, {})
            accession = item.get("accession", "")
            rows.append(
                {
                    "entrez_uid": uid,
                    "accession": accession,
                    "title": item.get("title", ""),
                    "gds_type": item.get("gdstype", ""),
                    "sample_count": str(item.get("n_samples", "")),
                    "organism": ";".join(item.get("taxon", []) or []),
                    "search_date": "2026-06-14",
                    "query": QUERY,
                    "screening_status": (
                        "Included in analysis"
                        if accession in {"GSE139061", "GSE30718", "GSE66494", "GSE180394"}
                        else "Search hit; not included after title/summary and eligibility screening"
                    ),
                }
            )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"count": len(rows), "output": str(OUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
