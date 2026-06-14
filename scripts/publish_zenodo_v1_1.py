"""Create and publish Zenodo version 1.1.0 from the existing v1.0.0 record."""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from requests.exceptions import RequestException
from urllib3.util.retry import Retry


ROOT = Path(__file__).resolve().parents[1]
PARTS_DIR = ROOT / "public_repository" / "zenodo_upload_parts_v1.1.0"
STATE = ROOT / "public_repository" / "zenodo_v1.1.0_state.json"
RESULT = ROOT / "public_repository" / "zenodo_v1.1.0_published_record.json"
API = "https://zenodo.org/api"
PREVIOUS_RECORD_ID = 20680932


def session(token: str) -> requests.Session:
    client = requests.Session()
    client.headers.update({"Authorization": f"Bearer {token}"})
    client.mount(
        "https://",
        HTTPAdapter(
            max_retries=Retry(
                total=6,
                connect=6,
                read=6,
                backoff_factor=2,
                status_forcelist=(429, 500, 502, 503, 504),
                allowed_methods=("GET", "POST", "PUT", "DELETE"),
            )
        ),
    )
    return client


def check(response: requests.Response, action: str) -> dict:
    if not response.ok:
        raise RuntimeError(f"{action} failed: HTTP {response.status_code}\n{response.text[:2000]}")
    return response.json() if response.content else {}


def upload(client: requests.Session, url: str, path: Path) -> None:
    last_error: Exception | None = None
    for attempt in range(1, 9):
        try:
            with path.open("rb") as handle:
                response = client.put(
                    url,
                    data=handle,
                    headers={
                        "Content-Type": "application/octet-stream",
                        "Content-Length": str(path.stat().st_size),
                    },
                    timeout=(90, 7200),
                )
            check(response, f"Upload {path.name}")
            return
        except RequestException as error:
            last_error = error
            wait = min(15 * attempt, 90)
            print(f"Upload attempt {attempt}/8 interrupted: {type(error).__name__}; retrying in {wait}s", flush=True)
            time.sleep(wait)
    raise RuntimeError(f"Upload failed after 8 attempts: {last_error}")


def main() -> None:
    token = os.environ.get("ZENODO_TOKEN", "").strip()
    if not token:
        sys.exit("ZENODO_TOKEN is not set.")
    upload_files = sorted(PARTS_DIR.glob("TIMP1_AKI_CKD_transcriptomics_zenodo_v1.1.0.zip.part*"))
    upload_files += [PARTS_DIR / "PARTS_MANIFEST.csv", PARTS_DIR / "README_REASSEMBLY.txt"]
    upload_files = [path for path in upload_files if path.exists()]
    if len(upload_files) < 3:
        sys.exit("v1.1.0 upload parts are missing. Run prepare_zenodo_v1_1.py first.")

    client = session(token)
    state = json.loads(STATE.read_text(encoding="utf-8")) if STATE.exists() else {}

    if not state.get("deposition_id"):
        previous = check(
            client.get(f"{API}/deposit/depositions/{PREVIOUS_RECORD_ID}", timeout=120),
            "Validate access to previous Zenodo record",
        )
        version_response = check(
            client.post(previous["links"]["newversion"], timeout=180),
            "Create Zenodo new version",
        )
        draft_url = version_response.get("links", {}).get("latest_draft")
        draft = check(client.get(draft_url, timeout=120), "Load new-version draft") if draft_url else version_response
        state = {
            "deposition_id": draft["id"],
            "bucket_url": draft["links"]["bucket"],
            "doi": draft.get("metadata", {}).get("prereserve_doi", {}).get("doi", ""),
            "uploaded_parts": [],
            "files_cleared": False,
            "published": False,
        }
        STATE.write_text(json.dumps(state, indent=2), encoding="utf-8")
        print(f"Created version draft {state['deposition_id']}; reserved DOI {state['doi']}", flush=True)

    draft = check(
        client.get(f"{API}/deposit/depositions/{state['deposition_id']}", timeout=120),
        "Load version draft",
    )
    if not state.get("files_cleared"):
        for item in draft.get("files", []):
            check(client.delete(item["links"]["self"], timeout=120), f"Delete inherited file {item['filename']}")
        state["files_cleared"] = True
        STATE.write_text(json.dumps(state, indent=2), encoding="utf-8")

    uploaded = set(state.get("uploaded_parts", []))
    for index, path in enumerate(upload_files, start=1):
        if path.name in uploaded:
            print(f"Skipping uploaded file {index}/{len(upload_files)}: {path.name}", flush=True)
            continue
        print(f"Uploading file {index}/{len(upload_files)}: {path.name}", flush=True)
        upload(client, f"{state['bucket_url']}/{path.name}", path)
        uploaded.add(path.name)
        state["uploaded_parts"] = sorted(uploaded)
        STATE.write_text(json.dumps(state, indent=2), encoding="utf-8")

    metadata = {
        "metadata": {
            "title": "Cross-Cohort Transcriptomic Analysis Associates TIMP1 with Tubular Injury-Repair, Extracellular Matrix Remodeling, and a Conserved Disease-State Program in Human Kidney Disease",
            "upload_type": "software",
            "description": (
                "Version 1.1.0 adds signature provenance and overlap audits, bootstrap confidence "
                "intervals, broad injury-factor adjustment, leave-one-discovery-cohort validation, "
                "a matched-random-module benchmark, revised figures, and updated reproducibility files."
            ),
            "creators": [
                {"name": "Ji, Yanzhao", "affiliation": "Shanxi Bethune Hospital, Shanxi Academy of Medical Sciences, Third Hospital of Shanxi Medical University, Tongji Shanxi Hospital"},
                {"name": "Gao, Zhihong", "affiliation": "Shanxi Bethune Hospital, Shanxi Academy of Medical Sciences, Third Hospital of Shanxi Medical University, Tongji Shanxi Hospital"},
            ],
            "keywords": ["TIMP1", "acute kidney injury", "chronic kidney disease", "maladaptive repair", "extracellular matrix", "fibrosis", "transcriptomics"],
            "license": "mit",
            "version": "1.1.0",
            "language": "eng",
            "related_identifiers": [
                {"identifier": "https://github.com/e6denjee-create/TIMP1-AKI-CKD-transcriptomics", "relation": "isSupplementTo", "scheme": "url"},
                {"identifier": "https://github.com/e6denjee-create/TIMP1-AKI-CKD-transcriptomics/releases/tag/v1.1.0", "relation": "isVersionOf", "scheme": "url"},
            ],
        }
    }
    check(
        client.put(f"{API}/deposit/depositions/{state['deposition_id']}", json=metadata, timeout=120),
        "Update v1.1.0 metadata",
    )

    if not state.get("published"):
        published = check(
            client.post(f"{API}/deposit/depositions/{state['deposition_id']}/actions/publish", timeout=180),
            "Publish v1.1.0",
        )
        state.update(
            {
                "published": True,
                "doi": published["doi"],
                "doi_url": published["links"]["doi"],
                "record_url": published["links"]["record_html"],
            }
        )
        STATE.write_text(json.dumps(state, indent=2), encoding="utf-8")
        RESULT.write_text(json.dumps(published, indent=2), encoding="utf-8")

    print(f"DOI={state['doi']}")
    print(f"DOI_URL={state['doi_url']}")
    print(f"RECORD_URL={state['record_url']}")


if __name__ == "__main__":
    main()
