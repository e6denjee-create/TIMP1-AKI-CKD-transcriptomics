"""Publish the prepared TIMP1 research compendium to Zenodo.

The access token is read from ZENODO_TOKEN and is never written to disk.
Progress state excludes credentials and allows safe retries before publication.
"""

from __future__ import annotations

import hashlib
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
ARCHIVE = (
    ROOT
    / "public_repository"
    / "TIMP1_AKI_CKD_transcriptomics_zenodo_v1.0.0.zip"
)
PARTS_DIR = ROOT / "public_repository" / "zenodo_upload_parts"
STATE = ROOT / "public_repository" / "zenodo_deposition_state.json"
RESULT = ROOT / "public_repository" / "zenodo_published_record.json"
API = "https://zenodo.org/api"


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
                allowed_methods=("GET", "POST", "PUT"),
            )
        ),
    )
    return client


def check(response: requests.Response, action: str) -> dict:
    if not response.ok:
        raise RuntimeError(
            f"{action} failed: HTTP {response.status_code}\n{response.text[:2000]}"
        )
    if not response.content:
        return {}
    return response.json()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def upload_archive(
    client: requests.Session, upload_url: str, archive: Path
) -> dict:
    """Upload a large file with explicit rewind-and-retry semantics."""
    last_error: Exception | None = None
    for attempt in range(1, 9):
        try:
            with archive.open("rb") as handle:
                response = client.put(
                    upload_url,
                    data=handle,
                    headers={
                        "Content-Type": "application/octet-stream",
                        "Content-Length": str(archive.stat().st_size),
                    },
                    timeout=(90, 7200),
                )
            return check(response, f"Upload archive (attempt {attempt})")
        except RequestException as error:
            last_error = error
            wait_seconds = min(15 * attempt, 90)
            print(
                f"Upload attempt {attempt}/8 was interrupted: "
                f"{type(error).__name__}. Retrying in {wait_seconds} seconds.",
                flush=True,
            )
            time.sleep(wait_seconds)
    raise RuntimeError(f"Archive upload failed after 8 attempts: {last_error}")


def main() -> None:
    token = os.environ.get("ZENODO_TOKEN", "").strip()
    if not token:
        sys.exit("ZENODO_TOKEN is not set.")
    if not ARCHIVE.exists():
        sys.exit(f"Archive not found: {ARCHIVE}")
    upload_files = sorted(PARTS_DIR.glob(f"{ARCHIVE.name}.part*"))
    upload_files.extend(
        [
            PARTS_DIR / "PARTS_MANIFEST.csv",
            PARTS_DIR / "README_REASSEMBLY.txt",
        ]
    )
    upload_files = [path for path in upload_files if path.exists()]
    if len(upload_files) < 3:
        sys.exit(
            "Split upload parts are missing. Run scripts/split_zenodo_archive.py."
        )

    client = session(token)
    state = json.loads(STATE.read_text(encoding="utf-8")) if STATE.exists() else {}

    if state.get("deposition_id"):
        access_check = client.get(
            f"{API}/deposit/depositions/{state['deposition_id']}",
            timeout=120,
        )
        if access_check.status_code == 403:
            sys.exit(
                "TOKEN_ACCESS_DENIED: The current Zenodo token cannot edit "
                f"deposition {state['deposition_id']}. Use a token created by "
                "the same Zenodo account that created this draft, with "
                "deposit:write and deposit:actions permissions."
            )
        check(access_check, "Validate deposition access")

    if not state.get("deposition_id"):
        deposition = check(
            client.post(f"{API}/deposit/depositions", json={}, timeout=120),
            "Create deposition",
        )
        state = {
            "deposition_id": deposition["id"],
            "bucket_url": deposition["links"]["bucket"],
            "doi": deposition.get("metadata", {}).get("prereserve_doi", {}).get(
                "doi", ""
            ),
            "archive_sha256": sha256(ARCHIVE),
            "uploaded_parts": [],
            "uploaded": False,
            "published": False,
        }
        STATE.write_text(json.dumps(state, indent=2), encoding="utf-8")
        print(
            f"Created deposition {state['deposition_id']}; reserved DOI "
            f"{state.get('doi') or '[pending]'}"
        )

    if not state.get("uploaded"):
        uploaded_parts = set(state.get("uploaded_parts", []))
        for position, upload_file in enumerate(upload_files, start=1):
            if upload_file.name in uploaded_parts:
                print(
                    f"Skipping uploaded file {position}/{len(upload_files)}: "
                    f"{upload_file.name}"
                )
                continue
            print(
                f"Uploading file {position}/{len(upload_files)}: "
                f"{upload_file.name} ({upload_file.stat().st_size} bytes)",
                flush=True,
            )
            upload_url = f"{state['bucket_url']}/{upload_file.name}"
            upload_archive(client, upload_url, upload_file)
            uploaded_parts.add(upload_file.name)
            state["uploaded_parts"] = sorted(uploaded_parts)
            STATE.write_text(json.dumps(state, indent=2), encoding="utf-8")
            print(f"Uploaded {upload_file.name}", flush=True)
        state["uploaded"] = True
        STATE.write_text(json.dumps(state, indent=2), encoding="utf-8")
        print("Uploaded all archive parts")

    metadata = {
        "metadata": {
            "title": (
                "Cross-Cohort Transcriptomic Analysis Associates TIMP1 with "
                "Tubular Injury-Repair, Extracellular Matrix Remodeling, and "
                "a Conserved Disease-State Program in Human Kidney Disease"
            ),
            "upload_type": "software",
            "description": (
                "Reproducible research compendium containing analysis scripts, "
                "processed and intermediate data, statistical outputs, session "
                "information, source data, and publication figures for an "
                "integrative transcriptomic analysis of TIMP1 in human kidney "
                "injury and fibrosis. Results are association-based and do not "
                "establish diagnostic utility, kidney specificity, longitudinal "
                "prediction, or causality."
            ),
            "creators": [
                {
                    "name": "Ji, Yanzhao",
                    "affiliation": (
                        "Shanxi Bethune Hospital, Shanxi Academy of Medical "
                        "Sciences, Third Hospital of Shanxi Medical University, "
                        "Tongji Shanxi Hospital"
                    ),
                },
                {
                    "name": "Gao, Zhihong",
                    "affiliation": (
                        "Shanxi Bethune Hospital, Shanxi Academy of Medical "
                        "Sciences, Third Hospital of Shanxi Medical University, "
                        "Tongji Shanxi Hospital"
                    ),
                },
            ],
            "keywords": [
                "TIMP1",
                "acute kidney injury",
                "chronic kidney disease",
                "maladaptive repair",
                "extracellular matrix",
                "fibrosis",
                "transcriptomics",
            ],
            "license": "mit",
            "version": "1.0.0",
            "language": "eng",
            "related_identifiers": [
                {
                    "identifier": (
                        "https://github.com/e6denjee-create/"
                        "TIMP1-AKI-CKD-transcriptomics"
                    ),
                    "relation": "isSupplementTo",
                    "scheme": "url",
                },
                {
                    "identifier": (
                        "https://github.com/e6denjee-create/"
                        "TIMP1-AKI-CKD-transcriptomics/releases/tag/v1.0.0"
                    ),
                    "relation": "isVersionOf",
                    "scheme": "url",
                },
            ],
        }
    }
    check(
        client.put(
            f"{API}/deposit/depositions/{state['deposition_id']}",
            json=metadata,
            timeout=120,
        ),
        "Update metadata",
    )
    print("Updated deposition metadata")

    if not state.get("published"):
        published = check(
            client.post(
                f"{API}/deposit/depositions/{state['deposition_id']}/actions/publish",
                timeout=180,
            ),
            "Publish deposition",
        )
        state["published"] = True
        state["doi"] = published["doi"]
        state["doi_url"] = published["links"]["doi"]
        state["record_url"] = published["links"]["record_html"]
        STATE.write_text(json.dumps(state, indent=2), encoding="utf-8")
        RESULT.write_text(json.dumps(published, indent=2), encoding="utf-8")

    print(f"DOI={state['doi']}")
    print(f"DOI_URL={state['doi_url']}")
    print(f"RECORD_URL={state['record_url']}")


if __name__ == "__main__":
    main()
