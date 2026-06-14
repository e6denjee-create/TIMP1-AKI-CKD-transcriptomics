"""Create and split the Zenodo v1.1.0 archive."""

from __future__ import annotations

import csv
import hashlib
import shutil
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAYLOAD = ROOT / "public_repository" / "TIMP1_AKI_CKD_transcriptomics_zenodo_payload"
ARCHIVE = ROOT / "public_repository" / "TIMP1_AKI_CKD_transcriptomics_zenodo_v1.1.0.zip"
PARTS_DIR = ROOT / "public_repository" / "zenodo_upload_parts_v1.1.0"
PART_SIZE = 20 * 1024 * 1024


def digest(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def main() -> None:
    if not PAYLOAD.exists():
        raise SystemExit(f"Payload directory not found: {PAYLOAD}")

    ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
    if ARCHIVE.exists():
        ARCHIVE.unlink()
    with zipfile.ZipFile(ARCHIVE, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for path in sorted(PAYLOAD.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(PAYLOAD))
    print(f"Created {ARCHIVE.name} ({ARCHIVE.stat().st_size} bytes)", flush=True)

    if PARTS_DIR.exists():
        shutil.rmtree(PARTS_DIR)
    PARTS_DIR.mkdir(parents=True)

    rows = []
    with ARCHIVE.open("rb") as source:
        index = 1
        while chunk := source.read(PART_SIZE):
            part = PARTS_DIR / f"{ARCHIVE.name}.part{index:03d}"
            part.write_bytes(chunk)
            rows.append(
                {
                    "filename": part.name,
                    "size_bytes": part.stat().st_size,
                    "sha256": digest(part),
                }
            )
            print(f"Created {part.name}", flush=True)
            index += 1

    with (PARTS_DIR / "PARTS_MANIFEST.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["filename", "size_bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)

    (PARTS_DIR / "README_REASSEMBLY.txt").write_text(
        "Reassemble the versioned ZIP by concatenating all .part files in filename order.\n"
        f"Expected archive SHA-256: {digest(ARCHIVE)}\n",
        encoding="utf-8",
    )
    print(f"Prepared {len(rows)} upload parts in {PARTS_DIR}", flush=True)


if __name__ == "__main__":
    main()
