"""Split the Zenodo archive into proxy-friendly parts with checksums."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "public_repository"
    / "TIMP1_AKI_CKD_transcriptomics_zenodo_v1.0.0.zip"
)
PARTS_DIR = ROOT / "public_repository" / "zenodo_upload_parts"
PART_SIZE = 20 * 1024 * 1024


def digest(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def main() -> None:
    PARTS_DIR.mkdir(parents=True, exist_ok=True)
    for old_part in PARTS_DIR.glob(f"{SOURCE.name}.part*"):
        old_part.unlink()

    rows = []
    with SOURCE.open("rb") as source:
        index = 1
        while True:
            chunk = source.read(PART_SIZE)
            if not chunk:
                break
            part = PARTS_DIR / f"{SOURCE.name}.part{index:03d}"
            part.write_bytes(chunk)
            rows.append(
                {
                    "filename": part.name,
                    "size_bytes": part.stat().st_size,
                    "sha256": digest(part),
                }
            )
            print(f"Created {part.name} ({part.stat().st_size} bytes)")
            index += 1

    manifest = PARTS_DIR / "PARTS_MANIFEST.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["filename", "size_bytes", "sha256"]
        )
        writer.writeheader()
        writer.writerows(rows)

    instructions = PARTS_DIR / "README_REASSEMBLY.txt"
    instructions.write_text(
        "Reassemble on Windows PowerShell:\n"
        "$parts = Get-ChildItem "
        f"'{SOURCE.name}.part*' | Sort-Object Name\n"
        f"$out = [System.IO.File]::Create('{SOURCE.name}')\n"
        "try { foreach ($part in $parts) { "
        "$bytes = [System.IO.File]::ReadAllBytes($part.FullName); "
        "$out.Write($bytes, 0, $bytes.Length) } } finally { $out.Dispose() }\n\n"
        "Reassemble on Linux/macOS:\n"
        f"cat {SOURCE.name}.part* > {SOURCE.name}\n\n"
        f"Expected complete archive SHA-256: {digest(SOURCE)}\n",
        encoding="utf-8",
    )
    print(f"Wrote {len(rows)} parts to {PARTS_DIR}")


if __name__ == "__main__":
    main()
