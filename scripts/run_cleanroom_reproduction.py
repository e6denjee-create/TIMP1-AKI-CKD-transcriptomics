"""Run and verify the core TIMP1 analysis from a clean repository checkout."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "timp1_validation"
FIGURES = ROOT / "figures" / "timp1_validation"


def find_rscript() -> str:
    configured = os.environ.get("RSCRIPT", "").strip()
    candidates = [
        configured,
        shutil.which("Rscript") or "",
        r"C:\Program Files\R\R-4.6.0\bin\Rscript.exe",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return candidate
    raise SystemExit("Rscript was not found. Set the RSCRIPT environment variable.")


def run(command: list[str]) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def read_rows(name: str) -> list[dict[str, str]]:
    path = RESULTS / name
    if not path.exists():
        raise AssertionError(f"Expected output is missing: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def near(actual: str, expected: float, tolerance: float, label: str) -> None:
    value = float(actual)
    if abs(value - expected) > tolerance:
        raise AssertionError(
            f"{label}: expected {expected} +/- {tolerance}, observed {value}"
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_inputs() -> None:
    required = [
        *(ROOT / "data/processed" / f"{dataset}_normalized_expression.csv.gz"
          for dataset in ("GSE139061", "GSE30718", "GSE66494")),
        *(ROOT / "data/metadata" / f"{dataset}_metadata.csv"
          for dataset in ("GSE139061", "GSE30718", "GSE66494")),
        ROOT / "data/external/GSE180394/external_GSE180394_expression_gene_symbol.csv.gz",
        ROOT / "data/external/GSE180394/external_GSE180394_expression_all_samples_gene_symbol.csv.gz",
        ROOT / "data/external/GSE180394/external_GSE180394_metadata.csv",
        ROOT / "data/external/GSE180394/external_GSE180394_metadata_all_samples.csv",
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    if missing:
        raise AssertionError("Required clean-room inputs are missing:\n- " + "\n- ".join(missing))

    manifest_path = ROOT / "INPUT_MANIFEST.csv"
    if not manifest_path.exists():
        raise AssertionError("Required input manifest is missing: INPUT_MANIFEST.csv")
    with manifest_path.open(encoding="utf-8-sig", newline="") as handle:
        manifest = list(csv.DictReader(handle))
    if not manifest:
        raise AssertionError("INPUT_MANIFEST.csv is empty")
    for row in manifest:
        relative_path = row["relative_path"]
        input_path = ROOT / relative_path
        if not input_path.exists():
            raise AssertionError(f"Manifest input is missing: {relative_path}")
        observed_size = input_path.stat().st_size
        if observed_size != int(row["size_bytes"]):
            raise AssertionError(
                f"Input size mismatch for {relative_path}: "
                f"expected {row['size_bytes']}, observed {observed_size}"
            )
        observed_hash = sha256(input_path)
        if observed_hash != row["sha256"]:
            raise AssertionError(
                f"Input hash mismatch for {relative_path}: "
                f"expected {row['sha256']}, observed {observed_hash}"
            )


def reset_outputs() -> None:
    for directory in (RESULTS, FIGURES):
        if directory.exists():
            shutil.rmtree(directory)
        directory.mkdir(parents=True)


def validate_outputs() -> None:
    bulk = {row["dataset"]: row for row in read_rows("TIMP1_bulk_validation_statistics.csv")}
    expected_counts = {
        "GSE139061": (39, 9),
        "GSE30718": (28, 11),
        "GSE66494": (53, 8),
    }
    for dataset, (disease_n, control_n) in expected_counts.items():
        row = bulk[dataset]
        assert int(row["n_disease"]) == disease_n
        assert int(row["n_control"]) == control_n

    stringent = read_rows("stringent_TIMP1_correlated_module.csv")
    assert len(stringent) == 421, f"Expected 421 stringent genes, observed {len(stringent)}"

    external = read_rows("external_GSE180394_TIMP1_group_statistics.csv")[0]
    assert int(external["n_disease"]) == 44
    assert int(external["n_control"]) == 9
    near(external["hedges_g"], 1.5688851215, 1e-8, "External Hedges' g")
    near(external["p_value"], 8.006238931e-05, 1e-10, "External Wilcoxon P")

    benchmark = read_rows("module_random_matched_benchmark_v6.csv")[0]
    near(benchmark["observed_rho"], 0.8496124031, 1e-8, "Program-score rho")
    near(benchmark["empirical_one_sided_p"], 0.000999001, 1e-9, "Random benchmark P")

    unified = read_rows("TIMP1_expression_statistics_unified_bootstrap_v7.csv")
    primary = next(
        row for row in unified
        if row["dataset"] == "GSE180394"
        and row["contrast"] == "Disease_vs_Living_Donor"
    )
    near(primary["hedges_g_ci_lower"], 1.1104851622, 1e-8, "Primary bootstrap lower CI")
    near(primary["hedges_g_ci_upper"], 2.2403888384, 1e-8, "Primary bootstrap upper CI")

    sensitivity = {
        row["contrast"]: row
        for row in unified
        if row["dataset"] == "GSE180394"
        and row["contrast"] != "Disease_vs_Living_Donor"
    }
    legacy_sensitivity = {
        row["control_definition"]: row
        for row in read_rows("external_GSE180394_sensitivity_control_statistics.csv")
    }
    contrast_map = {
        "Disease_vs_Living_Donor_plus_Tumor_Nephrectomy": "extended_controls",
        "Disease_vs_Tumor_Nephrectomy": "tumor_nephrectomy_only",
    }
    for contrast, definition in contrast_map.items():
        authoritative = sensitivity[contrast]
        legacy = legacy_sensitivity[definition]
        for column in (
            "hedges_g",
            "hedges_g_ci_lower",
            "hedges_g_ci_upper",
            "p_value",
            "bh_adjusted_p_value",
        ):
            near(
                legacy[column],
                float(authoritative[column]),
                1e-12,
                f"{definition} unified {column}",
            )

    required_figures = [
        FIGURES / "v6_manuscript_figure_1_workflow_cohort_overview.png",
        FIGURES / "v6_manuscript_figure_4_GSE180394_program_score_scatter.png",
        FIGURES / "v6_module_random_matched_benchmark.png",
    ]
    missing_figures = [str(path.relative_to(ROOT)) for path in required_figures if not path.exists()]
    if missing_figures:
        raise AssertionError("Required figures are missing:\n- " + "\n- ".join(missing_figures))


def write_report() -> None:
    report = {
        "status": "CLEANROOM_REPRODUCTION_OK",
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "repository_root": ".",
        "input_manifest": "INPUT_MANIFEST.csv",
        "checks": [
            "required inputs present and SHA256 hashes matched",
            "bulk sample counts matched the frozen analysis",
            "stringent TIMP1 module contained 421 genes",
            "external GSE180394 effect size and P value matched",
            "program-score correlation and empirical P value matched",
            "unified bootstrap interval matched",
            "required publication figures were generated",
        ],
    }
    report_path = RESULTS / "cleanroom_reproduction_report.json"
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--keep-existing",
        action="store_true",
        help="Do not delete existing results and figures before reproduction.",
    )
    args = parser.parse_args()
    validate_inputs()
    if not args.keep_existing:
        reset_outputs()

    rscript = find_rscript()
    workflows = [
        [rscript, "--vanilla", "scripts/run_timp1_bulk_validation.R"],
        [rscript, "--vanilla", "scripts/run_timp1_module_network.R"],
        [rscript, "--vanilla", "scripts/run_timp1_module_sensitivity_v2.R"],
        [
            rscript,
            "--vanilla",
            "scripts/run_external_bulk_validation_template.R",
            "--dataset",
            "GSE180394",
            "--expression",
            "data/external/GSE180394/external_GSE180394_expression_gene_symbol.csv.gz",
            "--metadata",
            "data/external/GSE180394/external_GSE180394_metadata.csv",
            "--gene-id-type",
            "symbol",
        ],
        [rscript, "--vanilla", "scripts/run_external_GSE180394_sensitivity_v4.R"],
        [rscript, "--vanilla", "scripts/run_revision_robustness_v5.R"],
        [rscript, "--vanilla", "scripts/run_revision_robustness_v6.R"],
        [rscript, "--vanilla", "scripts/run_unified_bootstrap_v7.R"],
        [sys.executable, "scripts/build_revision_figures_v6.py"],
    ]
    for command in workflows:
        run(command)
    validate_outputs()
    write_report()
    print("CLEANROOM_REPRODUCTION_OK", flush=True)


if __name__ == "__main__":
    main()
