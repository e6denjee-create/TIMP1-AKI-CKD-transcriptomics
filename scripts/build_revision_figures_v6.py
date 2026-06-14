"""Build revised manuscript figures without modifying immutable MVP outputs."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "timp1_validation"
FIGURES = ROOT / "figures" / "timp1_validation"


def save_figure(fig: plt.Figure, stem: str) -> None:
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES / f"{stem}.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_workflow() -> None:
    fig, ax = plt.subplots(figsize=(12, 4.2))
    boxes = [
        (0.2, "#D9E8F5", "Discovery cohorts\n3 bulk datasets"),
        (2.4, "#DDEED9", "External tubular cohort\nGSE180394\n44 disease / 9 control"),
        (
            4.6,
            "#F8E3C5",
            "Robustness analyses\nMatched random programs\nLeave-one-cohort\nLeave-one-signature-out\ninjury PC",
        ),
        (
            6.8,
            "#E8DDF1",
            "Manuscript evidence\nExpression\nPrespecified signatures\nProgram score\nCore genes",
        ),
    ]
    for x, color, label in boxes:
        ax.add_patch(Rectangle((x, 0.25), 1.6, 1.0, facecolor=color, edgecolor="#333333"))
        ax.text(x + 0.8, 0.75, label, ha="center", va="center", fontsize=9.2)
    for start in (1.82, 4.02, 6.22):
        ax.add_patch(
            FancyArrowPatch(
                (start, 0.75),
                (start + 0.5, 0.75),
                arrowstyle="->",
                mutation_scale=18,
                linewidth=1.5,
                color="black",
            )
        )
    ax.text(
        4.3,
        1.58,
        "TIMP1 kidney injury transcriptomic validation framework",
        ha="center",
        va="center",
        fontsize=15,
        fontweight="bold",
    )
    ax.set_xlim(0, 8.6)
    ax.set_ylim(0, 1.8)
    ax.axis("off")
    save_figure(fig, "v6_manuscript_figure_1_workflow_cohort_overview")


def build_program_scatter() -> None:
    scores = pd.read_csv(RESULTS / "external_GSE180394_stringent_module_scores.csv")
    expression = pd.read_csv(
        RESULTS / "external_GSE180394_sensitivity_control_source_data.csv"
    )
    expression = expression[
        expression["control_definition"].eq("living_donor_only")
    ].drop_duplicates("sample")
    data = scores.merge(
        expression[["sample", "TIMP1_expression"]], on="sample", how="inner"
    )
    rho_all = spearmanr(
        data["TIMP1_expression"], data["stringent_TIMP1_module_score"]
    ).statistic
    disease = data[data["group"].eq("Disease")]
    rho_disease = spearmanr(
        disease["TIMP1_expression"], disease["stringent_TIMP1_module_score"]
    ).statistic

    fig, ax = plt.subplots(figsize=(7, 5.8))
    colors = {"Control": "#4C78A8", "Disease": "#C44E52"}
    for group, frame in data.groupby("group"):
        ax.scatter(
            frame["TIMP1_expression"],
            frame["stringent_TIMP1_module_score"],
            s=38,
            alpha=0.8,
            color=colors[group],
            label=group,
        )
    x = data["TIMP1_expression"].to_numpy()
    y = data["stringent_TIMP1_module_score"].to_numpy()
    slope, intercept = np.polyfit(x, y, 1)
    x_line = np.linspace(x.min(), x.max(), 100)
    ax.plot(x_line, intercept + slope * x_line, color="#333333", linewidth=1.6)
    ax.text(
        0.01,
        0.97,
        f"All samples rho = {rho_all:.2f}\nDisease-only rho = {rho_disease:.2f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11,
    )
    ax.set_title(
        "TIMP1 correlates with the discovery-derived program score in GSE180394",
        fontsize=14,
        fontweight="bold",
    )
    ax.set_xlabel("TIMP1 expression")
    ax.set_ylabel("Disease-state program score")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.01), ncol=2)
    ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, "v6_manuscript_figure_4_GSE180394_program_score_scatter")


if __name__ == "__main__":
    FIGURES.mkdir(parents=True, exist_ok=True)
    build_workflow()
    build_program_scatter()
