import gzip
import io
import math
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy import stats
from sklearn.metrics import auc, roc_curve

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

GSES = ["GSE139061", "GSE30718", "GSE66494"]
TARGET = "LOXL4"
FIBROSIS_GENES = ["COL1A1", "COL3A1", "FN1", "TGFB1"]

PALETTES = {
    "Control": "#4C78A8",
    "AKI": "#D65F5F",
    "CKD": "#B279A2",
}

GENE_SETS = {
    "Fibrosis core": [
        "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL5A1", "FN1", "ACTA2",
        "TGFB1", "TGFBR1", "TGFBR2", "CTGF", "POSTN", "SPARC", "MMP2",
        "TIMP1", "LOX", "LOXL1", "LOXL2", "LOXL4", "VIM", "SNAI1", "SNAI2",
    ],
    "ECM organization": [
        "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL6A1",
        "COL6A2", "COL6A3", "FN1", "LAMB1", "LAMC1", "ITGA5", "ITGB1",
        "MMP2", "MMP9", "TIMP1", "SPARC", "POSTN", "THBS1", "VCAN",
    ],
    "TGF beta signaling": [
        "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3",
        "SMAD4", "SMAD7", "CTGF", "SERPINE1", "ID1", "ID2", "JUN", "FOS",
        "SNAI1", "SNAI2", "ZEB1", "ZEB2",
    ],
    "Kidney injury repair": [
        "HAVCR1", "LCN2", "VCAM1", "PROM1", "SOX9", "VIM", "KRT8", "KRT18",
        "KRT19", "MMP7", "TIMP1", "IL6", "CCL2", "CXCL1", "CXCL8", "STAT3",
        "JUN", "FOS", "EGR1",
    ],
    "Inflammation": [
        "IL1B", "IL6", "TNF", "NFKB1", "NFKBIA", "CCL2", "CCL5", "CXCL8",
        "CXCL10", "ICAM1", "VCAM1", "TLR2", "TLR4", "STAT1", "STAT3", "IRF1",
        "HLA-DRA", "CD68", "LST1",
    ],
    "Hypoxia oxidative stress": [
        "HIF1A", "VEGFA", "SLC2A1", "LDHA", "ENO1", "NQO1", "HMOX1", "SOD2",
        "CAT", "GPX1", "TXN", "TXNRD1", "DDIT4", "BNIP3", "P4HA1", "PLOD2",
    ],
}


def mkdirs():
    for d in (DATA, RESULTS, FIGURES):
        d.mkdir(parents=True, exist_ok=True)


def download(url, path):
    if path.exists() and path.stat().st_size > 100:
        return
    response = requests.get(url, timeout=240, stream=True)
    response.raise_for_status()
    with open(path, "wb") as handle:
        for chunk in response.iter_content(1024 * 256):
            if chunk:
                handle.write(chunk)


def download_geo_files():
    for gse in GSES:
        prefix = gse[:-3] + "nnn"
        download(
            f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/matrix/{gse}_series_matrix.txt.gz",
            DATA / f"{gse}_series_matrix.txt.gz",
        )
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE139nnn/GSE139061/suppl/GSE139061_Eadon_processed_QN_101419.csv.gz",
        DATA / "GSE139061_Eadon_processed_QN_101419.csv.gz",
    )
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
        DATA / "GPL570.annot.gz",
    )
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6480/annot/GPL6480.annot.gz",
        DATA / "GPL6480.annot.gz",
    )


def split_geo_line(line):
    return [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]


def parse_series_matrix(gse):
    path = DATA / f"{gse}_series_matrix.txt.gz"
    meta = {}
    table_lines = []
    in_table = False
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            if in_table:
                table_lines.append(line)
            elif line.startswith("!Sample_"):
                key, *vals = split_geo_line(line)
                meta.setdefault(key, []).append(vals)
            elif line.startswith("!Series_platform_id"):
                meta["platform"] = split_geo_line(line)[1:]

    expr = pd.read_csv(io.StringIO("".join(table_lines)), sep="\t", index_col=0)
    expr = expr.apply(pd.to_numeric, errors="coerce")

    samples = meta["!Sample_geo_accession"][0]
    sample_meta = pd.DataFrame({"sample": samples})
    titles = meta.get("!Sample_title", [[""] * len(samples)])[0]
    sources = meta.get("!Sample_source_name_ch1", [[""] * len(samples)])[0]
    sample_meta["title"] = titles
    sample_meta["source"] = sources
    for i, vals in enumerate(meta.get("!Sample_characteristics_ch1", []), start=1):
        sample_meta[f"characteristics_{i}"] = vals
    return expr, sample_meta


def parse_platform_symbols(gpl):
    path = DATA / f"{gpl}.annot.gz"
    rows = []
    in_table = False
    with gzip.open(path, "rt", errors="replace") as handle:
        header = None
        for line in handle:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table and header is None:
                header = line.rstrip("\n").split("\t")
                continue
            if in_table:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= len(header):
                    rows.append(parts[: len(header)])
    anno = pd.DataFrame(rows, columns=header)
    anno = anno[["ID", "Gene symbol"]].rename(columns={"Gene symbol": "gene_symbol"})
    anno["gene_symbol"] = anno["gene_symbol"].astype(str).str.split("///").str[0].str.strip()
    anno = anno[(anno["gene_symbol"] != "") & (anno["gene_symbol"] != "nan")]
    return anno.drop_duplicates("ID")


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    tmp = np.empty(n)
    tmp[order] = adj
    out[valid] = tmp
    return out


def quantile_normalize(df):
    values = df.to_numpy(dtype=float)
    sorted_values = np.sort(values, axis=0)
    means = np.nanmean(sorted_values, axis=1)
    ranks = np.apply_along_axis(lambda x: stats.rankdata(x, method="min") - 1, 0, values).astype(int)
    normalized = np.zeros_like(values, dtype=float)
    for col in range(values.shape[1]):
        normalized[:, col] = means[ranks[:, col]]
    return pd.DataFrame(normalized, index=df.index, columns=df.columns)


def collapse_to_gene(expr, anno=None):
    work = expr.copy()
    if anno is not None:
        work = work.reset_index().rename(columns={work.index.name or "index": "ID"})
        work = work.merge(anno, on="ID", how="inner")
        work = work.drop(columns=["ID"]).rename(columns={"gene_symbol": "Gene"})
    else:
        work = work.reset_index().rename(columns={work.index.name or "index": "Gene"})
    work["Gene"] = work["Gene"].astype(str).str.split("///").str[0].str.strip()
    sample_cols = [c for c in work.columns if c != "Gene"]
    work[sample_cols] = work[sample_cols].apply(pd.to_numeric, errors="coerce")
    return work.groupby("Gene", as_index=True)[sample_cols].mean()


def prepare_dataset(gse):
    if gse == "GSE139061":
        expr = pd.read_csv(DATA / "GSE139061_Eadon_processed_QN_101419.csv.gz", index_col=0)
        expr = expr.apply(pd.to_numeric, errors="coerce")
        expr = np.log2(expr + 1)
        expr = quantile_normalize(expr)
        meta = pd.DataFrame({"sample": expr.columns})
        meta["group"] = np.where(meta["sample"].str.startswith("AKI"), "AKI", "Control")
        meta["comparison"] = "AKI_vs_Control"
        gene_expr = collapse_to_gene(expr, None)
    elif gse == "GSE30718":
        expr, meta = parse_series_matrix(gse)
        anno = parse_platform_symbols("GPL570")
        char_text = meta.filter(like="characteristics").fillna("").agg(" ".join, axis=1).str.lower()
        source = meta["source"].str.lower()
        meta["group"] = np.where(char_text.str.contains("acute kidney injury"), "AKI",
                                 np.where(source.str.contains("protocol"), "Control", "Exclude"))
        meta = meta[meta["group"] != "Exclude"].copy()
        expr = expr.loc[:, meta["sample"]]
        if np.nanmax(expr.to_numpy()) > 100:
            expr = np.log2(expr + 1)
        expr = quantile_normalize(expr)
        meta["comparison"] = "AKI_vs_Control"
        gene_expr = collapse_to_gene(expr, anno)
    elif gse == "GSE66494":
        expr, meta = parse_series_matrix(gse)
        anno = parse_platform_symbols("GPL6480")
        char_text = meta.filter(like="characteristics").fillna("").agg(" ".join, axis=1).str.lower()
        meta["group"] = np.where(char_text.str.contains("chronic kidney disease"), "CKD",
                                 np.where(char_text.str.contains("normal kidney"), "Control", "Exclude"))
        meta = meta[meta["group"] != "Exclude"].copy()
        expr = expr.loc[:, meta["sample"]]
        if np.nanmax(expr.to_numpy()) > 100:
            expr = np.log2(expr + 1)
        expr = quantile_normalize(expr)
        meta["comparison"] = "CKD_vs_Control"
        gene_expr = collapse_to_gene(expr, anno)
    else:
        raise ValueError(gse)

    gene_expr = gene_expr.loc[gene_expr.notna().sum(axis=1) >= max(3, gene_expr.shape[1] // 4)]
    gene_expr.to_csv(RESULTS / f"{gse}_normalized_gene_expression.csv")
    meta.to_csv(RESULTS / f"{gse}_sample_metadata.csv", index=False)
    return gene_expr, meta


def differential_expression(gse, expr, meta):
    samples_case = meta.loc[meta["group"].isin(["AKI", "CKD"]), "sample"].tolist()
    samples_ctrl = meta.loc[meta["group"] == "Control", "sample"].tolist()
    case = expr[samples_case]
    ctrl = expr[samples_ctrl]
    logfc = case.mean(axis=1) - ctrl.mean(axis=1)
    stat, pval = stats.ttest_ind(case.T, ctrl.T, equal_var=False, nan_policy="omit")
    res = pd.DataFrame({
        "Gene": expr.index,
        "logFC": logfc.values,
        "P.Value": pval,
    })
    res["adj.P.Val"] = bh_adjust(res["P.Value"])
    res["abs_logFC"] = res["logFC"].abs()
    res = res.sort_values(["P.Value", "abs_logFC"], ascending=[True, False]).drop(columns="abs_logFC")
    res.to_csv(RESULTS / f"{gse}_differential_expression.csv", index=False)
    res.loc[res["Gene"].eq(TARGET)].to_csv(RESULTS / f"{gse}_LOXL4_differential_expression.csv", index=False)
    return res


def roc_analysis(gse, expr, meta):
    target_expr = expr.loc[TARGET, meta["sample"]].astype(float)
    y = meta["group"].isin(["AKI", "CKD"]).astype(int).to_numpy()
    scores = target_expr.to_numpy()
    fpr, tpr, thresholds = roc_curve(y, scores)
    roc_auc = auc(fpr, tpr)
    roc_df = pd.DataFrame({"FPR": fpr, "TPR": tpr, "threshold": thresholds})
    roc_df["AUC"] = roc_auc
    roc_df.to_csv(RESULTS / f"{gse}_LOXL4_ROC.csv", index=False)
    return roc_df, roc_auc


def correlation_analysis(gse, expr, meta):
    rows = []
    for gene in FIBROSIS_GENES:
        if gene not in expr.index or TARGET not in expr.index:
            rows.append({"Dataset": gse, "Gene": gene, "rho": np.nan, "P.Value": np.nan, "n": 0})
            continue
        sub = expr.loc[[TARGET, gene], meta["sample"]].T.dropna()
        rho, pval = stats.spearmanr(sub[TARGET], sub[gene])
        rows.append({"Dataset": gse, "Gene": gene, "rho": rho, "P.Value": pval, "n": len(sub)})
    out = pd.DataFrame(rows)
    out["adj.P.Val"] = bh_adjust(out["P.Value"])
    out.to_csv(RESULTS / f"{gse}_LOXL4_fibrosis_correlations.csv", index=False)
    return out


def enrichment_score(ranked_genes, gene_set):
    genes = np.array(ranked_genes["Gene"])
    scores = np.abs(ranked_genes["score"].to_numpy(dtype=float))
    hits = np.isin(genes, list(gene_set))
    nh = hits.sum()
    if nh == 0 or nh == len(hits):
        return np.nan
    hit_weights = np.where(hits, scores / scores[hits].sum(), 0)
    miss_weights = np.where(~hits, 1 / (len(hits) - nh), 0)
    running = np.cumsum(hit_weights - miss_weights)
    max_es = running.max()
    min_es = running.min()
    return max_es if abs(max_es) >= abs(min_es) else min_es


def preranked_gsea(gse, de):
    ranked = de[["Gene", "logFC", "P.Value"]].dropna().copy()
    ranked["score"] = np.sign(ranked["logFC"]) * -np.log10(ranked["P.Value"].clip(lower=1e-300))
    ranked = ranked.sort_values("score", ascending=False).drop_duplicates("Gene")
    rng = np.random.default_rng(20260607)
    rows = []
    for name, genes in GENE_SETS.items():
        genes = set(genes)
        overlap = len(genes.intersection(ranked["Gene"]))
        es = enrichment_score(ranked[["Gene", "score"]], genes)
        null = []
        for _ in range(1000):
            sample = set(rng.choice(ranked["Gene"].to_numpy(), size=max(overlap, 1), replace=False))
            null.append(enrichment_score(ranked[["Gene", "score"]], sample))
        null = np.array([x for x in null if np.isfinite(x)])
        if not np.isfinite(es) or len(null) == 0:
            nes, pval = np.nan, np.nan
        elif es >= 0:
            pos = null[null >= 0]
            nes = es / np.mean(pos) if len(pos) else np.nan
            pval = (np.sum(pos >= es) + 1) / (len(pos) + 1) if len(pos) else np.nan
        else:
            neg = null[null < 0]
            nes = es / abs(np.mean(neg)) if len(neg) else np.nan
            pval = (np.sum(neg <= es) + 1) / (len(neg) + 1) if len(neg) else np.nan
        rows.append({"Dataset": gse, "Pathway": name, "ES": es, "NES": nes, "P.Value": pval, "Overlap": overlap})
    out = pd.DataFrame(rows)
    out["adj.P.Val"] = bh_adjust(out["P.Value"])
    out.to_csv(RESULTS / f"{gse}_GSEA_preranked.csv", index=False)
    return out


def save_figure(fig, name):
    fig.savefig(FIGURES / f"{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_outputs(all_expr, all_meta, all_de, all_rocs, all_corr, all_gsea):
    sns.set_theme(style="whitegrid", font="Arial")

    # Figure 2: LOXL4 expression and differential effect.
    expr_rows = []
    for gse, expr in all_expr.items():
        meta = all_meta[gse]
        for _, row in meta.iterrows():
            expr_rows.append({"Dataset": gse, "Group": row["group"], "LOXL4": expr.loc[TARGET, row["sample"]]})
    fig2_df = pd.DataFrame(expr_rows)
    fig2_df.to_csv(RESULTS / "Figure2_LOXL4_expression.csv", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(9.4, 4.6), sharey=False)
    for ax, gse in zip(axes, GSES):
        sub = fig2_df[fig2_df["Dataset"] == gse]
        order = ["Control", "AKI"] if "AKI" in set(sub["Group"]) else ["Control", "CKD"]
        sns.boxplot(data=sub, x="Group", y="LOXL4", hue="Group", order=order, hue_order=order, palette=PALETTES,
                    width=0.58, fliersize=0, legend=False, ax=ax)
        sns.stripplot(data=sub, x="Group", y="LOXL4", order=order, color="#222222", alpha=0.52, size=3.1, jitter=0.13, ax=ax)
        ax.set_title(gse)
        ax.set_xlabel("")
        ax.set_ylabel("LOXL4 expression (normalized log2)" if gse == GSES[0] else "")
    fig.suptitle("Figure 2. LOXL4 expression in AKI and CKD datasets", y=1.02)
    save_figure(fig, "Figure2_LOXL4_differential_expression")

    # Figure 3: ROC curves.
    fig, ax = plt.subplots(figsize=(5.7, 5.2))
    summary_rows = []
    colors = ["#D65F5F", "#7A5195", "#4C78A8"]
    for color, (gse, roc_df) in zip(colors, all_rocs.items()):
        roc_auc = roc_df["AUC"].iloc[0]
        ax.plot(roc_df["FPR"], roc_df["TPR"], lw=2.2, color=color, label=f"{gse} AUC={roc_auc:.3f}")
        summary_rows.append({"Dataset": gse, "AUC": roc_auc})
    ax.plot([0, 1], [0, 1], ls="--", lw=1.2, color="#888888")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("Figure 3. Diagnostic ROC for LOXL4")
    ax.legend(frameon=False, loc="lower right")
    pd.DataFrame(summary_rows).to_csv(RESULTS / "Figure3_LOXL4_ROC_summary.csv", index=False)
    save_figure(fig, "Figure3_LOXL4_ROC")

    # Figure 4: correlations.
    corr_df = pd.concat(all_corr.values(), ignore_index=True)
    corr_df.to_csv(RESULTS / "Figure4_LOXL4_fibrosis_correlations.csv", index=False)
    fig, ax = plt.subplots(figsize=(7.6, 4.5))
    pivot = corr_df.pivot(index="Dataset", columns="Gene", values="rho").loc[GSES, FIBROSIS_GENES]
    sns.heatmap(pivot, cmap="vlag", center=0, vmin=-1, vmax=1, annot=True, fmt=".2f",
                linewidths=0.8, linecolor="white", cbar_kws={"label": "Spearman rho"}, ax=ax)
    ax.set_title("Figure 4. LOXL4 correlation with fibrosis genes")
    ax.set_xlabel("")
    ax.set_ylabel("")
    save_figure(fig, "Figure4_LOXL4_fibrosis_correlation_heatmap")

    # Figure 5: GSEA dot plot.
    gsea_df = pd.concat(all_gsea.values(), ignore_index=True)
    gsea_df.to_csv(RESULTS / "Figure5_GSEA_summary.csv", index=False)
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    plot_df = gsea_df.copy()
    plot_df["neg_log10_FDR"] = -np.log10(plot_df["adj.P.Val"].clip(lower=1e-300))
    sns.scatterplot(data=plot_df, x="Dataset", y="Pathway", hue="NES", size="neg_log10_FDR",
                    palette="coolwarm", hue_norm=(-2.5, 2.5), sizes=(55, 260), edgecolor="#222222", linewidth=0.4, ax=ax)
    ax.axhline(-0.5, color="none")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Figure 5. Preranked GSEA across LOXL4-associated disease contrasts")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    save_figure(fig, "Figure5_GSEA_dotplot")

    # Figure 6: integrated evidence.
    rows = []
    auc_map = {gse: all_rocs[gse]["AUC"].iloc[0] for gse in all_rocs}
    mean_corr = corr_df.groupby("Dataset")["rho"].mean().to_dict()
    for gse, de in all_de.items():
        lox = de.loc[de["Gene"].eq(TARGET)].iloc[0]
        rows.append({
            "Dataset": gse,
            "LOXL4_logFC": lox["logFC"],
            "LOXL4_P.Value": lox["P.Value"],
            "LOXL4_adj.P.Val": lox["adj.P.Val"],
            "LOXL4_AUC": auc_map[gse],
            "Mean_fibrosis_rho": mean_corr[gse],
        })
    integrated = pd.DataFrame(rows)
    integrated.to_csv(RESULTS / "Figure6_integrated_LOXL4_evidence.csv", index=False)
    long = integrated.melt(id_vars="Dataset", value_vars=["LOXL4_logFC", "LOXL4_AUC", "Mean_fibrosis_rho"],
                           var_name="Metric", value_name="Value")
    metric_labels = {"LOXL4_logFC": "logFC", "LOXL4_AUC": "AUC", "Mean_fibrosis_rho": "mean rho"}
    long["Metric"] = long["Metric"].map(metric_labels)
    fig, axes = plt.subplots(1, 3, figsize=(9, 4.3), sharey=True)
    for ax, metric in zip(axes, ["logFC", "AUC", "mean rho"]):
        sub = long[long["Metric"] == metric]
        sns.barplot(data=sub, y="Dataset", x="Value", color="#5B8C85", ax=ax)
        ax.axvline(0 if metric != "AUC" else 0.5, color="#555555", lw=1, ls="--")
        ax.set_title(metric)
        ax.set_xlabel("")
        ax.set_ylabel("")
    fig.suptitle("Figure 6. Integrated evidence for LOXL4 in AKI-to-CKD transition", y=1.02)
    save_figure(fig, "Figure6_integrated_LOXL4_evidence")


def main():
    mkdirs()
    download_geo_files()
    all_expr, all_meta, all_de, all_rocs, all_corr, all_gsea = {}, {}, {}, {}, {}, {}
    for gse in GSES:
        expr, meta = prepare_dataset(gse)
        if TARGET not in expr.index:
            raise RuntimeError(f"{TARGET} not found in {gse}")
        de = differential_expression(gse, expr, meta)
        roc_df, _ = roc_analysis(gse, expr, meta)
        corr = correlation_analysis(gse, expr, meta)
        gsea = preranked_gsea(gse, de)
        all_expr[gse] = expr
        all_meta[gse] = meta
        all_de[gse] = de
        all_rocs[gse] = roc_df
        all_corr[gse] = corr
        all_gsea[gse] = gsea
    plot_outputs(all_expr, all_meta, all_de, all_rocs, all_corr, all_gsea)
    print("Analysis complete.")
    print(f"Results: {RESULTS}")
    print(f"Figures: {FIGURES}")


if __name__ == "__main__":
    main()
