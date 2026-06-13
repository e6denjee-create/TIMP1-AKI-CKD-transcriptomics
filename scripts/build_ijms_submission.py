"""Build the IJMS manuscript, supplements, and public-repository payload."""

from __future__ import annotations

import csv
import json
import re
import shutil
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import PageBreak, Paragraph, SimpleDocTemplate, Table, TableStyle


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DOCX = Path(r"C:\Users\dell\Downloads\TIMP1_manuscript_submission_ready.docx")
OUT = ROOT / "submission_package" / "ijms_submission_2026-06-13"
REPO = ROOT / "public_repository" / "TIMP1_AKI_CKD_transcriptomics"
ZENODO = ROOT / "public_repository" / "TIMP1_AKI_CKD_transcriptomics_zenodo_payload"

TITLE = (
    "Cross-Cohort Transcriptomic Analysis Associates TIMP1 with Tubular "
    "Injury-Repair, Extracellular Matrix Remodeling, and a Conserved "
    "Disease-State Program in Human Kidney Disease"
)
AFFILIATION = (
    "Shanxi Bethune Hospital, Shanxi Academy of Medical Sciences, Third Hospital "
    "of Shanxi Medical University, Tongji Shanxi Hospital, Taiyuan 030032, China"
)
ADDRESS = (
    "No. 99 Longcheng Avenue, Xiaodian District, Taiyuan City, "
    "Shanxi Province 030032, China"
)
REPOSITORY_PLACEHOLDER = "Repository URL and archival DOI will be inserted after public release."

RECENT_PMIDS = [
    "39298548",
    "36142787",
    "32571916",
    "39110788",
    "36932062",
    "33176333",
    "35491858",
    "36265491",
    "38513647",
    "35064106",
    "38977708",
    "35922662",
    "38126209",
    "35513123",
    "41286490",
    "35709763",
    "38351505",
    "37130011",
    "26768243",
    "39113797",
    "32303283",
    "33073587",
    "31000567",
    "39016438",
    "39572154",
    "38609039",
    "34560077",
    "38086793",
    "36229672",
    "33411069",
]

CLASSIC_REFERENCES = [
    "Venkatachalam, M.A.; Weinberg, J.M.; Kriz, W.; Bidani, A.K. Failed Tubule Recovery, AKI-CKD Transition, and Kidney Disease Progression. J. Am. Soc. Nephrol. 2015, 26, 1765-1776. https://doi.org/10.1681/ASN.2015010006.",
    "Ferenbach, D.A.; Bonventre, J.V. Mechanisms of Maladaptive Repair after AKI Leading to Accelerated Kidney Ageing and CKD. Nat. Rev. Nephrol. 2015, 11, 264-276. https://doi.org/10.1038/nrneph.2015.3.",
    "Bonventre, J.V.; Yang, L. Cellular Pathophysiology of Ischemic Acute Kidney Injury. J. Clin. Investig. 2011, 121, 4210-4221. https://doi.org/10.1172/JCI45161.",
    "Humphreys, B.D.; Valerius, M.T.; Kobayashi, A.; Mugford, J.W.; Soeung, S.; Duffield, J.S.; McMahon, A.P.; Bonventre, J.V. Intrinsic Epithelial Cells Repair the Kidney after Injury. Cell Stem Cell 2008, 2, 284-291. https://doi.org/10.1016/j.stem.2008.01.014.",
    "Liu, Y. Cellular and Molecular Mechanisms of Renal Fibrosis. Nat. Rev. Nephrol. 2011, 7, 684-696. https://doi.org/10.1038/nrneph.2011.149.",
    "Yang, L.; Besschetnova, T.Y.; Brooks, C.R.; Shah, J.V.; Bonventre, J.V. Epithelial Cell Cycle Arrest in G2/M Mediates Kidney Fibrosis after Injury. Nat. Med. 2010, 16, 535-543. https://doi.org/10.1038/nm.2144.",
    "Canaud, G.; Brooks, C.R.; Kishi, S.; Taguchi, K.; Nishimura, K.; Magassa, S.; Scott, A.; Hsiao, L.L.; Ichimura, T.; Terzi, F.; et al. Cyclin G1 and TASCC Regulate Kidney Epithelial Cell G2-M Arrest and Fibrotic Maladaptive Repair. Sci. Transl. Med. 2019, 11, eaav4754. https://doi.org/10.1126/scitranslmed.aav4754.",
    "Rabb, H.; Griffin, M.D.; McKay, D.B.; Swaminathan, S.; Pickkers, P.; Rosner, M.H.; Kellum, J.A.; Ronco, C. Inflammation in AKI: Current Understanding, Key Questions, and Knowledge Gaps. J. Am. Soc. Nephrol. 2016, 27, 371-379. https://doi.org/10.1681/ASN.2015030261.",
    "Brew, K.; Nagase, H. The Tissue Inhibitors of Metalloproteinases (TIMPs): An Ancient Family with Structural and Functional Diversity. Biochim. Biophys. Acta 2010, 1803, 55-71. https://doi.org/10.1016/j.bbamcr.2010.01.003.",
    "Lake, B.B.; Menon, R.; Winfree, S.; Hu, Q.; Melo Ferreira, R.; Kalhor, K.; Barwinska, D.; Otto, E.A.; Ferkowicz, M.; Diep, D.; et al. An Atlas of Healthy and Injured Cell States and Niches in the Human Kidney. Nature 2023, 619, 585-594. https://doi.org/10.1038/s41586-023-05769-3.",
    "Meng, X.M.; Nikolic-Paterson, D.J.; Lan, H.Y. TGF-beta: The Master Regulator of Fibrosis. Nat. Rev. Nephrol. 2016, 12, 325-338. https://doi.org/10.1038/nrneph.2016.48.",
]


def clean_text(value: object) -> str:
    return str(value).replace("¦Ā", "beta").replace("β", "beta").strip()


def format_recent_reference(row: pd.Series) -> str:
    author_items = [x.strip() for x in str(row["authors"]).split(";") if x.strip()]
    formatted = []
    for item in author_items[:10]:
        bits = item.split()
        if len(bits) > 1:
            formatted.append(f"{bits[0]}, {' '.join(bits[1:])}")
        else:
            formatted.append(item)
    authors = "; ".join(formatted)
    if len(author_items) > 10:
        authors += "; et al."
    title = clean_text(row["title"]).rstrip(".")
    journal = clean_text(row["journal_abbrev"])
    year = str(row["year"])
    volume = "" if pd.isna(row["volume"]) else clean_text(row["volume"])
    pages = "" if pd.isna(row["pages_or_elocation"]) else clean_text(row["pages_or_elocation"])
    doi = "" if pd.isna(row["doi"]) else clean_text(row["doi"])
    details = ", ".join(x for x in [year, volume, pages] if x)
    suffix = f" https://doi.org/{doi}." if doi else f" PMID: {row['pmid']}."
    return f"{authors}. {title}. {journal} {details}.{suffix}"


def load_references() -> list[str]:
    refs = list(CLASSIC_REFERENCES)
    library = pd.read_csv(OUT / "references_verified.csv", dtype={"pmid": str})
    by_pmid = library.set_index("pmid")
    for pmid in RECENT_PMIDS:
        if pmid in by_pmid.index:
            refs.append(format_recent_reference(by_pmid.loc[pmid]))
    return refs


def set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def add_table(doc: Document, headers: list[str], rows: list[list[object]]) -> None:
    table = doc.add_table(rows=1, cols=len(headers))
    table.style = "Table Grid"
    table.autofit = True
    for idx, value in enumerate(headers):
        cell = table.rows[0].cells[idx]
        cell.text = value
        set_cell_shading(cell, "D9EAF7")
        cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
        for run in cell.paragraphs[0].runs:
            run.bold = True
            run.font.size = Pt(8)
    for row in rows:
        cells = table.add_row().cells
        for idx, value in enumerate(row):
            cells[idx].text = str(value)
            cells[idx].vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            for paragraph in cells[idx].paragraphs:
                for run in paragraph.runs:
                    run.font.size = Pt(8)
    doc.add_paragraph()


def section_paragraphs(source: Document, start: str, end: str | None) -> list[tuple[str, str]]:
    collecting = False
    items: list[tuple[str, str]] = []
    for paragraph in source.paragraphs:
        if paragraph.text == start:
            collecting = True
        if collecting:
            if end and paragraph.text == end:
                break
            items.append((paragraph.style.name, paragraph.text))
    return items


def add_section_content(doc: Document, items: list[tuple[str, str]], citation_updates: dict[str, str]) -> None:
    for style, original in items:
        text_value = citation_updates.get(original, original)
        if not text_value.strip():
            continue
        if style.startswith("Heading"):
            level = int(style.split()[-1])
            if text_value == "Methods":
                text_value = "Materials and Methods"
            doc.add_heading(text_value, level=level)
        else:
            paragraph = doc.add_paragraph(text_value)
            paragraph.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
            if original.startswith("Table 1. Public human kidney transcriptomic cohorts"):
                add_table(
                    doc,
                    ["Cohort", "Role", "Disease samples", "Control samples", "Tissue/platform", "Primary analytic use"],
                    [
                        ["GSE139061", "Discovery", "39 AKI", "9 reference nephrectomies", "Human kidney biopsy RNA-seq", "Expression, signatures, module discovery"],
                        ["GSE30718", "Discovery", "28 AKI/transplant injury", "11 controls", "Affymetrix microarray", "Expression, signatures, module discovery"],
                        ["GSE66494", "Discovery", "53 CKD/fibrosis", "8 controls", "Agilent microarray", "Expression, signatures, module discovery"],
                        ["GSE180394", "External validation", "44 kidney disease", "9 living donors; 6 tumor-nephrectomy sensitivity controls", "Microdissected tubules; Affymetrix", "External validation and robustness"],
                    ],
                )


def build_manuscript() -> None:
    source = Document(SOURCE_DOCX)
    refs = load_references()
    doc = Document()
    section = doc.sections[0]
    section.top_margin = Inches(0.8)
    section.bottom_margin = Inches(0.8)
    section.left_margin = Inches(0.9)
    section.right_margin = Inches(0.9)
    styles = doc.styles
    styles["Normal"].font.name = "Arial"
    styles["Normal"].font.size = Pt(10)
    for name in ["Heading 1", "Heading 2"]:
        styles[name].font.name = "Arial"
        styles[name].font.color.rgb = None
    styles["Heading 1"].font.size = Pt(14)
    styles["Heading 2"].font.size = Pt(12)

    title = doc.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = title.add_run(TITLE)
    run.bold = True
    run.font.name = "Arial"
    run.font.size = Pt(16)
    authors = doc.add_paragraph()
    authors.alignment = WD_ALIGN_PARAGRAPH.CENTER
    authors.add_run("Yanzhao Ji 1,* and Zhihong Gao 1").bold = True
    aff = doc.add_paragraph()
    aff.alignment = WD_ALIGN_PARAGRAPH.CENTER
    aff.add_run(f"1 {AFFILIATION}")
    corr = doc.add_paragraph()
    corr.alignment = WD_ALIGN_PARAGRAPH.CENTER
    corr.add_run(f"* Correspondence: jiyanzhao@sxbqeh.com.cn; {ADDRESS}")

    abstract = section_paragraphs(source, "Abstract", "Introduction")
    abstract_updates = {
        "Background: Maladaptive repair after kidney injury is accompanied by persistent tubular stress, inflammation, cellular senescence, and extracellular matrix (ECM) remodeling. Tissue inhibitor of metalloproteinases 1 (TIMP1) responds to tissue injury in diverse contexts, but its reproducible transcriptomic context across human kidney disease cohorts has not been systematically defined.": (
            "Maladaptive repair after kidney injury is accompanied by persistent tubular stress, "
            "inflammation, cellular senescence, and extracellular matrix (ECM) remodeling. We "
            "performed an integrative, platform-specific transcriptomic association study using "
            "three discovery cohorts and an independent microdissected tubular cohort. TIMP1 was "
            "directionally elevated in all discovery cohorts and was higher in GSE180394 disease "
            "samples than in living-donor controls. Nine prespecified signatures, calculated after "
            "excluding TIMP1, showed reproducible positive disease-only correlations with TIMP1. "
            "A 421-gene stringent disease-only module was externally reproduced, and all 12 core "
            "genes retained positive external correlations. Leave-one-diagnosis-out, alternative "
            "control, and diagnosis-adjusted analyses supported robustness. These findings associate "
            "TIMP1 with tubular injury-repair, ECM remodeling, inflammation, senescence, and fibrosis-"
            "related programs across human kidney disease cohorts. TIMP1 is nominated as a candidate "
            "component of a broader maladaptive injury-response state; the cross-sectional design "
            "does not establish prediction, diagnostic utility, kidney specificity, or causality."
        ),
    }
    doc.add_heading("Abstract", level=1)
    doc.add_paragraph(abstract_updates[next(iter(abstract_updates))])
    doc.add_paragraph(
        "Keywords: TIMP1; acute kidney injury; chronic kidney disease; tubular injury; "
        "maladaptive repair; extracellular matrix remodeling; fibrosis; transcriptomics"
    )

    citation_updates = {
        source.paragraphs[15].text: source.paragraphs[15].text.replace("[1,2]", "[1,2,12,13,29-32]").replace("[2-4]", "[2-4,14,15]").replace("[2,5]", "[2,5,16-18]"),
        source.paragraphs[16].text: source.paragraphs[16].text.replace("[6,7]", "[6,7,18,19,33,39,40]").replace("[8]", "[8,20,23,24]"),
        source.paragraphs[17].text: source.paragraphs[17].text + " Recent single-cell and spatial atlases further show that injury-associated epithelial, immune, endothelial, and stromal states are spatially and molecularly heterogeneous [10,14-20,26-28,34].",
        source.paragraphs[18].text: source.paragraphs[18].text.replace("[9]", "[9,25]") + " Experimental evidence linking TIMP1 expression with kidney fibrosis susceptibility supports investigation of this association while not proving direct causality [25].",
        source.paragraphs[19].text: source.paragraphs[19].text + " The overall validation workflow is summarized in Figure 1.",
        source.paragraphs[31].text: source.paragraphs[31].text.replace("Supplementary Table S1", "Supplementary Table S3"),
        source.paragraphs[68].text: source.paragraphs[68].text.replace("[1-4]", "[1-4,12-15,29-33]").replace("[6,7]", "[6,7,18,19,23,24,39,40]"),
        source.paragraphs[69].text: source.paragraphs[69].text.replace("[9]", "[9,11,16,21,22,25,35-38]"),
        source.paragraphs[71].text: source.paragraphs[71].text + " Comparable compartment-resolved kidney atlases demonstrate the value of validating such programs across epithelial and stromal contexts [10,14-20,26-28,34].",
        source.paragraphs[74].text: source.paragraphs[74].text.replace("[10]", "[10,14-20,26-28,34]"),
    }

    add_section_content(doc, section_paragraphs(source, "Introduction", "Methods"), citation_updates)
    add_section_content(doc, section_paragraphs(source, "Results", "Discussion"), citation_updates)
    add_section_content(doc, section_paragraphs(source, "Discussion", "Conclusions"), citation_updates)
    add_section_content(doc, section_paragraphs(source, "Methods", "Results"), citation_updates)
    add_section_content(doc, section_paragraphs(source, "Conclusions", "Declarations"), citation_updates)

    doc.add_heading("Supplementary Materials", level=1)
    doc.add_paragraph(
        "The following supporting information can be downloaded with the article: "
        "Table S1, cohort characteristics and analytic roles; Table S2, TIMP1 expression "
        "statistics; Table S3, prespecified signatures and TIMP1-signature correlations; "
        "Table S4, module and core-gene validation; Table S5, robustness analyses, feature "
        "coverage, and missing-data log."
    )
    doc.add_heading("Author Contributions", level=1)
    doc.add_paragraph(
        "Conceptualization, Y.J.; methodology, Y.J. and Z.G.; software, Y.J.; validation, "
        "Z.G.; formal analysis, Y.J.; investigation, Y.J. and Z.G.; data curation, Y.J. and "
        "Z.G.; visualization, Y.J.; writing-original draft preparation, Y.J.; writing-review "
        "and editing, Y.J. and Z.G.; supervision, Y.J.; project administration, Y.J.; funding "
        "acquisition, Y.J. All authors have read and agreed to the published version of the manuscript."
    )
    doc.add_heading("Funding", level=1)
    doc.add_paragraph(
        "This research was supported by the Fundamental Research Program of Shanxi Province "
        "(No. 202403021222407) and the Shanxi Bethune Hospital Project (No. 2023RC06)."
    )
    doc.add_heading("Institutional Review Board Statement", level=1)
    doc.add_paragraph(
        "Ethical review and approval were waived for this study because it was a secondary "
        "analysis of publicly available, de-identified transcriptomic datasets and involved "
        "no new participant recruitment or biospecimen collection."
    )
    doc.add_heading("Informed Consent Statement", level=1)
    doc.add_paragraph("Not applicable.")
    doc.add_heading("Data Availability Statement", level=1)
    doc.add_paragraph(
        "The public transcriptomic datasets are available from the NCBI Gene Expression "
        "Omnibus under accession numbers GSE139061, GSE30718, GSE66494, GSE180394, "
        "GSE210622, and GSE267242. Analysis code, derived tables, intermediate outputs, "
        "session information, source data for figures, and publication-quality figures "
        f"will be archived in a public GitHub repository and Zenodo record. {REPOSITORY_PLACEHOLDER}"
    )
    doc.add_heading("Acknowledgments", level=1)
    doc.add_paragraph(
        "The authors thank the investigators and participants of the public GEO studies "
        "analyzed in this work."
    )
    doc.add_heading("Conflicts of Interest", level=1)
    doc.add_paragraph("The authors declare no conflict of interest.")

    doc.add_heading("References", level=1)
    for index, reference in enumerate(refs, start=1):
        paragraph = doc.add_paragraph(f"{index}. {reference}")
        paragraph.paragraph_format.first_line_indent = Inches(-0.2)
        paragraph.paragraph_format.left_indent = Inches(0.2)
        paragraph.paragraph_format.space_after = Pt(3)

    doc.add_heading("Figure Legends", level=1)
    for paragraph in source.paragraphs[107:112]:
        doc.add_paragraph(paragraph.text)

    footer = doc.sections[0].footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.CENTER
    footer.add_run("TIMP1 integrative transcriptomic analysis - IJMS submission draft")
    doc.save(OUT / "TIMP1_IJMS_Manuscript.docx")


def dataframe_with_source(df: pd.DataFrame, source: str) -> pd.DataFrame:
    output = df.copy()
    output.insert(0, "source_file", source)
    return output


def build_supplement() -> dict[str, list[tuple[str, pd.DataFrame]]]:
    cohorts = pd.DataFrame(
        [
            ["GSE139061", "Discovery", "39 AKI", "9 controls", "Kidney biopsy RNA-seq", "Expression, signatures, module discovery"],
            ["GSE30718", "Discovery", "28 AKI/transplant injury", "11 controls", "Affymetrix microarray", "Expression, signatures, module discovery"],
            ["GSE66494", "Discovery", "53 CKD/fibrosis", "8 controls", "Agilent microarray", "Expression, signatures, module discovery"],
            ["GSE180394", "External validation", "44 kidney disease", "9 living donors; 6 tumor-nephrectomy sensitivity controls", "Microdissected tubules; Affymetrix", "External validation and robustness"],
            ["GSE210622", "Exploratory localization", "One locally available AKI donor", "Not applicable", "Single-cell RNA-seq", "Exploratory localization only"],
        ],
        columns=["dataset", "role", "disease_samples", "control_samples", "tissue_platform", "analytic_use"],
    )
    s1 = [("Cohort characteristics", cohorts)]

    bulk = pd.read_csv(ROOT / "results/timp1_validation/TIMP1_bulk_validation_statistics.csv")
    external = pd.read_csv(ROOT / "results/timp1_validation/external_GSE180394_TIMP1_group_statistics.csv")
    controls = pd.read_csv(ROOT / "results/timp1_validation/external_GSE180394_sensitivity_control_statistics.csv")
    s2 = [
        ("Discovery expression statistics", dataframe_with_source(bulk, "TIMP1_bulk_validation_statistics.csv")),
        ("External expression statistics", dataframe_with_source(external, "external_GSE180394_TIMP1_group_statistics.csv")),
        ("Alternative control sensitivity", dataframe_with_source(controls, "external_GSE180394_sensitivity_control_statistics.csv")),
    ]

    gene_sets = pd.read_csv(ROOT / "results/timp1_validation/signature_gene_sets_used.csv")
    correlations = pd.read_csv(ROOT / "results/timp1_validation/TIMP1_signature_correlations.csv")
    external_corr = pd.read_csv(ROOT / "results/timp1_validation/external_GSE180394_signature_correlations.csv")
    s3 = [
        ("Prespecified signature membership", dataframe_with_source(gene_sets, "signature_gene_sets_used.csv")),
        ("Discovery correlations", dataframe_with_source(correlations, "TIMP1_signature_correlations.csv")),
        ("External correlations", dataframe_with_source(external_corr, "external_GSE180394_signature_correlations.csv")),
    ]

    core = pd.read_csv(ROOT / "results/timp1_validation/core_module_genes_summary.csv")
    stringent = pd.read_csv(ROOT / "results/timp1_validation/stringent_TIMP1_correlated_module.csv")
    module_external = pd.read_csv(ROOT / "results/timp1_validation/external_GSE180394_stringent_module_correlations.csv")
    core_external = pd.read_csv(ROOT / "results/timp1_validation/external_GSE180394_core_module_correlations.csv")
    s4 = [
        ("Core module summary", dataframe_with_source(core, "core_module_genes_summary.csv")),
        ("Stringent 421-gene module", dataframe_with_source(stringent, "stringent_TIMP1_correlated_module.csv")),
        ("External module validation", dataframe_with_source(module_external, "external_GSE180394_stringent_module_correlations.csv")),
        ("External core-gene validation", dataframe_with_source(core_external, "external_GSE180394_core_module_correlations.csv")),
    ]

    lodo = pd.read_csv(ROOT / "results/timp1_validation/external_GSE180394_sensitivity_leave_one_diagnosis_out_summary.csv")
    robust = pd.read_csv(ROOT / "results/timp1_validation/external_GSE180394_sensitivity_robust_correlations.csv")
    coverage = pd.read_csv(ROOT / "results/timp1_validation/external_GSE180394_input_feature_coverage.csv")
    missing = pd.read_csv(ROOT / "results/timp1_validation/missing_data_log.csv")
    s5 = [
        ("Leave-one-diagnosis-out summary", dataframe_with_source(lodo, "external_GSE180394_sensitivity_leave_one_diagnosis_out_summary.csv")),
        ("Diagnosis-adjusted robustness", dataframe_with_source(robust, "external_GSE180394_sensitivity_robust_correlations.csv")),
        ("External feature coverage", dataframe_with_source(coverage, "external_GSE180394_input_feature_coverage.csv")),
        ("Missing-data log", dataframe_with_source(missing, "missing_data_log.csv")),
    ]
    sections = {
        "S1_Cohorts": s1,
        "S2_Expression": s2,
        "S3_Signatures": s3,
        "S4_Modules": s4,
        "S5_Sensitivity": s5,
    }

    workbook = Workbook()
    workbook.remove(workbook.active)
    for sheet_name, tables in sections.items():
        sheet = workbook.create_sheet(sheet_name)
        sheet.sheet_view.showGridLines = False
        row = 1
        for title, frame in tables:
            sheet.cell(row=row, column=1, value=title)
            sheet.cell(row=row, column=1).font = Font(bold=True, size=13, color="FFFFFF")
            sheet.cell(row=row, column=1).fill = PatternFill("solid", fgColor="1F4E78")
            sheet.merge_cells(start_row=row, start_column=1, end_row=row, end_column=max(1, len(frame.columns)))
            row += 1
            for col_idx, column in enumerate(frame.columns, start=1):
                cell = sheet.cell(row=row, column=col_idx, value=str(column))
                cell.font = Font(bold=True, color="FFFFFF")
                cell.fill = PatternFill("solid", fgColor="5B9BD5")
                cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
            row += 1
            for values in frame.itertuples(index=False, name=None):
                for col_idx, value in enumerate(values, start=1):
                    if pd.isna(value):
                        value = None
                    cell = sheet.cell(row=row, column=col_idx, value=value)
                    cell.alignment = Alignment(vertical="top", wrap_text=True)
                row += 1
            row += 2
        sheet.freeze_panes = "A3"
        sheet.auto_filter.ref = sheet.dimensions
        for col_idx in range(1, sheet.max_column + 1):
            max_len = max(len(str(sheet.cell(r, col_idx).value or "")) for r in range(1, min(sheet.max_row, 250) + 1))
            sheet.column_dimensions[get_column_letter(col_idx)].width = min(max(max_len + 2, 11), 34)
    xlsx = OUT / "Supplementary_Tables_S1-S5.xlsx"
    workbook.save(xlsx)

    styles = getSampleStyleSheet()
    title_style = ParagraphStyle("SuppTitle", parent=styles["Heading1"], alignment=TA_CENTER, fontSize=14)
    cell_style = ParagraphStyle("Cell", parent=styles["BodyText"], fontSize=5.5, leading=6.3)
    pdf = SimpleDocTemplate(
        str(OUT / "Supplementary_Tables_S1-S5.pdf"),
        pagesize=landscape(A4),
        leftMargin=8 * mm,
        rightMargin=8 * mm,
        topMargin=8 * mm,
        bottomMargin=8 * mm,
    )
    story = []
    for sheet_name, tables in sections.items():
        story.append(Paragraph(sheet_name.replace("_", " "), title_style))
        for table_title, frame in tables:
            story.append(Paragraph(table_title, styles["Heading2"]))
            display = frame.copy()
            max_rows = 80
            if len(display) > max_rows:
                display = display.head(max_rows)
                story.append(Paragraph(f"PDF preview shows first {max_rows} rows; complete table is in the Excel file.", styles["Italic"]))
            data = [[Paragraph(clean_text(c), cell_style) for c in display.columns]]
            for values in display.astype(object).where(pd.notnull(display), "").itertuples(index=False, name=None):
                data.append([Paragraph(clean_text(v), cell_style) for v in values])
            usable = landscape(A4)[0] - 16 * mm
            widths = [usable / max(1, len(display.columns))] * len(display.columns)
            table = Table(data, colWidths=widths, repeatRows=1)
            table.setStyle(
                TableStyle(
                    [
                        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1F4E78")),
                        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#B7C9D6")),
                        ("VALIGN", (0, 0), (-1, -1), "TOP"),
                        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F3F7FA")]),
                    ]
                )
            )
            story.append(table)
            story.append(Paragraph(" ", styles["BodyText"]))
        story.append(PageBreak())
    pdf.build(story)
    return sections


def build_cover_letter() -> None:
    doc = Document()
    doc.styles["Normal"].font.name = "Arial"
    doc.styles["Normal"].font.size = Pt(11)
    doc.add_paragraph("13 June 2026")
    doc.add_paragraph("Editorial Office\nInternational Journal of Molecular Sciences")
    doc.add_paragraph("Dear Editors,")
    doc.add_paragraph(
        f"We are pleased to submit the original article entitled \"{TITLE}\" for consideration "
        "in the International Journal of Molecular Sciences."
    )
    doc.add_paragraph(
        "This integrative transcriptomic study evaluates TIMP1 across three discovery kidney "
        "cohorts and an independent microdissected tubular cohort. The work emphasizes "
        "cross-dataset reproducibility, disease-only pathway associations, a conserved gene "
        "module, and robustness to diagnosis and control definitions. The manuscript treats "
        "TIMP1 as a candidate injury and extracellular-matrix-remodeling-associated gene and "
        "does not claim kidney specificity, diagnostic utility, or causality."
    )
    doc.add_paragraph(
        "The manuscript is original, is not under consideration elsewhere, and has been approved "
        "by both authors. The authors declare no conflict of interest. The analysis used publicly "
        "available de-identified datasets and did not involve new participant recruitment."
    )
    doc.add_paragraph("Thank you for considering our manuscript.")
    doc.add_paragraph(
        "Sincerely,\nYanzhao Ji, Corresponding Author\n"
        f"{AFFILIATION}\n{ADDRESS}\nEmail: jiyanzhao@sxbqeh.com.cn"
    )
    doc.save(OUT / "TIMP1_IJMS_Cover_Letter.docx")


def copy_tree_contents(source: Path, destination: Path, patterns: tuple[str, ...] | None = None) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        if item.is_file() and (patterns is None or item.suffix.lower() in patterns):
            shutil.copy2(item, destination / item.name)


def build_repository() -> None:
    if REPO.exists():
        shutil.rmtree(REPO)
    if ZENODO.exists():
        shutil.rmtree(ZENODO)
    REPO.mkdir(parents=True)
    ZENODO.mkdir(parents=True)

    copy_tree_contents(ROOT / "scripts", REPO / "scripts", (".r", ".py"))
    copy_tree_contents(ROOT / "results/timp1_validation", REPO / "results/timp1_validation")
    copy_tree_contents(ROOT / "figures/timp1_validation", REPO / "figures/timp1_validation", (".pdf", ".png"))
    copy_tree_contents(ROOT / "TIMP1_AKI_CKD_project/data/processed", REPO / "data/processed", (".gz", ".csv"))
    copy_tree_contents(ROOT / "TIMP1_AKI_CKD_project/data/metadata", REPO / "data/metadata", (".csv",))
    shutil.copy2(ROOT / "AGENTS.md", REPO / "ANALYSIS_GUARDRAILS.md")

    for folder in ["scripts", "results", "figures", "data"]:
        shutil.copytree(REPO / folder, ZENODO / folder)
    large_rds = ROOT / "TIMP1_AKI_CKD_project/data/processed/GSE210622_GSM6433706_seurat.rds"
    if large_rds.exists():
        target = ZENODO / "data/processed"
        target.mkdir(parents=True, exist_ok=True)
        shutil.copy2(large_rds, target / large_rds.name)

    readme = f"""# TIMP1 kidney injury transcriptomic analysis

This repository contains the reproducible code, derived tables, intermediate
outputs, session information, source data, and figures for:

**{TITLE}**

Authors: Yanzhao Ji and Zhihong Gao.

The analysis positions TIMP1 as a candidate kidney injury and extracellular
matrix remodeling-associated gene. Results are association-based and do not
establish diagnostic utility, kidney specificity, longitudinal prediction, or
causality.

## Contents

- `scripts/`: R and Python workflows.
- `results/timp1_validation/`: statistical tables, intermediate results,
  session information, interpretations, and missing-data records.
- `figures/timp1_validation/`: PDF and PNG figure files.
- `data/processed/`: processed expression matrices small enough for GitHub.
- `data/metadata/`: sample metadata.

The complete archive, including the >100 MB exploratory single-cell RDS object,
is intended for Zenodo because GitHub rejects individual files above 100 MB.

## Public datasets

NCBI GEO: GSE139061, GSE30718, GSE66494, GSE180394, GSE210622, and GSE267242.

## Reproduction

R 4.6.0 was used for the principal analyses. Run the scripts in numeric or
workflow order from the repository root. Package versions are recorded in the
`sessionInfo*.txt` files under `results/timp1_validation/`.

## Missing data

Unavailable optional datasets or resources are recorded in
`results/timp1_validation/missing_data_log.csv`; affected units are skipped
without silently terminating the remaining workflow.
"""
    (REPO / "README.md").write_text(readme, encoding="utf-8")
    (ZENODO / "README.md").write_text(readme, encoding="utf-8")
    license_text = (
        "MIT License\n\nCopyright (c) 2026 Yanzhao Ji and Zhihong Gao\n\n"
        "Permission is hereby granted, free of charge, to any person obtaining a copy "
        "of this software and associated documentation files (the \"Software\"), to deal "
        "in the Software without restriction, including without limitation the rights "
        "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies "
        "of the Software, subject to inclusion of this notice.\n\nTHE SOFTWARE IS PROVIDED "
        "\"AS IS\", WITHOUT WARRANTY OF ANY KIND."
    )
    (REPO / "LICENSE").write_text(license_text, encoding="utf-8")
    (ZENODO / "LICENSE").write_text(license_text, encoding="utf-8")
    citation = f"""cff-version: 1.2.0
message: "If you use this repository, please cite the archived release."
title: "{TITLE}"
type: software
authors:
  - family-names: Ji
    given-names: Yanzhao
    email: jiyanzhao@sxbqeh.com.cn
  - family-names: Gao
    given-names: Zhihong
version: 1.0.0
date-released: 2026-06-13
license: MIT
repository-code: "TO_BE_ASSIGNED"
"""
    (REPO / "CITATION.cff").write_text(citation, encoding="utf-8")
    (ZENODO / "CITATION.cff").write_text(citation, encoding="utf-8")
    zenodo_json = {
        "title": TITLE,
        "upload_type": "software",
        "description": "Reproducible code, derived data, intermediate outputs, and figures for an integrative TIMP1 kidney transcriptomic analysis.",
        "creators": [
            {"name": "Ji, Yanzhao", "affiliation": AFFILIATION},
            {"name": "Gao, Zhihong", "affiliation": AFFILIATION},
        ],
        "keywords": ["TIMP1", "acute kidney injury", "chronic kidney disease", "fibrosis", "transcriptomics"],
        "license": "mit",
    }
    (ZENODO / ".zenodo.json").write_text(json.dumps(zenodo_json, indent=2), encoding="utf-8")

    manifest_rows = []
    for path in sorted(ZENODO.rglob("*")):
        if path.is_file():
            manifest_rows.append(
                {
                    "relative_path": path.relative_to(ZENODO).as_posix(),
                    "size_bytes": path.stat().st_size,
                }
            )
    with (ZENODO / "FILE_MANIFEST.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["relative_path", "size_bytes"])
        writer.writeheader()
        writer.writerows(manifest_rows)


def build_readme() -> None:
    text = f"""IJMS submission package
Prepared: 2026-06-13

Files
- TIMP1_IJMS_Manuscript.docx
- TIMP1_IJMS_Cover_Letter.docx
- Supplementary_Tables_S1-S5.xlsx
- Supplementary_Tables_S1-S5.pdf
- references_verified.csv
- Figure_1 to Figure_5 in PDF and PNG formats

Authors
- Yanzhao Ji: first and corresponding author
- Zhihong Gao
- Affiliation: {AFFILIATION}
- Correspondence: jiyanzhao@sxbqeh.com.cn

Repository status
- A GitHub-ready package is under public_repository/TIMP1_AKI_CKD_transcriptomics.
- A complete Zenodo payload, including the large single-cell RDS object, is
  under public_repository/TIMP1_AKI_CKD_transcriptomics_zenodo_payload.
- The manuscript intentionally does not contain a fabricated URL or DOI.
"""
    (OUT / "README_Submission_Files.txt").write_text(text, encoding="utf-8")


def copy_figures() -> None:
    mapping = {
        "v4_manuscript_figure_1_workflow_cohort_overview": "Figure_1_TIMP1_Validation_Workflow",
        "v4_manuscript_figure_2_TIMP1_expression_cohorts": "Figure_2_TIMP1_Expression_Cohorts",
        "v4_manuscript_figure_3_signature_correlation_heatmap": "Figure_3_TIMP1_Signature_Correlations",
        "v4_manuscript_figure_4_GSE180394_stringent_module_scatter": "Figure_4_TIMP1_Stringent_Module_GSE180394",
        "v4_manuscript_figure_5_core_gene_correlation_heatmap": "Figure_5_TIMP1_Core_Gene_Correlations",
    }
    for source_name, target_name in mapping.items():
        for suffix in [".pdf", ".png"]:
            shutil.copy2(
                ROOT / "figures/timp1_validation" / f"{source_name}{suffix}",
                OUT / f"{target_name}{suffix}",
            )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    build_supplement()
    build_manuscript()
    build_cover_letter()
    copy_figures()
    build_repository()
    build_readme()
    print(f"Built IJMS package in {OUT}")
    print(f"Built GitHub package in {REPO}")
    print(f"Built Zenodo payload in {ZENODO}")


if __name__ == "__main__":
    main()
