"""Format the generated manuscript for PeerJ and embed all publication figures."""
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs/submission/Periplasmic_surface_charge_manuscript.docx"
OUTPUT = ROOT / "docs/submission/Acidophile_PeerJ_Line_Numbered_Manuscript.docx"
FIGURES = ROOT / "results/broad2/figures"
NEW_URL = "https://github.com/ajayvellanki0720-boop/acidophile-ph-limits-submission"
ORCID = "0009-0008-7994-152X"


def next_heading_one_index(paragraphs, start):
    for i in range(start + 1, len(paragraphs)):
        if paragraphs[i].style.name == "Heading 1":
            return i
    return len(paragraphs)


def remove_paragraph(paragraph):
    paragraph._element.getparent().remove(paragraph._element)


def add_figure_after(doc, anchor, image_name, label, caption_text, alt_text, width=6.0):
    image_p = doc.add_paragraph()
    image_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    image_p.paragraph_format.space_before = Pt(6)
    image_p.paragraph_format.space_after = Pt(4)
    image_p.paragraph_format.keep_with_next = True
    shape = image_p.add_run().add_picture(str(FIGURES / image_name), width=Inches(width))
    shape._inline.docPr.set("descr", alt_text)
    shape._inline.docPr.set("title", label)
    anchor._p.addnext(image_p._p)

    caption = doc.add_paragraph()
    caption.paragraph_format.line_spacing = 1.15
    caption.paragraph_format.space_after = Pt(8)
    caption.paragraph_format.keep_together = True
    lead = caption.add_run(label + " ")
    lead.bold = True
    caption.add_run(caption_text)
    for run in caption.runs:
        run.font.name = "Times New Roman"
        run.font.size = Pt(10)
        run.font.color.rgb = RGBColor(0, 0, 0)
    image_p._p.addnext(caption._p)
    return caption


doc = Document(SOURCE)

# Remove the journal-specific Impact Statement and move Data availability to the declarations block.
paragraphs = doc.paragraphs
impact_i = next(i for i, p in enumerate(paragraphs) if p.text.strip() == "Impact Statement")
impact_end = next_heading_one_index(paragraphs, impact_i)
for p in paragraphs[impact_i:impact_end]:
    remove_paragraph(p)

paragraphs = doc.paragraphs
data_i = next(i for i, p in enumerate(paragraphs) if p.text.strip() == "Data Summary")
data_end = next_heading_one_index(paragraphs, data_i)
data_block = paragraphs[data_i:data_end]
data_block[0].text = "Data availability"
conflicts = next(p for p in paragraphs if p.text.strip() == "Conflicts of interest")
conflicts.text = "Competing interests"
for p in data_block:
    conflicts._p.addprevious(p._p)

# Verified author identity and requested front-matter cleanup.
author = next(p for p in doc.paragraphs if p.text.strip().startswith("Authors: Ajay Vellanki."))
author.text = f"Authors: Ajay Vellanki. ORCID iD: {ORCID}"
running_title = next(p for p in doc.paragraphs if p.text.strip().startswith("Running title:"))
remove_paragraph(running_title)

for p in doc.paragraphs:
    if "https://github.com/ajayvellanki0720-boop/acidophile-ph-limits" in p.text:
        for run in p.runs:
            run.text = run.text.replace(
                "https://github.com/ajayvellanki0720-boop/acidophile-ph-limits", NEW_URL
            )
            # The source-generation step may already have inserted NEW_URL. Since
            # the legacy URL is a prefix of it, normalize the accidental doubled
            # suffix produced by the broad replacement above.
            run.text = run.text.replace(NEW_URL + "-submission", NEW_URL)
    if "(Fig. 4)." in p.text:
        for run in p.runs:
            run.text = run.text.replace("(Fig. 4).", "(Fig. 4b).")

# Page system and continuous Word line numbering.
section = doc.sections[0]
section.left_margin = Inches(1.15)
section.right_margin = Inches(1.0)
section.top_margin = Inches(1.0)
section.bottom_margin = Inches(1.0)
ln = section._sectPr.find(qn("w:lnNumType"))
if ln is None:
    raise RuntimeError("native line numbering disappeared from source manuscript")
ln.set(qn("w:countBy"), "1")
ln.set(qn("w:start"), "1")
ln.set(qn("w:restart"), "continuous")
ln.set(qn("w:distance"), "432")

# Reference-derived typography.
normal = doc.styles["Normal"]
normal.font.name = "Times New Roman"
normal.font.size = Pt(12)
normal.font.color.rgb = RGBColor(0, 0, 0)
normal.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
normal.paragraph_format.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
normal.paragraph_format.space_after = Pt(12)

for style_name, size in (("Heading 1", 14), ("Heading 2", 12)):
    style = doc.styles[style_name]
    style.font.name = "Times New Roman"
    style.font.size = Pt(size)
    style.font.bold = True
    style.font.color.rgb = RGBColor(0, 0, 0)
    style.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    style.paragraph_format.space_before = Pt(12)
    style.paragraph_format.space_after = Pt(4)
    style.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    style.paragraph_format.keep_with_next = True

for p in doc.paragraphs:
    if p.style.name == "Normal":
        p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
        p.paragraph_format.space_after = Pt(12)
    for run in p.runs:
        run.font.name = "Times New Roman"
        run._element.get_or_add_rPr().rFonts.set(qn("w:eastAsia"), "Times New Roman")
        run.font.color.rgb = RGBColor(0, 0, 0)

title = doc.paragraphs[0]
title.alignment = WD_ALIGN_PARAGRAPH.CENTER
title.paragraph_format.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
title.paragraph_format.space_after = Pt(12)
for run in title.runs:
    run.bold = True
    run.font.size = Pt(16)

for p in doc.paragraphs[1:6]:
    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
    p.paragraph_format.space_after = Pt(4)
doc.paragraphs[5].paragraph_format.space_after = Pt(12)

abstract = next(p for p in doc.paragraphs if p.text.strip() == "Abstract")
abstract.style = doc.styles["Normal"]
abstract.paragraph_format.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
abstract.paragraph_format.space_after = Pt(4)
abstract.paragraph_format.keep_with_next = True
for run in abstract.runs:
    run.bold = True
    run.font.size = Pt(12)

# Remove the detached legend section after retaining its scientific wording below.
legend_heading = next(p for p in doc.paragraphs if p.text.strip() == "Figure legends")
legend_i = next(i for i, p in enumerate(doc.paragraphs) if p.text.strip() == "Figure legends")
for p in list(doc.paragraphs[legend_i:]):
    remove_paragraph(p)

# Embed each panel near the section that interprets it.
anchor1 = next(p for p in doc.paragraphs if p.text.strip().startswith("The composite of net acidic residue share"))
add_figure_after(
    doc, anchor1, "fig1_compartment_contrast.png", "Fig. 1.",
    "The association is confined to the periplasmic-facing surface. Coefficients and 95% confidence intervals are shown for the net acidic composite computed on four position classes of the same 32 reference proteins in the same phylogenetic model, with genome-wide composition and source controlled. n = 2,768 except the soluble-surface control (n = 2,764). Positive values indicate a higher minimum growth pH and therefore less acid tolerance.",
    "Forest plot comparing exposed periplasmic, exposed cytoplasmic, transmembrane, and soluble-surface coefficients. Only the exposed periplasmic coefficient is clearly positive; the soluble-surface control is negative.",
)

anchor2 = next(p for p in doc.paragraphs if p.text.strip().startswith("On the genomes whose optimum is known"))
add_figure_after(
    doc, anchor2, "fig2_optimum_decomposition.png", "Fig. 2.",
    "Part of the association is independent of growth optimum and exposed histidine. Coefficients and 95% confidence intervals are shown for the periplasmic composite under four nested models, all fitted on the same 1,655 genomes for which the growth optimum is known. Percentages give the share of the composition-only coefficient retained.",
    "Forest plot of four nested models. The periplasmic surface-composition coefficient remains positive after controlling for growth optimum and histidine, retaining 39 percent of its composition-only magnitude.",
)

anchor3 = next(p for p in doc.paragraphs if p.text.strip().startswith("How many of the panel's reference families"))
add_figure_after(
    doc, anchor3, "fig3_per_reference_bh.png", "Fig. 3.",
    "Detectable presence of particular reference families predicts the limit. Per-reference coefficients and 95% confidence intervals for detection are sorted by effect; the seven families surviving Benjamini-Hochberg correction at q < 0.10 are shown with filled markers and the remainder with open markers. Every fit controls composition, source, mean alignment score, log genome size and assembly level.",
    "Forest plot of 22 reference-family detection coefficients. Seven coefficients survive multiple-testing correction, with effects in both directions and the largest negative effect for cytochrome bo3 ubiquinol oxidase subunit 1.",
    width=5.7,
)

anchor4 = next(p for p in doc.paragraphs if p.text.strip().startswith("The reported minima do carry a testing-grid signature"))
cap4a = add_figure_after(
    doc, anchor4, "fig4a_censoring_sensitivity.png", "Fig. 4a.",
    "Sensitivity of the headline coefficient to differential censoring aligned with the surface-composition predictor. The fitted association reaches zero only when the censoring displacement changes by 0.185 pH units per standard deviation of the predictor.",
    "Line plot showing the headline coefficient declining as differential censoring increases. The coefficient reaches zero at 0.185 pH units of censoring displacement per predictor standard deviation.",
)
add_figure_after(
    doc, cap4a, "fig4b_censoring_headroom.png", "Fig. 4b.",
    "No observable censoring pattern has enough physical headroom to account for the association. For each threatening proxy, the orange bar is the mean displacement required to drive the coefficient to zero and the grey bar is the displacement physically available before every true limit would fall at pH 0.",
    "Horizontal bars compare the censoring depth required to erase the result with the depth physically available. Both threatening proxies require more depth than exists.",
    width=5.7,
)

doc.core_properties.title = (
    "Periplasmic surface charge composition predicts the lower pH limit of bacterial growth independently of growth optimum"
)
doc.core_properties.subject = "PeerJ submission manuscript"
doc.core_properties.author = "Ajay Vellanki"
from peerj_revision import revise
revise(doc)
doc.save(OUTPUT)
print(OUTPUT)
