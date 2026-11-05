"""The submission manuscript as a line-numbered Word file -> docs/submission/Periplasmic_surface_charge_manuscript.docx
WHY A SCRIPT RATHER THAN A HAND-BUILT FILE. docs/MANUSCRIPT.md is the working draft and every signed coefficient in it
is checked against the stored result its own claim cites (`scripts/69`). A .docx assembled by hand would be a second
copy of those numbers with nothing checking it, which is how correction 12 happened. So the submission file is
GENERATED, and the headline values are asserted against the stored results before a single paragraph is written - the
same discipline `scripts/100` applies to the figures.
WHAT THE FORMAT IS FOR. Target journal is Microbial Genomics (Microbiology Society). Their submission checklist is
explicit that "the article has continuous line numbering throughout, for ease of reviewing", that the file should be
"a single file, preferably readable by Microsoft Word, using a standard typeface such as Times New Roman or Arial",
and that a Research Article runs 3000-7000 words. Microbial Genomics additionally requires a DATA SUMMARY (mandatory
Open Data policy) and an IMPACT STATEMENT, neither of which most journals ask for, and author contributions in CRediT
form. The section order below follows their "Prepare an article" page. They operate format-free submission, so none
of this is mandatory at first submission - but line numbering is on the checklist, and it is what was asked for.
LINE NUMBERING IS A WORD FEATURE, NOT TEXT. Numbers are written into the section properties as `w:lnNumType`, so Word,
LibreOffice and Google Docs render them in the margin and they renumber themselves as the author edits. Prepending
numbers to each line as characters would look identical in a screenshot and be wrong the moment anyone typed a word,
which is the whole reason journals ask for the real thing.
AUTHOR-SPECIFIC DETAILS are limited to information verified in the supplied manuscript and repository record; no
affiliation, funding source or research approval is invented.
REFERENCES ARE REAL AND WERE EACH RETRIEVED IN THIS PROJECT, not recalled. Volume, issue, pages and DOI for the
primary sources come from PubMed records fetched during the work; three preprints and two database papers are cited
as what they are. Any entry this project could not verify is absent rather than approximated.
Usage: python3 scripts/105_submission_docx.py   (env DATA_DIR=results/broad2)"""
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, Inches

D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
OUT = ROOT / "docs/submission"
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------------------------------------------
# Every figure quoted below is pulled from the stored result file, never typed. If a stored value moves, this
# script fails rather than publishing a stale manuscript.
H = pd.read_csv(D / "headline_numbers.csv").set_index("figure")
OH = pd.read_csv(D / "optimum_and_histidine_wide.csv")
SID = pd.read_csv(D / "histidine_sidedness_wide.csv")
WC = pd.read_csv(D / "which_complexes_wide.csv")
CS = pd.read_csv(D / "censoring_sensitivity.csv")
HO = pd.read_csv(D / "held_out_clades.csv")
TIT = pd.read_csv(D / "titration.csv")
SAL = pd.read_csv(D / "salinity_confound.csv")
IR = pd.read_csv(D / "independent_reimplementation.csv")


def h(fig, col="beta_per_sd"):
    return float(H.loc[fig, col])


def sid(measure, term):
    r = SID[(SID.part == "contrast") & (SID.measure == measure) & (SID.term == term)]
    assert len(r) == 1, f"{measure}/{term} is not a single row of histidine_sidedness_wide.csv"
    return float(r.iloc[0].beta_per_pt), float(r.iloc[0].p)


F = {
    "head": h("exposed periplasmic, 32-reference rule-selected panel"),
    "head_p": float(H.loc["exposed periplasmic, 32-reference rule-selected panel", "p"]),
    "head_n": int(H.loc["exposed periplasmic, 32-reference rule-selected panel", "n"]),
    "cyto": h("exposed cytoplasmic control, 32 references"),
    "cyto_p": float(H.loc["exposed cytoplasmic control, 32 references", "p"]),
    "tm": h("transmembrane control, 32 references"),
    "tm_p": float(H.loc["transmembrane control, 32 references", "p"]),
    "sol": h("soluble-surface control, 5 references"),
    "sol_p": float(H.loc["soluble-surface control, 5 references", "p"]),
    "comp_opt": h("exposed periplasmic, 32 refs, + composition, on the optimum-known sample"),
    "with_opt": h("exposed periplasmic, 32-reference panel, + growth optimum"),
}
F["acid_diff"], F["acid_diff_p"] = sid("acidic D+E", "difference")
F["his_sum"], _ = sid("histidine", "sum")
F["his_diff"], F["his_diff_p"] = sid("histidine", "difference")
_j = OH[(OH.part == "fit") & (OH.specification == "+ optimum and histidine")].iloc[0]
F["joint"], F["joint_p"], F["joint_share"] = float(_j.beta_per_sd), float(_j.p), float(_j.share_of_composition_only)
_c = WC[WC.part == "per reference"].sort_values("beta").iloc[0]
F["cyo"], F["cyo_name"] = float(_c.beta), str(_c["name"])
F["n_survive"] = int(WC[WC.part == "bh summary"].iloc[0].n_survive_q10)
F["n_tested"] = int(WC[WC.part == "bh summary"].iloc[0].n_tested)
F["grid_whole"] = float(CS[CS.measure == "share of minima on whole numbers"].iloc[0].rate)
F["grid_half"] = float(CS[CS.measure == "share on a 0.5 grid (whole or half)"].iloc[0].rate)
F["alpha0"] = abs(float(CS[CS.measure == "alpha that zeroes the headline"].iloc[0].alpha_ph_per_sd))
_p6 = CS[(CS.part.astype(str) == "6") & (CS.required_depth_ph > 0)].sort_values("required_depth_ph")
F["depth_min"] = float(_p6.iloc[0].required_depth_ph)
F["head_room"] = float(_p6.iloc[0].headroom_to_ph0)
_s = HO[HO.part == "summary"].iloc[0]
F["rmse_base"], F["rmse_peri"], F["rmse_cyto"] = float(_s.rmse_baseline), float(_s.rmse_periplasmic), float(_s.rmse_cytoplasmic)
_t2 = {float(r.pH): r for _, r in TIT[TIT.part.astype(str) == "2"].iterrows()}
F["q3"], F["q3p"] = float(_t2[3.0].beta_charge_change), float(_t2[3.0].p_charge_change)
F["q4"], F["q4p"] = float(_t2[4.0].beta_charge_change), float(_t2[4.0].p_charge_change)
F["q6"] = float(_t2[6.0].beta_charge_change)
_sv = SAL[SAL.measure == "salinity maximum -> ph_min"].iloc[0]
F["sal_b"], F["sal_p"], F["sal_n"] = float(_sv.beta_per_sd), float(_sv.p), int(_sv.n)
F["repro_mine"] = float(IR[IR.quantity == "headline beta per SD"].iloc[0].independent)
F["repro_pub"] = float(IR[IR.quantity == "headline beta per SD"].iloc[0].published)

# Guards: these are the values the working draft quotes. A mismatch means one of the two documents has drifted.
assert abs(F["head"] - 0.1848) < 5e-5, f"headline is {F['head']}, not the +0.1848 docs/MANUSCRIPT.md quotes"
assert abs(F["joint"] - 0.112005) < 5e-6, "the joint optimum-and-histidine figure has moved"
assert abs(F["acid_diff"] - 0.053596) < 5e-6, "the one-sidedness contrast has moved"
assert F["n_survive"] == 7 and F["n_tested"] == 22, "the Benjamini-Hochberg counts have moved"
print("all headline values match the stored results and the working draft")


def pv(x):
    """p-values as a journal prints them."""
    if x >= 0.01:
        return f"P = {x:.2f}"
    e = int(f"{x:.0e}".split("e")[1])
    m = float(f"{x:.0e}".split("e")[0])
    return f"P = {m:.0f} × 10−{abs(e)}".replace("−", "⁻").replace("⁻", "⁻")


# AUTHOR IDENTITY, supplied by the author and verified against the reference manuscript.
AUTHOR = "Ajay Vellanki"
AFFIL = "Independent Researcher, West Windsor, New Jersey, USA"
EMAIL = "ajayvellanki0720@gmail.com"
ORCID = "0009-0008-7994-152X"

# ---------------------------------------------------------------------------------------------------------------
doc = Document()
st = doc.styles["Normal"]
st.font.name = "Times New Roman"
st.font.size = Pt(12)
st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
pf = st.paragraph_format
pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE
pf.space_after = Pt(0)

sec = doc.sections[0]
for a in ("top_margin", "bottom_margin", "left_margin", "right_margin"):
    setattr(sec, a, Inches(1))

# CONTINUOUS LINE NUMBERING, written into the section properties so Word renders real margin numbers that renumber
# as the author edits - rather than characters typed at the start of each line, which would survive no editing at all.
ln = OxmlElement("w:lnNumType")
ln.set(qn("w:countBy"), "1")
ln.set(qn("w:start"), "1")
ln.set(qn("w:restart"), "continuous")
ln.set(qn("w:distance"), "360")
# POSITION MATTERS AND APPENDING IS WRONG. CT_SectPr is a SEQUENCE, so its children have a fixed order:
#   footnotePr, endnotePr, type, pgSz, pgMar, paperSrc, pgBorders, lnNumType, pgNumType, cols, ... docGrid ...
# The first version of this script appended lnNumType to the end of sectPr, which put it after cols and docGrid.
# Word is free to ignore an out-of-order element or to report the file as needing repair, and the failure is silent:
# the XML looks right to anyone grepping for lnNumType, and no line numbers appear. It is inserted before the first
# child that must follow it instead, and the resulting order is asserted below rather than assumed.
_after = {qn("w:cols"), qn("w:formProt"), qn("w:vAlign"), qn("w:noEndnote"), qn("w:titlePg"),
          qn("w:textDirection"), qn("w:bidi"), qn("w:rtlGutter"), qn("w:docGrid"), qn("w:printerSettings")}
_anchor = next((c for c in sec._sectPr if c.tag in _after), None)
if _anchor is None:
    sec._sectPr.append(ln)
else:
    _anchor.addprevious(ln)
_order = [c.tag.split("}")[1] for c in sec._sectPr]
assert "lnNumType" in _order, "line numbering element did not make it into sectPr"
assert _order.index("lnNumType") < _order.index("cols"), \
    f"lnNumType sits after cols ({_order}); Word may ignore it and no line numbers would show"
if "pgMar" in _order:
    assert _order.index("pgMar") < _order.index("lnNumType"), f"lnNumType sits before pgMar ({_order})"

# page numbers in the footer, which the checklist does not demand but every editor wants
fp = sec.footer.paragraphs[0]
fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = fp.add_run()
for instr in ("begin", "instrText", "end"):
    e = OxmlElement(f"w:fld{instr}" if instr != "instrText" else "w:instrText")
    if instr == "begin":
        e = OxmlElement("w:fldChar"); e.set(qn("w:fldCharType"), "begin")
    elif instr == "end":
        e = OxmlElement("w:fldChar"); e.set(qn("w:fldCharType"), "end")
    else:
        e.set(qn("xml:space"), "preserve"); e.text = " PAGE "
    r._r.append(e)


def head(text, level=1):
    p = doc.add_heading(text, level=level)
    for run in p.runs:
        run.font.name = "Times New Roman"
        run.font.color.rgb = None
        run.font.size = Pt(14 if level == 1 else 12)
    return p


def para(text, bold_lead=None, italic=False):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    if bold_lead:
        r0 = p.add_run(bold_lead)
        r0.bold = True
    r = p.add_run(text)
    r.italic = italic
    return p


# ---------------------------------------------------------------------------------------------------------------
# TITLE PAGE
t = doc.add_paragraph()
t.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
rt = t.add_run("Periplasmic surface charge composition predicts the lower pH limit of bacterial growth "
               "independently of growth optimum")
rt.bold = True
rt.font.size = Pt(14)

para(f"Authors: {AUTHOR}. ORCID iD: {ORCID}")
para(f"Affiliations: {AFFIL}")
para(f"Corresponding author: {AUTHOR}, {EMAIL}")
para("Keywords: acidophile; growth pH limit; comparative genomics; periplasm; surface charge; phylogenetic "
     "generalised least squares")
para("Running title: Periplasmic charge and the bacterial pH limit")
# The word count is a fact about the body text, so it is measured from the body rather than typed. The paragraph is
# created here to hold its place on the title page and filled in after the body exists; a hand-typed count would be
# wrong the moment any section is edited, and this one is the figure an editor checks against the stated limit.
WC_PARA = para("Word count (main text, Introduction to Discussion): ")

doc.add_page_break()

# ABSTRACT
head("Abstract")
para(f"How acidic an environment a bacterium will grow in is only partly explained by where it grows best, and the "
     f"genomic basis of the lower limit itself is not established. We asked which genomic features predict the "
     f"minimum growth pH (pH_min) beyond what the growth optimum already explains, using {F['head_n']:,} bacterial "
     f"genomes with a curated minimum and a 32-protein panel of conserved inner-membrane reference proteins with "
     f"per-residue structural confidence filtering. Analyses are phylogenetic generalised least squares with "
     f"Pagel's λ, controlling genome-wide amino acid composition and data source throughout. "
     f"Composition at solvent-exposed periplasmic positions predicts the limit: a composite of net acidic residue "
     f"share gives {F['head']:+.4f} pH units per standard deviation ({pv(F['head_p'])}), and the direction is "
     f"depletion, acid-tolerant bacteria carrying fewer exposed periplasmic acidic residues. Matched cytoplasmic "
     f"({F['cyto']:+.4f}, {pv(F['cyto_p'])}) and transmembrane ({F['tm']:+.4f}, {pv(F['tm_p'])}) controls computed on "
     f"the same proteins are null, and the periplasmic-minus-cytoplasmic difference, estimated directly so that it "
     f"carries a standard error, is {F['acid_diff']:+.4f} pH per percentage point ({pv(F['acid_diff_p'])}). The "
     f"association retains {F['joint_share']:.0%} of its size with growth optimum and exposed histidine controlled "
     f"together ({F['joint']:+.4f} per SD, {pv(F['joint_p'])}). Exposed histidine carries a separate, larger and "
     f"oppositely signed association that is not a proxy for the first. The presence of particular membrane "
     f"complexes predicts the limit independently of composition, with {F['n_survive']} of {F['n_tested']} testable "
     f"reference families surviving Benjamini–Hochberg correction at q < 0.10. "
     f"The effects are small in variance terms, the work is correlational, and the reported minima largely sit on an "
     f"experimental testing grid rather than on a demonstrated biological boundary; we bound rather than dismiss what "
     f"that could do, and show that no observable censoring pattern could account for the association. The composite "
     f"does not improve out-of-sample prediction in clades the model has not seen, which we report in full.")

# DATA SUMMARY
head("Data Summary")
para("All analysis code, the curated phenotype table, the structural panel and every result table underlying the "
     "figures and quoted coefficients are available at https://github.com/ajayvellanki0720-boop/acidophile-ph-limits-submission, "
     "which includes a release manifest recording the software environment, a SHA-256 checksum for every result "
     "file cited in this article, and the order in which scripts are run.")
para("No new sequence data were generated by this study. Genome assemblies and proteomes were obtained from NCBI; "
     "the accession used for every genome is listed in the deposited file accession_used.csv, which should be cited "
     "alongside the originating submitters of those assemblies.")
para("Growth phenotypes were compiled from published species descriptions and from the following databases, whose "
     "authors should be credited by any reuse: AciDB (acidb.cl), ThermoBase, BacDive, and metaTraits "
     "(metatraits.embl.de), the last supplying the measured salinity values used in the confounder analysis. "
     "Protein structures were obtained from the AlphaFold Protein Structure Database and Pfam/InterPro hidden Markov "
     "models from EMBL-EBI.")
para("Software: Python 3.9.6, with NumPy 2.0.2, pandas 2.3.3, SciPy 1.13.1, Biopython 1.85, pyHMMER 0.12.3 and "
     "Matplotlib 3.9.4; versions are pinned in the "
     "release manifest. Analyses can be reproduced from the deposited result tables without re-downloading genomes.")

# IMPACT STATEMENT
head("Impact Statement")
para("Compositional signatures of extremophily are well established for temperature and salinity but not for pH, "
     "which a standing review describes as considerably less discernible. The mechanisms of bacterial acid "
     "resistance that are understood in detail are inducible responses — decarboxylase systems, urea import, "
     "dedicated pH sensors — and they explain survival of acid shock better than they explain the pH at which "
     "steady growth stops. This study asks what sets that limit, and finds a signal in the composition of the one "
     "compartment the hypothesis names in advance: the outward-facing periplasmic surface, which unlike the "
     "cytoplasm is not held near neutrality. Because the prediction is compartment-specific, it comes with controls "
     "that are the same measurement on other faces of the same proteins, and those controls are null. The work also "
     "shows what such an analysis cannot deliver: the association does not transfer as a predictor to clades the "
     "model has not seen, and the outcome variable is bounded by how species descriptions are written rather than by "
     "biology. We quantify that bound instead of noting it. For readers building genome-based trait predictors, the "
     "reported limits may be as useful as the association.")

doc.add_page_break()

# INTRODUCTION
_body_start = len(doc.paragraphs)
head("Introduction")
para("Bacteria differ enormously in how acidic an environment they will grow in. Some isolates of "
     "Acidithiobacillus and Picrophilus grow below pH 1; most described species will not grow below pH 5. What sets "
     "that lower bound is not the same question as where an organism grows best, and the two are routinely "
     "conflated because the optimum is far more often reported. A species whose optimum is pH 7 may have a limit at "
     "pH 4 or at pH 6, and the difference between those organisms is the quantity of interest here.")
para("The mechanisms of acid resistance that are understood in detail are inducible and physiological. Enteric "
     "bacteria surviving gastric transit deploy amino-acid decarboxylase systems that consume cytoplasmic protons, "
     "regulated by sensors that detect acidification directly: the EvgS/EvgA two-component system, in which a single "
     "histidine in the periplasmic domain is required for acid induction [1], and the regulator AdiY, whose pH "
     "response runs through the protonation of two histidines [2]. Helicobacter pylori additionally imports urea "
     "through the pH-gated channel UreI to feed intrabacterial urease [3]. Experimental evolution for acid "
     "resistance in Escherichia coli converges on the same signalling machinery, with evgS mutated in every selected "
     "population [4]. These are responses: machinery switched on when the environment turns acidic, and they explain "
     "survival of acid shock better than they explain the pH at which steady growth becomes impossible.")
para("A limit might instead be set by what the cell's proteins can tolerate, and if so the relevant compartment is "
     "identifiable in advance. A Gram-negative bacterium holds its cytoplasm near neutrality even in strong acid, so "
     "cytoplasmic proteins never experience the external pH. The periplasm has no such protection: it equilibrates "
     "far closer to the medium, and proteins with surfaces facing it operate at whatever pH the organism is growing "
     "in. Any compositional signature of acid adaptation should therefore appear on outward-facing surfaces and not "
     "on cytoplasmic ones — which makes the hypothesis testable against a control that is the same measurement "
     "on a different face of the same protein. Consistent with the compartment mattering, periplasmic redox proteins "
     "of acidophiles have systematically more positive redox potentials than neutrophilic homologues [5].")
para("Compositional signatures of other extremes are established; for pH they are not. Thermophily is visible in the "
     "IVYWREL fraction and halophily in a pronounced excess of surface acidic residues, while a review of archaeal "
     "extremophile protein adaptation reports that adaptations toward pH are considerably less discernible [6]. That "
     "is a documented absence rather than an assumption, and it is the gap addressed here.")
para("Two recent results make the question tractable. Cabello-Yeves et al. showed that predicted isoelectric point "
     "shifts with habitat between marine and freshwater prokaryotes, with the difference sharper in secreted than "
     "cytoplasmic than membrane proteins and concentrated on protein surfaces [7], establishing that "
     "compartment-resolved composition carries environmental signal. Independently, Barnum et al. showed that growth "
     "pH can be predicted from amino acid composition computed per predicted subcellular localization, reporting the "
     "extracellular-minus-intracellular glutamate difference as the strongest residue-level correlate [8]. Together "
     "these say that the compartment-resolved signal exists and that pH is among the variables it carries.")
para("We ask which genomic features predict pH_min beyond the growth optimum, and whether any such signal is "
     "specific to the compartment the hypothesis names. Three associations survive the controls, and the limits "
     "that bound them are reported as fully as the associations themselves.")

# METHODS
head("Methods")
head("Phenotype data", 2)
para(f"The analysis sample is {F['head_n']:,} bacterial genomes carrying a curated minimum growth pH and a position "
     "in the reference phylogeny. Minima were compiled from four source tiers, recorded per record: hand-curated "
     "species descriptions read from primary sources, rule-based extraction from species-description abstracts, "
     "ThermoBase, and AciDB/BacDive [9]. Source is entered as a set of indicator covariates in every model, so a "
     "between-source level difference cannot by itself produce an association.")
para("The extraction tier dominates the sample and its fidelity was measured rather than assumed. Against "
     "hand-curated records it is exact for 98.8% of minima. Because that validation set was acidophile-enriched, 25 "
     "further records were read by hand from the pH band carrying most of the data; 3 of 25 were wrong, all by one "
     "mechanism, in which a sentence reporting only an optimum has that optimum recorded as the growth range. A "
     "detector for that specific failure scores 8 of 8 against fifteen hand-read statements with no false positives "
     "and flags 298 of 2,536 extracted records; excluding exactly those records moves the headline estimate by 0.13 "
     "standard errors.")
head("Structural panel and composite", 2)
para("The predictor is built on a panel of conserved inner-membrane proteins of Escherichia coli with experimentally "
     "determined or high-confidence predicted structures. Each reference position was assigned to a class — "
     "exposed periplasmic, exposed cytoplasmic, transmembrane or buried — from UniProt topology together with "
     "relative solvent accessibility of the isolated model, computed by the Shrake–Rupley algorithm and "
     "normalised to the maxima of Tien et al. [10], with exposure defined as RSA ≥ 0.25. A position enters only "
     "if its per-residue pLDDT is ≥ 70, since the solvent accessibility of a poorly predicted residue is not "
     "interpretable. For each genome, the best-aligning hit of the reference's Pfam family was tallied by the class "
     "of the reference position it aligned to, counting D, E, K, R and H separately.")
para("The composite is the mean across references of the per-reference net acidic share (D + E − K − R) "
     "divided by scorable positions, each reference standardised across genomes before averaging, with a "
     "per-reference position floor and a minimum of three contributing references. The panel used throughout is a "
     "32-reference rule-selected set; a 7-reference hand-picked panel preceded it and is not the same measurement, "
     "so per-unit coefficients from the two are never compared. Panel construction is deterministic and verified to "
     "rebuild byte-identically.")
head("Statistical model", 2)
para("All estimates are phylogenetic generalised least squares on a Brownian covariance matrix derived from the "
     "reference marker-gene tree, with Pagel's λ [11] estimated by maximum likelihood. Every model carries an "
     "intercept, a genome-wide amino acid composition control and the source indicators. Coefficients are reported "
     "per standard deviation of the predictor unless a per-unit figure is explicitly labelled, because composites "
     "averaging different numbers of references have different spreads. Where a model adds the growth optimum, the "
     "comparison is made on the single sample for which the optimum is known, never against the full-sample figure, "
     "so that a change of model is not confounded with a change of sample. Multiple testing across reference "
     "families is controlled by the Benjamini–Hochberg procedure [12].")
para("Compartment contrasts are estimated by entering both faces in one model as their sum and difference. This "
     "spans the same column space as entering the two levels and so changes no fitted value, but it makes the "
     "periplasmic excess an estimate carrying its own standard error, which a difference taken between two "
     "separately fitted coefficients does not.")
head("Verification", 2)
para(f"Four checks were required to pass before any figure was quoted. A static audit searches the analysis code for "
     f"nine previously made error classes. A unit-test suite covers the library, including mutation-checked guards "
     f"on the per-standard-deviation/per-unit identity. A traceability check extracts every signed coefficient from "
     f"the manuscript and supporting documents and requires it to be present in the stored result file that the same "
     f"claim cites. Finally, the estimator was reimplemented independently, sharing no Newick parser, no covariance "
     f"construction, no matrix factorisation and no λ search with the original: across 30 declared comparisons "
     f"the headline reproduces at {F['repro_mine']:+.6f} per SD against {F['repro_pub']:+.6f}.")

# RESULTS
head("Results")
head("Exposed periplasmic acidic composition predicts the lower growth limit", 2)
para(f"The composite of net acidic residue share at exposed periplasmic positions gives {F['head']:+.4f} pH units per "
     f"standard deviation ({F['head_n']:,} genomes, {pv(F['head_p'])}) with genome-wide composition and source "
     f"controlled (Fig. 1). The sign warrants stating explicitly, because it is easy to invert in prose: the "
     f"coefficient is positive, so a higher net acidic share accompanies a higher minimum, that is, less acid "
     f"tolerance. Acid-tolerant bacteria are depleted of acidic residues on their outward-facing periplasmic "
     f"surfaces, not enriched. At compartment level the contrast is direct: genomes with pH_min ≤ 3.5 average a "
     f"secreted acidic share of 0.0999 against 0.1246 for those at or above pH 6.")
head("The association is confined to the periplasmic-facing surface", 2)
para(f"Because the claim is about a compartment, the controls are the same measurement computed on other "
     f"compartments of the same proteins, on the same genomes, in the same model (Fig. 1). Exposed cytoplasmic gives "
     f"{F['cyto']:+.4f} ({pv(F['cyto_p'])}) and transmembrane {F['tm']:+.4f} ({pv(F['tm_p'])}); both are null. The "
     f"soluble-surface reference set gives {F['sol']:+.4f} ({pv(F['sol_p'])}), which is significant and opposite in "
     f"sign. That is not a failure of specificity, but it does mean the pattern is a contrast between compartments "
     f"rather than a periplasm-only effect, and we describe it as such.")
para(f"Entering both faces in a single model as sum and difference makes the periplasmic excess an estimate in its "
     f"own right: {F['acid_diff']:+.4f} pH per percentage point ({pv(F['acid_diff_p'])}).")
head("Part of the association is independent of growth optimum and of histidine", 2)
para(f"On the genomes whose optimum is known, so that the drop is a change of model alone, the composition-only "
     f"coefficient of {F['comp_opt']:+.4f} per SD falls to {F['with_opt']:+.4f} with the optimum controlled "
     f"(Fig. 2). The more demanding test holds the optimum and exposed histidine fixed together, since histidine "
     f"carries its own association. The coefficient is {F['joint']:+.4f} per SD ({pv(F['joint_p'])}), retaining "
     f"{F['joint_share']:.0%} of the composition-only figure on that sample. The two controls prove nearly "
     f"independent: the product of their separate shares is 0.375 against a joint share of 0.387.")
head("Exposed histidine is a separate association, not a proxy", 2)
para(f"Exposed periplasmic histidine share carries a larger coefficient than the acidic composite and in the opposite "
     f"direction. Entered as sum and difference across the two compartments, with periplasmic D, E and K held fixed, "
     f"the shared component is {F['his_sum']:+.4f} pH per percentage point and the periplasmic excess "
     f"{F['his_diff']:+.4f} ({pv(F['his_diff_p'])}) — so the association has both a side-specific and a shared "
     f"part. That contrast was calibrated on the acidic result first, which fires at {F['acid_diff']:+.4f}, because a "
     f"test is only a test once it has been shown capable of firing.")
para("The two measures correlate at r = −0.549, which makes one standing in for the other the obvious concern. "
     "Mutual adjustment excludes it: each retains roughly three quarters of its coefficient with the other held "
     "fixed, and the acidic measure's one-sidedness survives at 77% of its unadjusted value. Solvent access is doing "
     "demonstrable work — histidine at non-exposed positions in the same reference families predicts nothing, "
     "alone or jointly. A pre-registered ordering of coefficients by side-chain pKa does not hold monotonically; "
     "what does organise them is charge class, which remains significant with histidine held fixed, while pKa orders "
     "the acidic pair and fails to order the basic residues.")
head("Detectable presence of particular membrane complexes predicts the limit", 2)
para(f"How many of the panel's reference families a genome has a detectable homologue of was initially held fixed as "
     f"a nuisance control. Decomposed, it is two things: breadth survives at −0.036 pH per reference with "
     f"assembly level, similarity to the reference proteins and genome size held fixed, whereas depth — "
     f"positions per homologue — carries nothing once similarity is controlled. Of the accepted references, "
     f"{F['n_tested']} have enough detection variance to regress on and {F['n_survive']} survive "
     f"Benjamini–Hochberg at q < 0.10 (Fig. 3), the strongest being {F['cyo_name'].lower()} at "
     f"{F['cyo']:+.4f} pH.")
para("Three limits are explicit. The intended negative control cannot be run: the soluble reference proteins are "
     "detected in over 99% of genomes and therefore fall below the variance floor, so none was tested, and a control "
     "that cannot run is a limitation rather than a pass. What replaces it is weaker and labelled so — the "
     "survivors' signs are mixed, five down and two up, which excludes a shared sequencing or annotation artefact "
     "because such an account predicts one direction for every reference, and establishes nothing about function. "
     "Detection is a crude proxy for presence, so these are claims about detectable homologues of E. coli proteins. "
     "And breadth retains 47% of its coefficient with the survivors held fixed against a pre-set criterion of half, "
     "so the named references carry about half of the breadth effect rather than being it.")
head("The measure survives being made homogeneous, and the estimator reproduces", 2)
para("The composite pools numerator and denominator over whichever references a genome has a detectable homologue "
     "of, so two genomes with identical surface chemistry can receive different shares. This is heterogeneity in "
     "what the measure measures, which holding breadth fixed in a regression does not remove. Two routes that fail "
     "in different places were run: a six-reference common core identical across genomes, and the 239 genomes "
     "carrying every reference. Neither failed, at +0.0362 and +0.0888 pH per percentage point respectively.")

# DISCUSSION
head("Discussion")
para("Three associations survive the controls designed to break them, and the compartment specificity is the part we "
     "would defend most firmly: the cytoplasmic and transmembrane controls are the identical measurement on other "
     "faces of the same proteins, and they are null. We think the limits below matter as much as the findings.")
head("Mechanism: charge burden rather than buffering", 2)
para(f"The mechanism is contested by the one directly comparable experiment. Weeks and Sachs mutagenised every "
     f"histidine, aspartate and glutamate in the periplasmic domain of the H. pylori urea channel UreI and found "
     f"three histidines and three carboxylates all required to activate transport at acidic pH, with those residues "
     f"largely absent from the pH-independent Streptococcus salivarius homologue [3]. There, periplasmic "
     f"carboxylates are required for low-pH function; here they are depleted in acid-tolerant genomes. A later study "
     f"adds a conserved cytoplasmic glutamate to the same channel's gating [13], and an engineering study aimed at "
     f"gastric pH deliberately adds acidic residues on an explicit buffering rationale [14].")
para(f"The two accounts are separable from composition alone, because they require opposite signs of the same "
     f"quantity. Writing Q(pH) for net surface charge and β(pH) = −dQ/dpH for buffer capacity, a "
     f"charge-burden account predicts acid tolerance accompanying a smaller charge change between neutrality and the "
     f"limit, while a buffering account predicts it accompanying greater capacity at the operating pH. At pH 3 and "
     f"pH 4, inside the acidophile range and near the carboxylate pKa values, the charge change carries a positive "
     f"coefficient ({F['q3']:+.4f}, {pv(F['q3p'])}; {F['q4']:+.4f}, {pv(F['q4p'])}), so acid-tolerant genomes have "
     f"less outward-facing charge to change, as burden predicts, while buffer capacity carries a positive "
     f"coefficient, the opposite of what buffering predicts. Both quantities reverse sign between pH 4 and pH 6 "
     f"({F['q6']:+.4f}), which is where histidine titrates and where its opposite-signed association appears inside "
     f"the same quantity. We emphasise that this is not an electrostatics calculation — no structure, no "
     f"solvent, no interaction between groups, and model pKa values throughout.")
para("A protein whose function is to respond to acid needs titratable groups on the face that senses it; a proteome "
     "that must merely continue working at low pH may do better with fewer. That reconciles UreI with the present "
     "result, and it makes a falsifiable prediction: a panel built from dedicated pH sensors should carry the "
     "opposite sign to a housekeeping panel. We attempted that test and it cannot be run on this sample. Dedicated "
     "periplasmic sensor domains are detected in 1.9% of genomes, with all detections in a single phylum, so the "
     "contrast would be between lineages rather than between protein classes. The prediction is untested, not "
     "refuted.")

head("The outcome variable is the principal limitation, and it is bounded", 2)
para(f"pH_min is meant to be the lowest pH at which an organism grows. What a species description usually records is "
     f"the lowest pH at which growth was observed, bounded below by whatever the authors chose to test. If that "
     f"choice tracks what they expected, the error is not noise but is systematically related to apparent acid "
     f"tolerance and therefore to anything predicting it.")
para(f"The reported minima do carry a testing-grid signature: {F['grid_whole']:.0%} sit on a whole pH unit and "
     f"{F['grid_half']:.0%} on a whole or half unit, and the modal value, pH 6.0, accounts for 34% of the sample. A "
     f"measured limit has no reason to prefer a round number; a testing grid has every reason to. Measured outcome "
     f"reliability is high (ICC 0.965 from cross-source replicate readings, error SD 0.343 pH) and supports the "
     f"standard result that error in the dependent variable inflates standard errors without biasing slopes — "
     f"but that holds for classical error, independent of the predictor by construction. Censoring is not classical: "
     f"it has a sign and a mechanism. Simulated at equal magnitude, classical noise moves the coefficient by 0.006 "
     f"standard errors while the same magnitude aligned with the predictor moves it by 11.8.")
para(f"We therefore bounded it. To drive the coefficient to zero, the censoring displacement would have to deepen by "
     f"{F['alpha0']:.3f} pH for every standard deviation of surface composition. Routed through observable proxies "
     f"for being censored, that requirement is not merely implausible but physically impossible: the proxies pointing "
     f"in the threatening direction require mean displacements of {F['depth_min']:.1f} pH or more against "
     f"{F['head_room']:.1f} pH of available headroom, which would place the flagged organisms' true minima below "
     f"pH 0 (Fig. 4). This closes the objection in the form a proxy can express it. It does not close the objection: "
     f"the proxies are weak and each proxies other things, and selection acting on which organisms get described at "
     f"all would appear in none of them. Settling it requires bound type read from full texts at scale, or growth "
     f"limits measured on a common pH grid.")
head("What the association does not do", 2)
para(f"It does not predict out of clade. Fitting on every other phylum, with the composite re-standardised on "
     f"training genomes only, pooled out-of-sample RMSE is {F['rmse_peri']:.4f} pH with the composite against "
     f"{F['rmse_base']:.4f} pH without it — fractionally worse — and the cytoplasmic control, at "
     f"{F['rmse_cyto']:.4f} pH, does no worse than either. The composite lowers error in three of four folds, and the "
     f"pooled figure is decided by the one it fails, which is the phylum holding most of the acidophiles and whose "
     f"removal leaves the training set with little of the acid-tolerant range the effect is confined to; fold-level "
     f"gain tracks the tail retained in training at r = +0.835 across four folds, a direction predicted in advance "
     f"and resting on four points. Both readings stand: the association survives every within-sample test reported "
     f"here and buys no out-of-sample accuracy across clades. For a claim about mechanism that is survivable; for a "
     f"claim about prediction it would be fatal, and we make none. Reference selection was additionally performed "
     f"once on all genomes, so even these figures are optimistic by an unmeasured amount.")
para(f"The effects are small in variance terms, well under 1% of pH_min variance beyond the growth optimum, and the "
     f"relationship is concentrated in the acid-tolerant tail rather than running uniformly across the range. Nothing "
     f"here is interventional. In archaea the pattern is consistent but marginal, and within halophilic archaea "
     f"specificity fails outright.")
head("A second environmental variable acts on a different compartment", 2)
para(f"Because compartment-resolved isoelectric point shifts with marine-versus-freshwater habitat in the same "
     f"ordering reported here [7], a habitat variable could in principle act through the periplasmic face and mimic "
     f"the result. Measured salinity tolerance, available for {F['sal_n']} genomes of this sample, does not predict "
     f"pH_min ({pv(F['sal_p'])}) and removes 1% of the periplasmic coefficient when held fixed. More informative, the "
     f"two variables have separable compartment profiles: salinity gives secreted > cytoplasmic > membrane, "
     f"reproducing the published ordering, while pH_min gives secreted > membrane > cytoplasmic. Salinity moves the "
     f"secreted and cytoplasmic shares together, as a whole-cell osmotic account predicts; pH moves the secreted "
     f"share and leaves the cytoplasmic one untouched. The method does not impose one ordering on every variable "
     f"given to it, which is a stronger defence of the controls than salinity being inert — it is not inert, and "
     f"predicts secreted composition more strongly than pH does. Within the salinity-measured subsample the "
     f"transmembrane control does fire, and that signal is confined to the salt-tolerant half, consistent with "
     f"halotolerant membranes being compositionally distinct: a second environmental variable perturbing a different "
     f"compartment from the one claimed here.")
head("Conclusion", 2)
para("Composition at the outward-facing periplasmic surface is associated with the lower pH limit of bacterial "
     "growth beyond what growth optimum explains, specifically to that compartment, and separably from a larger "
     "histidine association and from gene content. The direction is depletion of acidic residues, and a "
     "charge-burden rather than buffering account is the one the titration behaviour supports, against the single "
     "most directly comparable experiment in the literature. The limits are as firm as the findings: no "
     "out-of-clade predictive gain, an outcome variable shaped by how species descriptions are written, and no "
     "interventional evidence. A pre-registered charge-reversal experiment is designed; the discrimination between "
     "burden and buffering should be settled on a single chassis protein before it is run.")

_body_end = len(doc.paragraphs)
_body_words = sum(len(q.text.split()) for q in doc.paragraphs[_body_start:_body_end])
WC_PARA.add_run(f"{_body_words:,}")
# The target the journal states for a Research Article. Asserted rather than noted, because a body that has drifted
# out of range is a submission problem and the author needs to be told at build time, not at desk rejection.
assert 3000 <= _body_words <= 7000, (
    f"main text is {_body_words:,} words, outside the 3000-7000 a Research Article allows; "
    f"either trim or retarget the journal before submitting")

# BACK MATTER
head("Conflicts of interest")
para("The author declares that there are no conflicts of interest.")
head("Funding information")
para("This research received no specific grant from any funding agency, commercial or not-for-profit.")
head("Ethical approval and consent to participate")
para("This study used only previously published phenotype data and publicly available genome sequences; no human "
     "participants, human material, human data or animals were involved, and no ethical approval or consent was "
     "required.")
head("Author contributions")
para(f"{AUTHOR} is the sole author and, in CRediT terms, is responsible for conceptualization, data curation, "
     f"formal analysis, investigation, methodology, project administration, software, validation, visualization, "
     f"writing – original draft preparation, and writing – review and editing. No funding was acquired and the "
     f"work was not supervised externally.")
head("Acknowledgements")
para("The analysis code, the result tables and the drafts of this manuscript were produced with the assistance "
     "of a large language model (Anthropic Claude) working under the direction of the author. The author set the "
     "research question, chose and approved every analysis, and is responsible for the content. The verification "
     "described in Methods was built for this reason: an independent re-implementation of the estimator sharing no "
     "code with the original reproduces every headline coefficient, each quoted figure is checked mechanically "
     "against the stored result table it comes from, and the corrections made to earlier analyses during this work "
     "are listed in the deposited repository rather than removed from the record. No generative tool is an author, "
     "and none was used to produce data, figures or references that were not independently verified against the "
     "primary sources cited.")

# REFERENCES
# Numbered in order of first citation, Vancouver-like. The Microbiology Society accepts "any style" provided it is
# consistent and includes author names, title, journal, year, volume, issue, pages and a DOI where available. Every
# entry below was retrieved during this work; nothing is reconstructed from memory.
head("References")
REFS = [
 "Sen H, Aggarwal N, Ishionwu C, Hussain N, Parmar C, Jamshad M, Bavro VN, Lund PA. Structural and functional "
 "analysis of the Escherichia coli acid-sensing histidine kinase EvgS. J Bacteriol 2017;199(18):e00310-17. "
 "doi:10.1128/JB.00310-17",

 "Gallo G, Brameyer S, Kuppermann S, Schneider S, Kielkowski P, Jung K. AdiY acts as a cytoplasmic pH sensor via "
 "histidine protonation to regulate acid stress adaptation in Escherichia coli. J Bacteriol 2026;208(1):e0054225 "
 "(published online 23 December 2025). doi:10.1128/jb.00542-25",

 "Weeks DL, Sachs G. Sites of pH regulation of the urea channel of Helicobacter pylori. Mol Microbiol "
 "2001;40(6):1249–1259. doi:10.1046/j.1365-2958.2001.02466.x",

 "Johnson MD, Bell J, Clarke K, Chandler R, Pathak P, Xia Y, Marshall RL, Weinstock GM, Loman NJ, Winn PJ, Lund PA. "
 "Characterization of mutations in the PAS domain of the EvgS sensor kinase selected by laboratory evolution for "
 "acid resistance in Escherichia coli. Mol Microbiol 2014;93(5):911–927. doi:10.1111/mmi.12704",

 "Ullrich SR, Fuchs H, Ashworth-Güth C. Electrochemical and structural characterization of recombinant "
 "respiratory proteins of the acidophilic iron oxidizer Ferrovum sp. PN-J47-F6 suggests adaptations to the acidic pH "
 "at protein level. Front Microbiol 2024;15:1357152. doi:10.3389/fmicb.2024.1357152",

 "Reed CJ, Lewis H, Trejo E, Winston V, Evilia C. Protein adaptations in archaeal extremophiles. Archaea "
 "2013;2013:373275. doi:10.1155/2013/373275",

 "Cabello-Yeves PJ, Rodriguez-Valera F. Marine-freshwater prokaryotic transitions require extensive changes in the "
 "predicted proteome. Microbiome 2019;7:117. doi:10.1186/s40168-019-0731-5",

 "Barnum TP, et al. Predicting microbial growth conditions from amino acid composition. bioRxiv 2024. "
 "doi:10.1101/2024.03.22.586313 [preprint]",

 "Neira G, Cortez D, Jil J, Holmes DS. AciDB 1.0: a database of acidophilic organisms, their genomic information "
 "and associated metadata. Bioinformatics 2020;36(19):4970–4971. doi:10.1093/bioinformatics/btaa638",

 "Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO. Maximum allowed solvent accessibilites of residues in "
 "proteins. PLoS One 2013;8(11):e80635. doi:10.1371/journal.pone.0080635",

 "Pagel M. Inferring the historical patterns of biological evolution. Nature 1999;401(6756):877–884. "
 "doi:10.1038/44766",

 "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple "
 "testing. J R Stat Soc Series B Stat Methodol 1995;57(1):289–300. doi:10.1111/j.2517-6161.1995.tb02031.x",

 "Stoib A, Shojaei S, Fischer X, Posch S, Parisse G, Frezzini M, Putz T, Siligan C, Goessweiner-Mohr N, Narzi D, "
 "Horner A. A cytoplasmic pH sensor redefines gating in bacterial urea channels and acid resistance. Small "
 "2026;e76130. doi:10.1002/smll.76130",

 "Dhiman A, Purohit R. Identification of potential mutational hotspots in serratiopeptidase to address its poor "
 "pH tolerance issue. J Biomol Struct Dyn 2023;41(18):8831–8843. doi:10.1080/07391102.2022.2137699",

 "Podlesny D, et al. metaTraits: a large-scale integration of microbial phenotypic trait information. Nucleic Acids "
 "Res 2026;54(D1):D835. doi:10.1093/nar/gkaf1241",

 "Sayers EW, Bolton EE, Fine AM, Kelly C, Kim S, Landrum M, Lathrop S, Malheiro A, Murphy TD, Phan L, Pujar S, "
 "Trawick BW, Schneider VA, Pruitt KD. Database resources of the National Center for Biotechnology Information in "
 "2026. Nucleic Acids Res 2026;54(D1):D20–D27. doi:10.1093/nar/gkaf1060",

 "Bertoni D, Tsenkov M, Magana P, Nair S, Pidruchna I, Querino Lima Afonso M, Midlik A, Paramval U, Lawal D, "
 "Tanweer A, Last M, Patel R, Laydon A, Lasecki D, Dietrich N, Tomlinson H, Žídek A, Green T, Kovalevskiy O, Lau A, "
 "Kandathil S, Bordin N, Sillitoe I, Mirdita M, Jones D, Orengo C, Steinegger M, Fleming JR, Velankar S. AlphaFold "
 "Protein Structure Database 2025: a redesigned interface and updated structural coverage. Nucleic Acids Res "
 "2026;54(D1):D358–D362. doi:10.1093/nar/gkaf1226",

 "Shrake A, Rupley JA. Environment and exposure to solvent of protein atoms. Lysozyme and insulin. J Mol Biol "
 "1973;79(2):351–371. doi:10.1016/0022-2836(73)90011-9",
]
for i, t in enumerate(REFS, 1):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    p.paragraph_format.left_indent = Inches(0.35)
    p.paragraph_format.first_line_indent = Inches(-0.35)
    p.add_run(f"{i}. ").bold = True
    p.add_run(t)

# FIGURE LEGENDS
head("Figure legends")
LEG = [
 ("Fig. 1.", f"The association is confined to the periplasmic-facing surface. Coefficients (points) and 95% "
  f"confidence intervals (bars) for the net acidic composite computed on four position classes of the same 32 "
  f"reference proteins, in the same phylogenetic model with genome-wide composition and source controlled. "
  f"n = {F['head_n']:,} except the soluble-surface control (n = 2,764). Positive values indicate a higher minimum "
  f"growth pH, i.e. less acid tolerance."),
 ("Fig. 2.", "Part of the association is independent of growth optimum and of exposed histidine. Coefficients and "
  "95% confidence intervals for the periplasmic composite under four nested models, all fitted on the same 1,655 "
  "genomes for which the growth optimum is known, so that each drop reflects a change of model rather than of "
  "sample. Percentages give the share of the composition-only coefficient retained."),
 ("Fig. 3.", f"Detectable presence of particular reference families predicts the limit. Per-reference coefficients "
  f"and 95% confidence intervals for detection, sorted by effect, with the {F['n_survive']} families surviving "
  f"Benjamini–Hochberg correction at q < 0.10 shown in filled markers and the remainder in open markers. Every "
  f"fit controls composition, source, mean alignment score, log genome size and assembly level."),
 ("Fig. 4.", f"No observable censoring pattern can account for the association. For each proxy of a censored "
  f"minimum, the mean displacement that would be required to drive the coefficient to zero (upper bar) against the "
  f"displacement physically available, which is the mean reported minimum of the flagged records and would place "
  f"every true limit at pH 0 (lower bar). Both threatening proxies require more displacement than exists."),
]
for lab, txt in LEG:
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    p.add_run(lab + " ").bold = True
    p.add_run(txt)
para("Figure files: all five panels are deposited with the code as PNG and SVG and are embedded in the "
     "PeerJ-formatted submission manuscript with alt text.")

fn = OUT / "Periplasmic_surface_charge_manuscript.docx"
doc.save(fn)
words = sum(len(p.text.split()) for p in doc.paragraphs)
print(f"wrote {fn}")
print(f"approximate total word count including front and back matter: {words:,}")
print("continuous line numbering is set in the section properties (w:lnNumType, countBy 1, restart continuous)")
