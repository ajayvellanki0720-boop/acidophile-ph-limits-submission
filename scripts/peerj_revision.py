"""Apply submission revisions without changing stored analyses or coefficients."""
from docx.oxml.ns import qn
from docx.shared import Pt, RGBColor
import re
from copy import deepcopy

TITLE = "Extracytoplasmic surface charge composition is associated with reported minimum growth pH after adjustment for growth optimum"


def replace_paragraph(doc, start, text):
    p = next(p for p in doc.paragraphs if p.text.startswith(start))
    p.text = text
    return p


def before(anchor, text, style=None):
    p = anchor.insert_paragraph_before(text, style=style)
    return p


def revise(doc):
    # The submission portal adds its own headers and footers.
    for section in doc.sections:
        for part in (section.header, section.footer, section.first_page_header, section.first_page_footer):
            for paragraph in part.paragraphs:
                paragraph.clear()
    next(p for p in doc.paragraphs if p.text == 'Abstract').paragraph_format.page_break_before = True
    doc.paragraphs[0].text = TITLE
    doc.paragraphs[0].style = doc.styles['Title']
    for r in doc.paragraphs[0].runs:
        r.font.size = Pt(16)
        r.bold = True
    doc.core_properties.title = TITLE
    for p in doc.paragraphs:
        if p.text == 'Methods':
            p.text = 'Materials & Methods'
        if p.text == 'Conclusion':
            p.text = 'Conclusions'
            p.style = doc.styles['Heading 1']
        headings = {
            'Exposed periplasmic acidic composition predicts the lower growth limit': 'Exposed outward-facing composition is associated with reported minimum growth pH',
            'The association is confined to the periplasmic-facing surface': 'Associations differ between membrane position classes',
            'Detectable presence of particular membrane complexes predicts the limit': 'Detectable reference families are associated with reported minimum growth pH',
            'Mechanism: charge burden rather than buffering': 'Charge burden and buffering hypotheses',
            'The outcome variable is the principal limitation, and it is bounded': 'Reported minima and censoring sensitivity',
        }
        if p.text in headings:
            p.text = headings[p.text]
        if p.text.startswith('Word count ('):
            p.text = ''

    abstract = next(p for p in doc.paragraphs if p.text.startswith('How acidic an environment'))
    abstract_parts = [
        ('Background.', 'The genomic correlates of reported minimum growth pH may differ from those of growth optimum. We tested whether charge composition on the outward-facing surfaces of conserved membrane proteins is associated with reported minimum growth pH after adjustment for growth optimum.'),
        ('Methods.', 'We analysed 2,768 bacterial genomes using a 32-reference membrane-protein panel, structural confidence and solvent-accessibility filters, and phylogenetic generalised least squares with genome-wide composition and source controls. Matched cytoplasmic and transmembrane positions provided compartment controls. Models incorporating growth optimum used the same 1,655 genomes for each comparison.'),
        ('Results.', 'The extracytoplasmic net acidic composite was positively associated with reported minimum growth pH (+0.1848 pH units per standard deviation; P = 2 × 10⁻10). The cytoplasmic and transmembrane controls were not statistically significant. A separate pooled residue-share analysis adjusted for growth optimum and exposed histidine retained 39% of its own baseline association (+0.1120 per standard deviation; P = 4 × 10⁻3). Seven of 22 detectable reference families passed Benjamini–Hochberg correction at q < 0.10. The composite did not improve prediction in held-out phyla.'),
        ('Conclusions.', 'Acid-tolerant bacteria have lower net acidic composition on extracytoplasmic-facing positions of this reference panel. The result is an observational association; most reported minima are not experimentally demonstrated boundaries. A charge-burden interpretation is a hypothesis requiring experimental testing.'),
    ]
    for label, text in abstract_parts:
        p = before(abstract, '')
        p.add_run(label + ' ').bold = True
        p.add_run(text)
    abstract._element.getparent().remove(abstract._element)

    replace_paragraph(doc, 'Bacteria differ enormously',
        'Bacteria differ greatly in the acidity at which growth is observed. Some Acidithiobacillus isolates grow below pH 1; the archaeal genus Picrophilus also includes organisms capable of growth at very low pH. A reported minimum growth pH is distinct from the growth optimum. A species with an optimum near pH 7 may have reported growth extending to pH 4 or only to pH 6. Here we study that reported minimum, without assuming it is a demonstrated biological boundary.')

    intro_end = next(p for p in doc.paragraphs if p.text.startswith('We ask which genomic features'))
    intro_end.text = 'We test which genomic features are associated with reported minimum growth pH after adjustment for growth optimum, and whether the association differs between outward-facing and cytoplasmic positions of the same proteins.'
    before(intro_end, 'We use extracytoplasmic to describe the outward-facing side of the cytoplasmic membrane across the bacterial sample. In diderm bacteria these positions face the periplasm; in monoderms they face the cell wall and extracellular environment. Position classes are transferred from reference annotations by sequence alignment, rather than measured structurally in each sampled organism.')

    anchor = next(p for p in doc.paragraphs if p.text == 'Structural panel and composite')
    before(anchor, 'One assembly was retained per species. The bacterial phenotype table contains 2,771 assemblies; intersection with the reference phylogeny produces the 2,768-genome headline sample. Before that intersection, source tiers contain 2,536 automatically extracted records, 100 ThermoBase records, 90 hand-curated records, 30 AciDB records and 15 BacDive records. Sources can overlap; the deposited table records the selected source and assembly provenance. The hand-curated tier is enriched for acidophiles and is not a representative random validation set.')
    before(anchor, 'The 25-record central-range check found three optimum-as-range extraction errors. The targeted detector was assessed on 15 hand-read statements and flagged 298 extracted records. The primary analysis retains the deposited phenotype table; exclusion of flagged records is a sensitivity analysis, not a correction of those records. Only five species were recorded as having a strictly demonstrated lower limit in the project curation. Cross-source agreement largely measures agreement in transcription of shared literature and does not establish the precision of the original growth experiments.')

    anchor = next(p for p in doc.paragraphs if p.text == 'Statistical model')
    before(anchor, 'Reference selection used reviewed E. coli K-12 proteins with at least four transmembrane helices, at least three annotated outward-facing and three cytoplasmic topological domains, and exactly one Pfam family. One candidate per family was retained using the number of annotated topological domains. Screening used 150 sampled bacterial proteomes, an 80% detection threshold, median model pLDDT of at least 70, and at least 15 exposed outward-facing positions after residue-confidence filtering. Accepted reference identifiers and annotations are deposited in data/curated/structure_panel_refs.csv; five conserved soluble references serve as controls.')
    before(anchor, 'For the headline composite, a reference contributes only when the number of scorable positions reaches max(8, half the reference class length). Its net acidic share is (D + E − K − R)/n. Each reference share is centred and scaled across the available genomes, then available standardised shares are averaged; at least three references must contribute. Model rows are restricted to complete outcomes, predictors and covariates; missing growth optima are not imputed. The genome-wide composition control is acidic excess, (D + E − K − R), expressed as percentage points. Source indicators use hand-curated literature as the reference category.')
    before(anchor, 'A separate analysis pools counts across reference proteins, requires at least 200 outward-facing exposed positions, and models net residue share together with exposed histidine. It is a different predictor from the mean of standardised reference shares. Retained percentages are calculated within that pooled-share analysis and must not be computed by dividing its coefficient by a coefficient from the headline composite. The exact model samples, scales, coefficients and standard errors are recorded in headline_numbers.csv and optimum_and_histidine_wide.csv.')

    anchor = next(p for p in doc.paragraphs if p.text == 'Verification')
    before(anchor, 'Phylogeny construction and sensitivity', 'Heading 2')
    before(anchor, 'The phylogeny uses 29 ribosomal-protein Pfam marker families. Highest-scoring gathering-cutoff hits were aligned to their Pfam profiles with pyHMMER hmmalign using match columns and concatenated. Genomes with fewer than half the retained families and columns with less than 50% occupancy were removed. Gamma-corrected Poisson distances on shared non-gap sites were used for neighbour joining, followed by midpoint rooting. The marker list, assembly order and retained tree are deposited with the code. The Brownian covariance is shared root-to-tip path length; Pagel’s λ scales its off-diagonal entries and is estimated by maximum likelihood.')
    before(anchor, 'The headline association remained significant on 25 column-bootstrap trees, with coefficients ranging from +0.1751 to +0.1863 per standard deviation (median +0.1822). These trees share an alignment and inference method. No maximum-likelihood or approximate-maximum-likelihood tree sensitivity analysis was completed for the full 2,768-genome sample; shared systematic topology or branch-length error remains a limitation. Parametric null calibration and bootstrap analyses assess specified assumptions, rather than proving the phylogenetic model correct.')
    before(anchor, 'Panel provenance and analysis scope', 'Heading 2')
    before(anchor, 'The canonical stored panel predates deterministic paralogue tie-breaking. A corrected build sorts tied hits and output rows and reproduces byte-identically on repeated runs, but the canonical results have not all been migrated to that build. Refitting the authoritative headline table on the corrected panel changes +0.1848 to +0.1858 per standard deviation (0.035 standard errors), with no sign changes or p = 0.05 crossings among its 14 rows. This comparison bounds that measured impact; it does not make the original stored panel a byte-identical rebuild. Independent estimator replication begins from stored intermediates and does not independently validate phenotype curation, hit detection or structural annotation.')

    p = next(p for p in doc.paragraphs if 'Panel construction is deterministic' in p.text)
    p.text = p.text.replace('Panel construction is deterministic and verified to rebuild byte-identically.', 'Canonical and deterministic panel provenance is distinguished below; the canonical stored panel is not a byte-identical rebuild of the corrected code.')
    replace_paragraph(doc, 'On the genomes whose optimum is known', 'On the same 1,655 genomes with a recorded optimum, the headline mean-of-standardised-shares composite gives +0.2123 per standard deviation with composition and source controlled and +0.0904 with growth optimum also controlled. Figure 2 presents a separate pooled-share analysis on these genomes: adding growth optimum and exposed histidine together leaves +0.1120 per standard deviation (P = 4 × 10⁻3), or 39% of that pooled predictor’s own baseline coefficient. These coefficients come from different predictor definitions and should not be divided to estimate a retained percentage. The product of the pooled analysis’s separate retained shares is 0.375, compared with a joint share of 0.387.')
    replace_paragraph(doc, 'The composite pools numerator and denominator', 'The contribution of reference families varies between genomes, so holding reference breadth fixed does not guarantee that every genome’s composite measures identical positions. Sensitivity analyses use a six-reference common core and the 239 genomes carrying the complete reference set; the corresponding associations are +0.0362 and +0.0888 pH units per percentage point. These results support stability under those restrictions, without establishing that reference-level measurements are interchangeable.')
    p = next(p for p in doc.paragraphs if p.text.startswith('Three limits are explicit.'))
    p.text = p.text.replace('which excludes a shared sequencing or annotation artefact because such an account predicts one direction for every reference, and establishes nothing about function.', 'which argues against a uniform directional artefact but does not exclude family-specific detection bias or establish function.')
    replace_paragraph(doc, 'Four checks were required', 'Verification used a static audit for nine previously encountered code-error classes, 26 unit tests, numerical traceability checks on the working Markdown manuscript and supporting documents, and an independently implemented PGLS estimator. The traceability check verifies that quoted numbers occur in a cited result table, but does not establish that a number supports its attached interpretation; it does not directly parse this Word document. The independent estimator shares no Newick parser, covariance construction, matrix factorisation or λ search with the original. Across 30 declared comparisons, its headline estimate is +0.184838 per standard deviation compared with +0.184837 in the original. These checks assess specified computational properties rather than independently validating all biological inputs.')
    p = next(p for p in doc.paragraphs if p.text.startswith('We therefore bounded it.'))
    p.text = ('We evaluated specified observable censoring proxies. Erasing the headline slope requires censoring displacement to change by 0.185 pH units per predictor standard deviation (Fig. 4A). For the two proxies with effects in the threatening direction, the required mean displacements are at least 8.4 pH units, exceeding the roughly 6.0-unit headroom under an assumed lower bound of pH 0 (Fig. 4B). This bound applies only to those proxies and that assumed floor; pH 0 is not a universal physical lower limit. The proxies are weak and may capture other traits. Predictor-related censoring not represented by these proxies, or selection of organisms for description, remains unexcluded. Common-grid measurements or full-text classification of bound type would provide a stronger test.')
    p = next(p for p in doc.paragraphs if p.text.startswith('The two accounts are separable'))
    p.text = p.text.replace('The two accounts are separable from composition alone, because they require opposite signs of the same quantity.', 'We contrasted simplified charge-burden and buffering hypotheses using composition-derived titration summaries; these summaries do not identify a mechanism on their own.')
    p = next(p for p in doc.paragraphs if p.text.startswith('Composition at the outward-facing'))
    p.text = ('Extracytoplasmic-facing charge composition in the reference panel is associated with reported minimum growth pH after adjustment for growth optimum and histidine. The direction is lower net acidic composition in acid-tolerant organisms. Composition-derived titration summaries are consistent with a charge-burden hypothesis, but do not establish it over alternative mechanisms. The outcome is shaped by growth-testing conventions, the phylogeny and structural annotation impose additional assumptions, and the composite provides no out-of-clade predictive gain. An experimental charge-reversal protocol has been specified but not performed; causal interpretation awaits experimental evidence.')
    p = next(p for p in doc.paragraphs if p.text.startswith('It does not predict out of clade.'))
    p.text = p.text.replace('For a claim about mechanism that is survivable; for a claim about prediction it would be fatal, and we make none.', 'These results support an in-sample conditional association and do not demonstrate generalisable prediction or a causal mechanism.')

    anchor = next(p for p in doc.paragraphs if p.text == 'What the association does not do')
    before(anchor, 'Reference transfer and envelope architecture', 'Heading 2')
    before(anchor, 'The reference models are isolated chains and do not represent every oligomeric interface or species-specific structure. For alternative reference structures from the same Pfam family, per-genome measurements correlate weakly (median r = 0.126), although rescoring the same reference is reproducible. Topology-only annotation retains approximately 101% of the headline effect; solvent exposure concentrates signal per position but is not required for the aggregate association. The strongest interpretation is an outward-facing versus cytoplasmic contrast, rather than an exposure-specific causal effect.')
    before(anchor, 'The association is present in both diderm and monoderm subsets. In a ten-phylum meta-regression, the monoderm increment is +0.1239 per standard deviation (SE 0.0788; one-sided P = 0.058). Envelope classification is coarse and partly confounded with phylogeny. These comparisons motivate consistent use of extracytoplasmic terminology and do not establish that outer-membrane architecture modifies the effect.')

    p = next(p for p in doc.paragraphs if p.text.startswith('All analysis code, the curated'))
    p.text = ('The code, curated inputs, stored results and figures are available at https://github.com/ajayvellanki0720-boop/acidophile-ph-limits-submission. The code and reproducibility archive (version 1.0.0, commit 97257f22dde358b920c70a67aa921a08556b41fb) is preserved at https://doi.org/10.5281/zenodo.24147308. Supplemental Data S1 is the accompanying analysis-data archive, containing curated phenotype and accession tables, reference annotations, stored model inputs and results, phylogenies, environment specifications, provenance and SHA-256 checksums. It is supplied for direct upload as journal-hosted supplemental material. Oversized third-party raw genomes, proteomes and structures are accessed through the accession and retrieval instructions in DATA_AVAILABILITY.md; no new sequence data were generated. The repository documents which reconstruction steps require those inputs and manual curation.')
    p = next(p for p in doc.paragraphs if p.text.startswith('Software: Python'))
    p.text += ' Software versions describe the original analysis environment; an environment specification and package requirements accompany the deposited archive.'
    p = next(p for p in doc.paragraphs if p.text.startswith('Growth phenotypes were compiled'))
    p.text = p.text.replace('metaTraits (metatraits.embl.de)', 'metaTraits (metatraits.embl.de) [15]').replace('AlphaFold Protein Structure Database', 'AlphaFold Protein Structure Database [17]')
    p = next(p for p in doc.paragraphs if p.text.startswith('No new sequence data were'))
    p.text = p.text.replace('obtained from NCBI;', 'obtained from NCBI [16];')
    p = next(p for p in doc.paragraphs if p.text.startswith('The predictor is built'))
    p.text = p.text.replace('Shrake–Rupley algorithm', 'Shrake–Rupley algorithm [18]')
    p = next(p for p in doc.paragraphs if p.text.startswith('The analysis code, the result tables and the drafts'))
    p.text = ('Anthropic Claude Opus 5, as identified by the author, assisted development of the analysis code and earlier manuscript drafts. OpenAI Codex using GPT-6.1 assisted manuscript revision, reproducibility packaging and figure-export code. These tools were used to support coding, checks and editorial preparation. The scientific figures are programmatic plots of stored numerical results, not generated research imagery. The author is responsible for reviewing and approving the content, code, data visualisations and references. Earlier and revised drafts and the analysis scripts are retained in the project record.')
    p.text += ' The author confirms that the manuscript, references and figures have been checked for originality and accuracy, that publication rights and any required permissions are secured, and that the applicable AI tool terms have been checked for suitability for publication. The author takes full responsibility for the integrity of the complete article, including the references.'

    # Adopt cross-domain terminology and remove unsubstantiated registration claims.
    for p in doc.paragraphs:
        if p.text.startswith('Fig. '):
            continue
        for r in p.runs:
            r.text = r.text.replace('periplasmic', 'extracytoplasmic').replace('Periplasmic', 'Extracytoplasmic')
            r.text = r.text.replace('pre-registered ordering', 'previously specified ordering').replace('pre-registered charge-reversal experiment', 'specified charge-reversal experiment')
    # Journal upload copy: collect legends and remove the embedded figure paragraphs.
    legends = []
    for p in list(doc.paragraphs):
        if p._p.findall('.//' + qn('w:drawing')):
            p._p.getparent().remove(p._p)
        elif p.text.startswith('Fig. '):
            legends.append(p.text)
            p._p.getparent().remove(p._p)
    doc.add_heading('Figure legends', level=1)
    for text in legends[:3]:
        if text.startswith('Fig. 1.'):
            text = text.replace('The association is confined to the periplasmic-facing surface.', 'Associations differ between membrane position classes.')
        if text.startswith('Fig. 2.'):
            text = text.replace('Part of the association is independent of growth optimum and exposed histidine.', 'Pooled residue-share associations after adjustment for growth optimum and exposed histidine.')
        doc.add_paragraph(text.replace('periplasmic', 'extracytoplasmic'))
    doc.add_paragraph('Fig. 4. Sensitivity to specified differential-censoring models. (A) The headline coefficient reaches zero when displacement changes by 0.185 pH units per predictor standard deviation. (B) For the two proxies with effects in the threatening direction, required displacement exceeds available headroom under an assumed pH 0 floor. This floor is a modeling assumption, not a universal physical boundary; these calculations do not exclude unmeasured censoring mechanisms. Orange bars show required depth and grey bars headroom to the assumed floor.')
    for p in doc.paragraphs:
        if p.text.startswith('Fig. 3.'):
            p.text = p.text.replace('Detectable presence of particular reference families predicts the limit.', 'Detectable reference families are associated with reported minimum growth pH.')
        if p.text.startswith('Fig. 1.'):
            p.text = p.text.replace('computed on four position classes of the same 32 reference proteins', 'computed on three membrane position classes from the same 32 membrane references and a surface class from five separate soluble references')
        if p.text.startswith('Fig. 2.'):
            p.text = p.text.replace('extracytoplasmic composite under four nested models', 'pooled extracytoplasmic net residue share under four adjustment specifications')
        for r in p.runs:
            r.text = r.text.replace('accessibilites', 'accessibilities')
        if p.text.startswith('Fig. '):
            p.paragraph_format.keep_together = True
            p.paragraph_format.line_spacing = 1.15
            p.paragraph_format.space_after = Pt(8)
            for r in p.runs:
                r.font.size = Pt(12)
    for p in doc.paragraphs:
        for r in p.runs:
            r.text = r.text.replace('Fig. 4b', 'Fig. 4B').replace('Fig. 4a', 'Fig. 4A')
    # Neutralise Word's theme title rule and apply the reference's black serif
    # typography to newly inserted headings as well as existing paragraphs.
    for name in ['Title', 'Heading 1', 'Heading 2']:
        style = doc.styles[name]
        style.font.name = 'Times New Roman'
        style.font.color.rgb = RGBColor(0, 0, 0)
        for node in list(style.element.iter()):
            if node.tag == qn('w:pBdr'):
                node.getparent().remove(node)
    for p in doc.paragraphs:
        for node in list(p._p.iter()):
            if node.tag == qn('w:pBdr'):
                node.getparent().remove(node)
        for r in p.runs:
            r.font.name = 'Times New Roman'
            r.font.color.rgb = RGBColor(0, 0, 0)
            fonts = r._element.get_or_add_rPr().rFonts
            for attr in ['asciiTheme', 'hAnsiTheme', 'eastAsiaTheme', 'cstheme']:
                fonts.attrib.pop(qn('w:' + attr), None)
            fonts.set(qn('w:ascii'), 'Times New Roman')
            fonts.set(qn('w:hAnsi'), 'Times New Roman')
    species = r'Escherichia coli|E\. coli|Helicobacter pylori|H\. pylori|Streptococcus salivarius|S\. salivarius|Acidithiobacillus|Picrophilus'
    # Scientific notation uses a complete superscript exponent, including digits.
    superscripts = str.maketrans('0123456789', '⁰¹²³⁴⁵⁶⁷⁸⁹')
    for p in doc.paragraphs:
        for r in list(p.runs):
            r.text = re.sub(r'10⁻(\d+)', lambda m: '10⁻' + m.group(1).translate(superscripts), r.text)
            if not re.search(species, r.text):
                continue
            value = r.text
            r.text = ''
            anchor = r._r
            pos = 0
            for match in re.finditer(species, value):
                for text, italic in [(value[pos:match.start()], False), (match.group(), True)]:
                    if not text:
                        continue
                    new = p.add_run(text)
                    if r._r.rPr is not None:
                        new._r.insert(0, deepcopy(r._r.rPr))
                    if italic:
                        new.italic = True
                    anchor.addnext(new._r); anchor = new._r
                pos = match.end()
            if value[pos:]:
                new = p.add_run(value[pos:])
                if r._r.rPr is not None:
                    new._r.insert(0, deepcopy(r._r.rPr))
                anchor.addnext(new._r)
