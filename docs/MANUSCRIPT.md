# Periplasmic surface charge composition predicts the lower pH limit of bacterial growth, independently of growth optimum

**Status.** Working draft, assembled 2026-11-01. Every **signed coefficient** in this document is checked automatically against the stored result file its own claim cites (`scripts/69`), and all 34 currently trace; a figure here that disagrees with [`NUMBERS.md`](NUMBERS.md) means this document is wrong. **The automatic check does not cover unsigned quantities** — shares, counts, ratios, pH depths — because the sign is part of what gives it power, and §6.1 is written almost entirely in those. Two errors in §6.1 were found by reading it back against the stored file after the automatic run came up clean, which is the honest measure of that gap. Unsigned figures in this document have been checked by hand once, at assembly; they are not protected against going stale. Sections marked **[OPEN]** are not yet written and are listed in §9.

---

## Abstract

How acidic an environment a bacterium can grow in is only partly explained by where it grows best. We asked what genomic features set the lower pH limit (`ph_min`) beyond what the growth optimum predicts, using 2,768 bacterial genomes with a curated minimum and a 32-protein panel of conserved inner-membrane references with per-residue structural confidence filtering.

Three associations emerge. First, the composition of **solvent-exposed periplasmic** positions predicts the limit: a composite of net acidic residue share (D+E−K−R) gives **+0.1848** pH units per standard deviation with genome-wide amino acid composition controlled (p 2.1e−10, `scripts/57`), and the direction is depletion — acid-tolerant bacteria carry *fewer* exposed periplasmic acidic residues. The matched cytoplasmic (**+0.0212**, p 3.3e−01) and transmembrane (**-0.0430**, p 7.3e−02) controls are null (`scripts/57`), and the association is estimated directly as a periplasmic-minus-cytoplasmic difference of **+0.053596** pH per percentage point (`scripts/82 wide`). It keeps **39%** of its size with growth optimum and exposed histidine controlled together.

Second, **histidine** carries a separate and larger association that is not a proxy for the first: the two measures correlate at **r −0.549** and each retains roughly three quarters of its coefficient when the other is held fixed (`scripts/87 wide`). Histidine's effect is partly side-specific (**-0.080003** pH per point periplasmic excess, p 3.1e−03; `scripts/82 wide`) and partly shared with something that eleven candidate explanations fail to identify (`scripts/85 wide`).

Third, the *presence* of particular membrane complexes predicts the limit independently of composition. Detectable breadth across the panel gives **-0.03599** pH per reference (`scripts/91 wide`), and 7 of 22 testable references survive Benjamini–Hochberg correction at q<0.10, strongest cytochrome *bo*(3) ubiquinol oxidase subunit 1 at **-0.42859** pH (`scripts/92 wide`).

The effects are small in variance terms — well under 1% beyond the growth optimum — and the work is correlational. Its principal limitation is the outcome: most reported minima sit on an experimental testing grid rather than a measured biological boundary, and §6 quantifies rather than dismisses what that could do.

---

## 1. Introduction

Bacteria differ enormously in how acidic an environment they will grow in. Some *Acidithiobacillus* and *Picrophilus* isolates grow below pH 1; most described species will not grow below pH 5. What sets that lower bound is not the same question as where an organism grows *best*, and the two are routinely conflated because the optimum is far more often reported. A species whose optimum is pH 7 may have a limit at 4 or at 6, and the difference is the quantity of interest here.

**The mechanisms that are well understood are inducible and physiological.** Enteric bacteria surviving gastric transit deploy amino-acid decarboxylase systems that consume cytoplasmic protons — the glutamate-dependent system through GadBC, the arginine-dependent one through AdiA and AdiC — regulated by sensors that detect acidification directly: the EvgS/EvgA two-component system, and AdiY, whose pH response runs through the protonation of two histidines. *Helicobacter pylori* additionally imports urea through the pH-gated channel UreI to feed intrabacterial urease. These are *responses*: machinery switched on when the environment turns acidic, and they explain survival of acid shock better than they explain the pH at which steady growth becomes impossible.

**A limit might instead be set by what the cell's proteins can tolerate**, and if so the relevant compartment is identifiable in advance. A Gram-negative bacterium holds its cytoplasm near neutrality even in strong acid, so cytoplasmic proteins never experience the external pH. The periplasm does not have that protection: it equilibrates far closer to the medium, and proteins with surfaces facing it work at whatever pH the organism is growing in. Any compositional signature of acid adaptation should therefore appear on outward-facing surfaces and not on cytoplasmic ones — which makes the hypothesis testable by a control that is the same measurement on a different face of the same protein.

**Compositional signatures of other extremes are established; for pH they are not.** Thermophily is visible in the IVYWREL fraction; halophily in a pronounced excess of surface acidic residues. A 2013 review of archaeal extremophile protein adaptation reports those two cases clearly and notes that adaptations toward pH are considerably less discernible. That is a documented absence rather than an assumption, and it is the gap this paper addresses.

**Two recent results make the question tractable.** Cabello-Yeves and colleagues showed in 2019 that predicted isoelectric point shifts with habitat between marine and freshwater prokaryotes, with the difference sharper in secreted than cytoplasmic than membrane proteins and concentrated on protein surfaces — establishing that compartment-resolved composition carries environmental signal at all, and that the approach can be applied to phylogenetically matched pairs. Independently, Barnum and colleagues showed in 2024 that growth pH can be predicted from amino acid composition computed per predicted subcellular localization, reporting the extracellular-minus-intracellular glutamate difference as their strongest residue-level correlate. The two together say that the compartment-resolved signal exists and that pH is among the variables it carries.

**What is asked here.** Using 2,768 bacterial genomes with a curated lower growth limit and a 32-protein panel of conserved inner-membrane references with per-residue structural confidence filtering, we ask which genomic features predict `ph_min` beyond what the growth optimum already explains, and whether any such signal is specific to the compartment the hypothesis names. Three associations survive the controls — exposed periplasmic acidic composition, a separate and larger histidine association, and the presence of particular membrane complexes — and the paper is as much about what bounds them as about what they are. The outcome variable is the weakest part of the work, and §6.1 quantifies rather than excuses it.

---

## 2. Data and methods

### 2.1 Phenotypes

The analysis sample is 2,768 bacterial genomes carrying a curated `ph_min` and a place in the reference phylogeny. Minima come from four source tiers, recorded per row: hand-curated species descriptions read from primary sources, rule-based extraction from species-description abstracts, ThermoBase, and AciDB/BacDive. Source is entered as a set of indicator covariates in every model, so a between-source level difference cannot produce an association.

The extraction tier dominates the sample and its fidelity is measured rather than assumed. Against hand-curated rows it is exact for **98.8%** of minima. That validation was run on an acidophile-heavy subset, so the gap was closed by hand-reading 25 random rows from the band that carries the data, which found 3 of 25 wrong, all by the same mechanism: a sentence reporting only an optimum has that optimum recorded as the growth range. A detector for exactly that failure scores 8 of 8 against fifteen hand-read statements with no false positives, flags **298 of 2,536** parser rows, and dropping precisely those rows moves the headline by **0.13** standard errors (`scripts/60`, `scripts/61`).

### 2.2 The structural panel

The predictor is built on a panel of conserved inner-membrane proteins with experimentally determined or high-confidence predicted structures. For each reference and each genome, residues are assigned to a position class — exposed periplasmic, exposed cytoplasmic, transmembrane, buried — and counted. Positions enter only above a per-residue structural-confidence filter, and a reference contributes to a genome only above a position floor.

The panel used throughout is the 32-reference rule-selected set. A 7-reference hand-picked panel preceded it and is retained for comparison; the two are **not** the same measurement and their per-unit coefficients are never compared. Panel construction is deterministic and verified to rebuild byte-identically (`scripts/79`).

### 2.3 Model

All estimates are phylogenetic generalised least squares with Pagel's λ estimated by maximum likelihood, on a covariance matrix derived from the reference tree. Every model carries an intercept, a genome-wide amino acid composition control, and the source indicators. Coefficients are reported **per standard deviation of the predictor** unless a per-unit figure is explicitly labelled, because composites averaging different numbers of references have different spreads.

Where a model adds the growth optimum, the comparison is made on the single sample for which the optimum is known (1,655 genomes), never against the full-sample figure — mixing a change of model with a change of sample is an error this project has made and corrected.

### 2.4 Verification

Four layers, all of which must pass before a figure is quoted:

- a static audit of nine previously-made error classes in this project's own scripts (`scripts/77`);
- a unit-test suite over the library, including mutation-checked guards on the per-unit/per-SD identity and the optimum-as-range detector;
- a traceability check that extracts every signed coefficient from this and the other documents and looks for it in the result file the same claim cites (`scripts/69`);
- an **independent reimplementation** of the estimator sharing no code with it, described in §3.5.

---

## 3. Result 1 — exposed periplasmic acidic residues

### 3.1 The association and its direction

The composite of net acidic residue share at exposed periplasmic positions gives **+0.1848** pH units per standard deviation (n 2,768, p 2.1e−10, λ 0.573), with genome-wide composition and source controlled (`scripts/57`).

The sign is worth stating explicitly because it is easy to invert in prose. The coefficient is positive, so a *higher* net acidic share goes with a *higher* minimum — that is, with **less** acid tolerance. Acid-tolerant bacteria are **depleted** of acidic residues on their outward-facing periplasmic surfaces, not enriched. At the compartment level the contrast is direct: genomes with `ph_min` ≤ 3.5 average a secreted acidic share of 0.0999 against 0.1246 for those at or above pH 6.

### 3.2 Specificity (Figure 1)

The claim is about a compartment, so the controls are the same measurement computed on other compartments of the same proteins, on the same genomes and in the same model (`scripts/57`):

| Position class | per SD | p |
|---|---|---|
| exposed periplasmic | **+0.1848** | 2.1e-10 |
| exposed cytoplasmic | +0.0212 | 3.3e-01 |
| transmembrane | -0.0430 | 7.3e-02 |
| soluble-surface references | -0.0935 | 1.5e-02 |

The cytoplasmic and transmembrane controls are null. The soluble-surface control is significant and **opposite in sign**, which is not a failure of specificity but does mean the pattern is a contrast between compartments rather than a periplasm-only effect.

### 3.3 One-sidedness, estimated directly

Reporting a periplasmic coefficient beside a null cytoplasmic one does not establish that the two differ, because a difference between two separately-fitted coefficients has no standard error. Entering both compartments in a single model as **sum and difference** makes the periplasmic excess an estimate in its own right: **+0.053596** pH per percentage point (p 3.3e−06), on a reparametrisation that spans the same column space and therefore changes no fitted value (`scripts/82 wide`).

### 3.4 What survives the growth optimum, and histidine (Figure 2)

On the 1,655 genomes whose optimum is known, so that the drop is a change of model alone, the composition-only coefficient of **+0.2123** per SD falls to **+0.0904** with the optimum controlled — **43%** independent of the optimum (`scripts/57`).

The harder test holds the optimum and exposed histidine fixed together, since histidine carries its own association (§4). The coefficient is **+0.035022** per percentage point, **+0.112005** per SD, p 3.894060e−03, retaining **39%** of the composition-only figure on that sample. The two controls turn out to be nearly independent: the product of their separate shares is 0.375045 against a joint share of 0.387212 (`scripts/89 wide`).

### 3.5 Measurement homogeneity, and reproduction

Two threats to the measure itself were tested and neither broke it.

**The reference set is not the same for every genome.** The composite's numerator and denominator are summed over whichever references a genome has a detectable homologue of; detection runs from 54.6% to 100% and breadth averages 28.8 of 32, so two genomes with identical surface chemistry can receive different shares. Holding breadth fixed in a regression does not remove this, because it is heterogeneity in what the measure *measures*. Two routes that fail in different places were run: a **6-reference common core** identical across genomes gives net **+0.03620** per point (p 0.012), and the **239 genomes carrying every reference** give **+0.08880** (p 0.036). Neither failed (`scripts/93 wide`).

**The estimator had never been checked against anything but itself.** An independent implementation — its own recursive-descent Newick parser, covariance from shared root-to-tip path length by prefix comparison, symmetric eigendecomposition with explicit normal equations and a 4,001-point λ grid, no Cholesky and no Brent search — reproduces the headline at **+0.184838** per SD against the published **+0.184837**, both matched controls, and both sided results to within 2e−7 per point, across 30 declared comparisons (`scripts/95`).

### 3.6 Did it arise once or repeatedly?

A phylogenetic regression prevents shared ancestry being counted as independent evidence, but it does not count how many times the thing happened. The coefficient is equally consistent with a shift that accompanies acid tolerance wherever it arises and with one ancient shift in a clade that happens to contain most acidophiles — both real associations, meaning different things.

Treating acid tolerance as a tip state at `ph_min` ≤ 4.0, the 230 tolerant genomes fall into **151** maximal wholly-tolerant clades, of which **33** have a sister clade carrying at least three non-tolerant scorable tips. Comparing each tolerant clade against its own sister — composites residualised on genome-wide composition first, as the sister-pair work showed is necessary — the acidic measure points the predicted way in **19 of 33** contrasts with a median difference of **-0.0379**, which is not significant (signed-rank p 0.126). The informative figure is the power: detecting the effect at this sample size would need the predicted direction in **75%** of contrasts, and 58% was observed. **This is underpowered, not null** (`scripts/97`).

Across five thresholds the direction is stable — predicted-way at 5 of 5 — and significant at 4 of them, without the primary threshold being one of those four. The honest reading is that direction is consistent while significance is not, which is what an underpowered real effect looks like and also what a weak artefact looks like; adopting whichever threshold cleared would be selection, so the primary stands.

**[OPEN]** This is the weakest of the three findings' supporting tests and the one most worth strengthening. A proper ancestral-state reconstruction would separate gains from losses, which tip states cannot.

---

### 3.7 It does not predict in clades it has never seen

Leave-one-phylum-out refitting (§3.2's lineage work) asks whether any single lineage carries the coefficient; it does not ask whether the model transfers, because every such fit still had the held-out phylum contributing to panel selection, composite standardisation and covariate structure.

Fitting on every other phylum, with the composite re-standardised on training genomes only and λ estimated there, pooled out-of-sample RMSE is **0.9253** pH with the composite against **0.9248** without it — fractionally *worse* — and the cytoplasmic control, at **0.9238**, does no worse than either (`scripts/102`).

**Two readings must both be carried.** The composite lowers error in 3 of 4 folds, and the pooled figure is decided by the one it fails: Pseudomonadota, 49% of scored genomes and the phylum holding most of the acidophiles, whose removal leaves training with little of the acid-tolerant tail the effect is confined to. Fold-level gain tracks the tail retained in training at **r +0.835** (`scripts/102`) — a direction predicted in advance, on four points, and offered as nothing more. So the association survives every within-sample test in this paper and does not buy out-of-sample accuracy across clades. For a claim about mechanism that is survivable; for a claim about prediction it would be fatal, and this paper makes none.

Reference selection remains fitted on everything — re-running the panel choice per fold needs the structural pipeline, not the stored table — so even these figures are optimistic by an unmeasured amount.

---

---

## 4. Result 2 — histidine is a separate association, not a proxy

### 4.1 The association, and its sidedness

Exposed periplasmic histidine share carries a larger coefficient than the acidic composite and in the opposite direction. Entered as sum and difference across the two compartments, with periplasmic D, E and K held fixed, the shared component is **-0.131680** pH per percentage point and the periplasmic excess **-0.080003** (p 3.1e−03) — so the association has both a side-specific and a shared part. That contrast is calibrated on the acidic result first, which fires at **+0.053596**, because a test is only a test once it has been shown it can fire (`scripts/82 wide`).

### 4.2 Not a proxy for the acidic measure

The two measures correlate at **r -0.549173**, which is large enough that one standing in for the other is the obvious worry. Mutual adjustment rules it out: each retains roughly three quarters of its coefficient with the other held fixed, and the acidic measure's one-sidedness survives at **77%** of its unadjusted value (`scripts/87 wide`).

### 4.3 Not pKa, but charge class and solvent access

An ordering by side-chain pKa was pre-registered and does not hold monotonically. What does organise the coefficients is charge class, which remains real with histidine held fixed (**+0.08681** per point, p 1.6e−04), and pKa orders the acidic pair (**+0.03958** per point, p 3.9e−03) while failing to order the basic ones. Histidine itself runs **8.6×** its charge-class mates (`scripts/90 wide`).

Solvent access is doing work: histidine at **non-exposed** positions in the same 32 families predicts nothing, alone (p 0.37) or jointly (p 0.47). A histidine no proton can reach carries none of this, which is the falsifier the mechanism needed and clears (`scripts/83 wide`).

### 4.4 It also travels with the trait across independent clades

The transition test of §3.6 is the one place histidine outperforms the acidic measure. Across the same 33 independent tolerant-versus-sister contrasts it shifts in the predicted direction in **22 of 33**, median **+0.1772**, signed-rank p 0.046 (`scripts/97`) — the only measure in that test to reach significance at the primary threshold, and in the direction fixed in advance from its published sign.

Specificity is **not** established there, and the reason is worth stating plainly rather than rounding away: histidine's own matched cytoplasmic control comes in at p 0.0506, which is over the conventional line by an amount no specificity claim should rest on. The compartment claim for histidine rests on the directly-estimated sided contrast of §4.1, not on this test (`scripts/97`).

### 4.5 The shared component is unidentified

Correction 19 split the histidine association into a periplasmic excess and a larger shared component, and what that component is shared with is not known. Twelve dataset properties were listed, the eleven this dataset supports were tested against a criterion fixed before the run, and four predict `ph_min` in their own right — GC content, genome size, protein count, mean membrane-protein pI — while none accounts for the component: the largest reduction any live candidate manages leaves **85%** of it. Thermophily, the most plausible candidate since histidine is a thermal-adaptation residue, leaves it at **100%** by both the growth optimum and the IVYWREL proxy (`scripts/85 wide`).

### 4.6 Charge burden against buffering, measured

Correction 21 left the mechanism in question: UreI needs its periplasmic carboxylates, while acid-tolerant genomes shed theirs. The two accounts are separable from composition alone, because they want opposite signs of the same quantity. Writing Q(pH) for net surface charge and β(pH) = −dQ/dpH for buffer capacity, **burden** says acid tolerance should come with a smaller charge change between neutrality and the limit; **buffering** says it should come with greater capacity at the operating pH.

At pH 3 and 4 — inside the acidophile range and near the carboxylate pKa values — the charge change carries a **positive** coefficient (**+0.13994**, p 3.8e−08; **+0.07764**, p 1.1e−03), so acid-tolerant genomes have less outward-facing charge to change, as burden predicts. Buffer capacity carries a positive coefficient too, which is the opposite of what buffering predicts (`scripts/104`).

**Both quantities reverse sign between pH 4 and pH 6** (**-0.20915** at pH 6; `scripts/104`). That is where histidine titrates, and histidine's association runs opposite to the acidic one throughout this paper — so the reversal is the second mechanism appearing inside the same quantity rather than instability in the first.

**Two caveats, one of them a discarded analysis.** A first version evaluated both quantities at each genome's *own* `ph_min`, which puts the outcome inside the predictor; the cytoplasmic control — about which neither account predicts anything — fired at p 0, and that part is discarded rather than reported. And this is not an electrostatics calculation: no structure, no solvent, no interaction between groups, and model pKa values throughout. It is the first-order question a solver would refine, worth asking first because the two accounts differ in *sign*, which a first-order treatment can see.

So the reconciliation of §4.7 now has a measurement behind it and not only an argument: a gate that must respond to acid needs titratable groups, and a bulk housekeeping surface appears to do better with fewer.

### 4.7 Independent experimental support, and where it is absent

The mechanistic literature was searched for direct manipulations of the two residue classes, and the result is asymmetric in a way worth stating before any reader finds it themselves ([`EXPERIMENTAL_EVIDENCE.md`](EXPERIMENTAL_EVIDENCE.md)).

**Histidine is supported from three independent directions.** In *E. coli* EvgS, an acid-sensing histidine kinase, a single histidine in the periplasmic domain is required for acid induction. In *H. pylori* UreI, pH-dependent gating has been attributed to periplasmic histidines. In *E. coli* AdiY, replacing two histidines abolishes pH-dependent activation in vivo and in vitro. Two of the three are periplasmic, and all three establish histidine protonation as a working pH-sensing mechanism in bacteria.

This supports the residue, not the comparative claim. These are single residues in dedicated sensor proteins; this project measures a compositional share averaged over 32 housekeeping membrane proteins, and nothing in that literature says such a share should predict a growth limit. A dedicated pH sensor and a bulk surface property are different hypotheses that happen to involve the same amino acid.

**The acidic prediction has one direct experiment, and it confirms the compartment while contradicting the direction.** Weeks & Sachs mutagenised *every* histidine, aspartate and glutamate in the periplasmic domain of the *H. pylori* urea channel UreI — the exhaustive periplasmic titratable-residue scan this prediction concerns — and found His-123, His-131, His-193 **and Asp-129, Asp-140, Glu-138** all needed to activate transport at acidic pH, with those residues largely absent from the pH-independent *S. salivarius* homologue. There, periplasmic carboxylates are *required* for low-pH function; here they are *depleted* in acid-tolerant genomes.

**The reconciliation is specific and should not be glossed.** A protein whose job is to respond to acid needs titratable residues on the face that senses it; a proteome that must merely function at low pH may do better with fewer. UreI's residues are a dedicated gate in one specialised channel, while this project measures bulk composition over 32 housekeeping membrane proteins and relates it to a growth limit. Those are different quantities, and a trend in one does not contradict a mechanism in the other.

**What it does change is the mechanism, not the association.** The charge-burden story offered here — carboxylates protonate, outward-facing surface charge shifts — is not the account the direct evidence favours. UreI needs its carboxylates, and the one engineering study that targets gastric pH *adds* acidic residues on an explicit buffering rationale. Two further studies on soluble enzymes measure in vitro stability rather than a growth limit and are directionally mixed. So the association stands on its own evidence while its stated explanation is weaker than the data behind it, and that asymmetry belongs in the discussion rather than in a footnote ([`EXPERIMENTAL_EVIDENCE.md`](EXPERIMENTAL_EVIDENCE.md)).

**This yields a sharp follow-up the project cannot currently run.** If the sensors-versus-housekeeping reconciliation is right, a panel built from *dedicated pH sensors* — EvgS, CadC, UreI, none of which is in the 32 — should show the opposite sign to the housekeeping panel. That is a falsifiable prediction of the reconciliation, and it needs structures and detection for a sensor panel this project does not have.

---

## 5. Result 3 — particular membrane complexes

### 5.1 Gene content, separated from measurement artefact

How many of the panel's reference families a genome has a detectable homologue of was first held fixed as a nuisance control. Decomposed, it is two different things. **Breadth** survives at **-0.03599** pH per reference (p 6.7e−04) with assembly level, similarity to the reference proteins and genome size all held fixed; **depth** — positions per homologue — carries nothing once similarity is controlled (p 0.93). So the measurement half is correctly a control and the gene-content half is a secondary association (`scripts/91 wide`).

### 5.2 Which complexes (Figure 3)

Of the 37 accepted references, **22** have enough detection variance to regress on, and **7** survive Benjamini–Hochberg at q<0.10 with composition, source, mean alignment score, log genome size and assembly level all controlled. The strongest is cytochrome *bo*(3) ubiquinol oxidase subunit 1 at **-0.42859** pH, with four transporters the same way and two the other (`scripts/92 wide`).

### 5.3 Three limits, stated

**The intended negative control cannot run.** The five soluble references were to be the null — essential, universally present, predicting nothing — and they *are* universally present, at 99.1% to 100.0% detection, which is why every one falls below the variance floor and none was tested. A control that cannot run is a limitation, not a pass. What replaces it is weaker and labelled so: the survivors' signs are **mixed**, five down and two up, which rules out a shared sequencing-or-annotation artefact (that account predicts one direction for every reference) and establishes nothing about function.

**Detection is a proxy for presence.** A divergent homologue below the gathering threshold reads as absent, so these are claims about detectable homologues of *E. coli* proteins, which the alignment-score control bounds but does not remove.

**The concentration result sits at the edge of its own rule.** Breadth keeps **47%** of its coefficient with the survivors held fixed, against a criterion of half, so the named references carry about half of the breadth effect rather than being it (`scripts/92 wide`).

---

## 6. Limitations

### 6.1 The outcome is the weakest part, and it is quantified (Figures 4a, 4b)

`ph_min` is meant to be the lowest pH at which an organism grows. What a species description usually records is the lowest pH at which growth *was seen*, bounded below by whatever the authors chose to test. If that choice tracks what they expected, the error is not noise: it is related to the organism's apparent acid tolerance and therefore to anything predicting it.

**The reported minima do carry a testing-grid signature.** **74%** sit on a whole pH unit and **93%** on a whole or half unit, and the modal value — pH 6.0, **34%** of the sample — is a standard medium's lower end rather than a biological boundary. A measured limit has no reason to prefer a round number; a testing grid has every reason to (`scripts/96`).

**The field that appears to record bound type does not.** `ph_min_bound` reads `measured` for 2,732 rows, but the table-building script assigns that value unconditionally to every parser-extracted row; **97%** of the rows carrying it arrive that way, never having been read. Only the hand-curated tier had bound type determined from text, and 4 of the 90 such rows surviving into the analysis sample came back censored. Three strata reported as separate robustness checks elsewhere in this project — bound type, confidence, and source — agree on over 95% of rows and are therefore one partition reported three times (`scripts/96`).

**The existing reassurance addresses a different error.** Measured outcome reliability (ICC 0.965, error SD 0.343 pH) supports the standard result that error in the dependent variable inflates standard errors without biasing slopes — which is true of *classical* error, independent of the predictor by construction. Censoring is not classical: it has a sign and a mechanism. Simulated at equal magnitude, classical noise moves the coefficient by **0.006** standard errors (mean over 200 draws) while the same magnitude aligned with the predictor moves it by **11.8** (`scripts/96`).

**So the question is bounded rather than dismissed.** To drive the coefficient to zero, the censoring displacement would have to deepen by **0.185** pH for every standard deviation of surface composition. Routed through observable proxies for being censored, that requirement is not merely implausible but **physically impossible**: the two proxies pointing in the threatening direction require mean displacements of **8.36** and **16.21** pH against **6.04** and **6.00** pH of available headroom — the flagged organisms' true minima would have to lie below pH 0 (`scripts/96`).

**What that does and does not settle.** It closes the objection in the form a proxy can express it. It does not close the objection: these proxies are weak and each proxies other things too, and a selection process acting on *which organisms get described at all* would appear in none of them. Settling it needs bound type read from full texts at scale — which cannot be mined from a literature index, because negations are indexed as stopwords (`scripts/65`) — or growth limits measured on a common grid, which is bench work.

### 6.2 Effect size

The associations are small in variance terms: well under 1% of `ph_min` variance beyond what the growth optimum explains. The relationship is also **one-sided** in a second sense — concentrated in the acid-tolerant tail rather than running across the whole range — so the coefficient should not be read as a uniform slope.

### 6.3 A second environmental variable with the same compartment signature — tested

The compartment-resolved approach has a published precedent for a *different* variable, and it was the sharpest specificity threat found. Marine-to-freshwater prokaryotic transitions shift predicted isoelectric point, with the difference sharper in secreted than cytoplasmic than membrane proteins and concentrated on protein surfaces — the same ordering reported here for pH, by a different group, six years earlier. The matched controls used throughout this paper rule out *position-class* artefacts; by construction they cannot rule out a second environmental variable acting through the periplasmic face, because such a variable would move the periplasmic class and leave the controls alone, which is exactly the pattern reported here as evidence.

**It is now tested.** A measured-only trait export supplies `salinity maximum` in % NaCl for 5,199 species, 669 of which join the analysis sample. Salinity **does not predict `ph_min`** (p 3.6e−01), and holding it fixed removes **1%** of the periplasmic coefficient — measured against a baseline refitted on those same 669 genomes, never against the full-sample figure (`scripts/98`).

**The compartment profiles are separable, which matters more than the covariate adjustment.** Entered side by side as predictors of proteome compartment composition, salinity gives secreted > cytoplasmic > membrane — *the 2019 ordering, reproduced here independently* — while `ph_min` gives secreted > membrane > cytoplasmic. Salinity moves the secreted and cytoplasmic shares together, as a whole-cell osmotic account predicts; pH moves the secreted share and leaves the cytoplasmic one untouched, as a compartment-specific account predicts. The method does not impose one ordering on every variable it is given. That is a stronger defence of the controls than salinity being inert would have been — and salinity is not inert: it predicts secreted composition *more* strongly than pH does.

**Two limits, and one of them is a new open question.** The salinity values are species-level attached to strain genomes, the subsample is 24% of the sample and selected on being measured, and habitat is more than % NaCl — so a null here means "salinity as measured this crudely does not account for it". Separately, on these 669 genomes the transmembrane control *fires*, at **-0.14768** (p 3.0e−03; `scripts/98`), **1.06×** the periplasmic coefficient against 0.23× on the full sample. **That is salt adaptation, and it is now identified** (`scripts/103`): the signal is confined to the salt-tolerant half — **-0.21178** (p 1.4e−03) above 5% NaCl against **-0.06218** (p 3.4e−01) below it — while an independently-defined well-characterised subset does not reproduce it (optimum-known, n 1,655, p 0.44) and no single reference moves it more than 0.76 SE. Halotolerant membranes are compositionally distinct, which is what the archaeal halophile result already showed. A second environmental variable perturbs a *different* compartment from the one claimed here, which is a result in favour of the controls rather than against them. Across 300 random subsamples of the same size only **2.0%** reach a transmembrane coefficient this extreme, so it was never smaller-n noise.

### 6.4 Not causal

Nothing here is interventional. A pre-registered charge-reversal experiment is designed but not run, and the required replicate count rose from 64 to **199** once two separately-recorded rescalings of the predicted effect were multiplied, which had never been done ([`CHARGE_REVERSAL_DESIGN.md`](CHARGE_REVERSAL_DESIGN.md), `scripts/88`).

### 6.5 Archaea are consistent but marginal

On the wide panel with halophiles excluded (n 177) the exposed extracellular composite gives **+0.153** per SD (p 0.065) — not significant, but 0.4 archaeal standard errors from the bacterial estimate and unable to exclude zero either. Over all 318 with a halophile indicator it is **+0.237** (p 0.001) with every matched control null. Within the 141 halophiles specificity fails and remains confounded with alkaliphily (`scripts/43`).

### 6.6 Twenty-one corrections

Twenty errors in this project's own analyses have been found, recorded and propagated, and the great majority belong to two families: a residue share entered without the correlate it shares variance with, and a figure compared against one computed on a different sample. Both families are now checked mechanically.

A twenty-first is of a different kind and is the only one that bears on this document's literature claims. §4.6's first version asserted that no experiment had manipulated periplasmic acidic residues in a bacterial membrane protein and measured acid resistance. That rested on a semantic search returning three results per query; an exhaustive PubMed query on the same question returns 59 records, and one of them is the most directly relevant experiment in the project. An absence claim drawn from a three-result search is not an absence claim, and nothing in this project's checking apparatus covers literature the way it covers figures.

All twenty-one are listed in [`SYNTHESIS.md`](SYNTHESIS.md) §5 rather than summarised here, because a summary of one's own errors tends to be kind to them.

---

## 7. What would change the conclusions

Stated as falsifiers rather than future work:

- the periplasmic coefficient failing to survive a reference set that is identical across genomes — tested, and it survived by two routes (§3.5);
- the cytoplasmic or transmembrane control coming back comparable in size — tested, both null (§3.2);
- non-exposed histidine carrying the histidine association — tested, it carries none (§4.3);
- the breadth effect concentrating entirely in a few references, making it a gene-presence story rather than a composition one — partially the case, at about half (§5.3);
- a demonstrated-limit subset reversing the sign. **Not tested, and not testable from these sources**; this is the one falsifier with no result attached.

---

## 8. Reproducibility

The package is [`RELEASE.md`](RELEASE.md), generated by `scripts/101_release_manifest.py`, and it opens by stating what is *not* reproducible rather than burying it: the full pipeline does not re-run from the repository, because its inputs include 12 GB of genome assemblies and 4.3 GB of proteomes that are gitignored, and the phenotype curation involved hand-reading primary literature that no script replays.

What is reproducible is tiered. **Tier 1 needs nothing but the repository** and takes about fifteen seconds: it verifies that every signed coefficient in every document matches the stored result its own claim cites, that none of nine previously-made error classes is mechanically present, and that the library passes 26 tests including mutation-checked guards on the per-SD/per-unit identity. **Tier 2** re-runs every analysis script against the committed results. **Tier 3** — panel construction, HMM detection, phylogeny — needs re-fetched inputs, and each script that needs one carries its own fetch command in a header comment.

The manifest pins the environment, and records a SHA-256 for each result file the documents cite, so a reader checking a quoted coefficient can confirm the file they have is the file it was checked against. Of 108 declared result files, 89 are present and checksummed; the remaining 19 are declared by scripts whose outputs are not in this dataset and are listed as absent rather than omitted.

All five figures are redrawn from stored results by `scripts/100_figures.py`, and each panel asserts its headline values against this document's text before drawing, so a figure cannot drift from the prose it illustrates.

---

## 9. What is not yet written

| Section | What it needs |
|---|---|
| The sensor-panel prediction (§4.6) | **Attempted and unrunnable** (`scripts/99`). Dedicated periplasmic sensor domains are detected in 1.9% of the sample and 100% of those detections sit in one phylum, so a sensor-versus-housekeeping contrast would be a between-lineage contrast in disguise. The prediction is untested, not refuted. Needs a sensor set chosen for breadth rather than for being the textbook examples. |
