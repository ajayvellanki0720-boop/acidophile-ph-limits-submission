# Analysis plan: genomic determinants of the lower pH limit in acidophiles

Adapted from the halophile (S_max) plan. Everything pH-specific is called out in **Differences from the halophile plan**.

## Outcomes
- **Primary:** `ph_min`, lowest pH permitting growth (lower = more acid-tolerant).
- **Comparison covariate:** `ph_opt`. Every genomic model must beat the `ph_opt`-only model on held-out lineages.
- **Secondary:** `ph_max`, used with `ph_min` to separate specialists (narrow range) from acid-tolerant generalists; `acid_reach = ph_opt - ph_min`.
- **Mechanistic variable:** which acid-resistance strategy is in use, inferred from system content rather than assumed from taxonomy (strategies are not mutually exclusive):
  1. *Donnan / inside-positive potential:* K+ uptake (Kdp, Trk/Ktr) builds a reversed membrane potential that repels protons.
  2. *Low-permeability membrane:* tetraether (GDGT) lipids, hopanoids, cyclopropane fatty acids.
  3. *Proton consumption and extrusion:* Gad/Adi/Cad decarboxylases, urease, arginine deiminase, proton-pumping respiration, F/V/A-type ATPases.
  Acid-tolerant neutralophiles lean on (3); obligate acidophiles are thought to lean on (1) and (2). That is a hypothesis to test, not an assumption.

## Stage 1. Phenotype-genome dataset (`src/acidophile/phenotype.py`)
Per strain: ph_min/opt/max, each with a **bound type** (`measured`, `lowest_tested`, `highest_tested`, `extrapolated`, `unknown`); acid identity (H2SO4, HCl, organic, mixed); organic-acid identity and concentration; temperature **at which pH was measured**; growth temperature, oxygen, medium, metabolism, metal load; source and confidence; exact matched NCBI assembly. "Grew at pH 1.0, lower not tested" is a censored bound and is excluded from the primary set; a sensitivity check asks whether predictions contradict it (`censored_violation_rate`).
Acid chemistry replaces water activity: equal pH is not equal stress. Weak organic acids cross the membrane in the HA form and uncouple the proton gradient (`undissociated_fraction`), and low pH mobilises metals, so metal load is a covariate.
Source: BacDive (ranges and optima), with literature curation for bound type and acid identity. Acquisition is manual/API with credentials and is not scripted here.

## Stage 2. Broad genomic signatures (`proteome.py`)
Acidic excess (D+E)-(K+R) with **sign unassumed**, predicted pI (mean, median, fraction >7 and <5, plus membrane vs soluble subproteomes), GRAVY, Cys and His, GC, dinucleotide odds ratios, genome size, membrane-protein fraction (Kyte-Doolittle proxy, to be replaced by DeepTMHMM). IVYWREL enters as a **thermal control** because many acidophiles are thermoacidophiles. Accessory-genome size and codon usage need pangenome and CDS inputs and are deferred.
Question asked: do these form a useful acid-proteome index, without assuming one score explains the exact limit.

## Stage 3. Systems (`systems.py`)
Nine categories, 28 systems: Donnan K+ (Kdp, Trk/Ktr, Kef); proton handling (F/V/A ATPases, Complex I, oxidases, P-type ATPases); acid consumption (Gad, Adi, Cad, urease, ADI, carbonic anhydrase); membrane barrier (GDGT ring synthase, hopanoid SHC, CFA); envelope (S-layer, EPS); proteostasis and repair; oxidative and metal resistance (copper, arsenic, CDF); acid-only energetics (rusticyanin/Cyc2 iron oxidation, sulfur oxidation); weak-acid catabolism.
Detection is annotation-keyword based and **provisional**. Before any system is used in a claim: HMM or reference-sequence validation with a coverage check (the thermophile audit showed keyword and single-hit calls produce false positives).

## Stage 4. Models (`models.py`)
`ph_min ~ ph_opt + proteome composition + pathway/gene features + acid chemistry + temperature + domain`.
M1 ph_opt only; M2 acidic proteome; M3 composition; M4 gene content; M5 pathway; M6 combined. M2 to M6 are also run without ph_opt, asking whether genome can substitute for it. Ridge with in-fold imputation and scaling. Validation is **family-blocked** and **leave-one-lineage-out**; the code asserts no group appears on both sides of any split. ΔRMSE vs M1 gets a bootstrap CI over whole groups. Bacteria/Archaea and acid class are analysed separately; strata below the minimum size are reported as skipped, not pooled.

## Stage 5. Lineage-specific limits (`lineage.py`)
Implemented now: matched close-relative contrasts (same genus or family, |Δph_opt| ≤ 0.5, |Δph_min| ≥ 1), averaged per taxon so no genus dominates, sign test and Wilcoxon with Benjamini-Hochberg across features. Deferred until a marker tree exists: reconstruction of independent increases in acid tolerance and testing whether the same changes recur.
Positive-control lineages to include: *Acidithiobacillus* and *Leptospirillum* (bacterial iron/sulfur oxidisers), *Sulfolobales* and *Thermoplasmatales* (thermoacidophilic archaea; *Picrophilus* sits at the extreme), *Acidiphilium*, and acid-tolerant *Escherichia/Lactobacillus* for the proton-consumption route.

## Criteria for a convincing result
A system contributes to pH_min only if it (1) improves held-out pH_min over ph_opt, (2) survives phylogenetic, temperature, acid-chemistry and metal controls, (3) recurs across independent transitions or is strong within a lineage, (4) predicts the boundary rather than just identifying acidophiles, and (5) has a plausible failure mechanism at low pH.

## Differences from the halophile plan
| Halophile plan | Here | Why |
|---|---|---|
| S_max, upper limit | pH_min, lower limit | acid tolerance is the low-pH boundary; pH is already logarithmic so no transform |
| Salt identity and water activity | Acid identity, weak-acid HA fraction, metals | the stress is proton gradient plus uncoupling plus metal mobilisation |
| Salt-in vs salt-out | Donnan vs membrane barrier vs proton consumption | strategies are overlapping, so inferred from systems, not labelled |
| Acidic proteome expected | Sign of acidic excess unassumed | direction in acidophiles is not established |
| none | IVYWREL thermal control, temperature of pH measurement | thermoacidophiles confound acid and heat |
