**The analysis for this experiment is pre-registered in [`docs/PREREGISTRATION.md`](PREREGISTRATION.md), written before any strain was built. Where this document and that one disagree about an analysis choice, the pre-registration governs.**

# Charge-reversal scan: experimental design

Test of the one robust finding in this project: across 2,768 bacteria, acidic residues on the **exposed periplasmic** surfaces of conserved membrane proteins track a **higher** lower growth-pH limit (less acid tolerance), while buried, cytoplasmic-facing, transmembrane and soluble-surface positions do not. The association survives family-mean-centring, the hand-curated evidence tier, and controls for genome-wide composition, GC, genome size and growth optimum. It is still an association. This document specifies the experiment that would make it causal, and states in advance what would refute it.

Supporting analyses: `scripts/36a`–`36d` (panel, composition, distribution, robustness), `scripts/37a`–`37b` (per-position scan, target selection). Target list: `results/broad2/mutagenesis_targets.csv`.

## 1. Direction of the prediction

The slope is **positive**: more acidic excess on exposed periplasmic surfaces → higher pH_min. Therefore:

- **Removing** negative charge there (D/E → K) should **lower** pH_min — *more* acid tolerant.
- **Adding** negative charge (K/R → D/E) should **raise** pH_min — *less* acid tolerant.

Both arms are testable in the same organism, which matters: a gain-of-function and a loss-of-function result together are far stronger than either alone.

## 2. Expected effect size, and why the naive number is wrong

A single D/E→K reversal changes that protein's exposed-periplasmic acidic excess by `2/n_positions`. Applying the per-protein comparative slope directly gives a much larger prediction than the composite basis does — reversing all 5 selected NuoN positions would shift pH_min by about 0.2 pH units. The comparative slope is estimated across species that change **many** proteins in a correlated way, so it absorbs the effect of the whole coordinated shift; one protein out of thousands should deliver a small fraction.

**Corrected 2026-10-31 (synthesis correction 15).** This paragraph used to say 1.3 pH units and call the naive transfer "absurd", and the table's naive row said −1.9 for the primary arm. Those figures multiplied a share change in percentage points by a coefficient that is per *standard deviation* of that share, so they were inflated by roughly that SD — about sixfold (`scripts/74`). On the correct scale the naive transfer gives about **−0.23** for the 14 reversals against the composite basis's −0.05, so the gap that justifies working through a composite is about **five- to sixfold, not thirtyfold**. The argument is unchanged and so is the conclusion; the rhetoric was carrying a units error. **The committed predicted effects are not affected**: the chain below converts points to a z-shift by dividing by the share's own SD before applying the per-SD composite slope, so it never multiplies by a per-point coefficient, and both committed figures and the *n* ≥ 64 derived from them (`scripts/70`) stand as they are.

The defensible estimate uses the composite predictor from `scripts/36d`, which averages z-scored acidic excess over ~6 references. Changing one protein moves the composite by about one-sixth of that protein's own z-shift.

**The predicted per-reversal effect is not identified by the comparative data.** This is the honest statement, and it replaces the single figure this section used to give. The problem: a composite that averages *k* references moves by only 1/*k* of one protein's own z-shift when that protein is mutated, while the slope per SD of that composite *rises* with *k* as averaging removes measurement error. The product does not converge, and *k* is an arbitrary modelling choice.

| Basis | Slope per SD | NuoN, per reversal | NuoL + NuoN, 14 reversals |
|---|---|---|---|
| Per-reference slope applied directly, no composite | — (+0.0336 per SD of the share, i.e. +0.0053 per point) | −0.024 | **−0.23** (was −1.9; see the correction above) |
| 7-reference composite | +0.101 | −0.012 | −0.16 |
| 32-reference composite | +0.185 | −0.004 | −0.05 |

**These two tables' per-reversal and combined-variant figures are derived in this document**, by multiplying a stored slope by position counts, so no result file contains them and `scripts/69` reports them as unsourced — correctly. The inputs are cited: the slopes in the row below, the composite in the paragraph above, the position counts in Section 5, and the first row's +0.0336 per SD and +0.0053 per point in `scripts/74`. Citing the table ROWS would have been worse than leaving them bare — a row carries both its stored slope and its derived prediction, so the citation verifies the slope and reports the prediction as a miss. A citation here would point at a file that does not hold these numbers, which is worse than none.

The first row is the naive transfer rejected below; the last two bracket the defensible range. **Design for −0.05, treat −0.16 as the optimistic end, and require the assay to resolve 0.02.** A measured effect far above this range would point to folding or expression damage rather than to the hypothesis.

**Rescaling pending, 2026-11-01 (synthesis correction 20, `scripts/88 wide`).** These figures are computed from a coefficient that does not hold exposed-periplasmic histidine fixed, because the panel carried no histidine column until 31 October. The experiment should inherit the controlled coefficient, and the reason is what the edit does: reversing an aspartate or glutamate to lysine changes that class's acidic and basic shares and leaves histidine exactly as it was, so the causal quantity is the acidic effect with histidine held fixed. On the composite this table scales through, holding it fixed takes `net D+E−K−R` from **+0.07623 to +0.05453 per point**, a factor of **0.7153** — a **28.5%** reduction, with both fits highly significant, so this is a rescaling and not a threat to the premise. Applied as a ratio rather than a substitution (the committed figures descend from a different model, and swapping absolute coefficients would be correction 12), the rows above become: per reversal **−0.003** and **−0.009**, primary strain **−0.036**, optimistic end **−0.114**. The primary strain stays above the 0.02 resolution target but at **1.8×** rather than 2.5×, and because replicates scale as the inverse square of the effect, `scripts/70`'s **n ≥ 64** becomes **n ≥ 126**. **Not adopted in this table**: amending committed predictions is the deviations procedure's decision, and `scripts/88 wide` exists to give it the number. Note the direction — Section 11 declines a 1.55× *increase* partly because raising a predicted effect is the direction that flatters the experiment; this lowers it, which removes the usual reason for suspicion without earning any discount on scrutiny.

**On the slope itself: +0.101 to +0.185 pH per SD, depending on the panel** (`scripts/39a`, `scripts/57`, corrected 2026-10-26; the endpoints were previously rounded to +0.10 and +0.19, which left them matching no stored value and so unverifiable — `scripts/69`'s own advice is that quoting more decimals is what gives it power). The 32-protein panel gives +0.185 and the 7-protein panel +0.101; the difference is attenuation, not disagreement. The apparent conflict between a "within-family +0.312" and a "between-species +0.101" was an error of mine, not a real disagreement: the first was a slope per *unit* of the composite and the second per *standard deviation*, and the composite's SD is 0.53, not 1. On one scale every estimator agrees — within-family OLS +0.112, PGLS +0.107, PGLS with family fixed effects +0.088 per SD, all with genome-wide composition controlled. There is no upper-bound variant to design against, so the single slope below is used throughout. Reversals are also capped by how many acidic positions each protein actually has, which the first version of this table ignored:

Per-protein, with every Δ column **derived here** from the slope above and the position counts — like Section 2's predictions, no result file holds them, so `scripts/69` reports them unsourced rather than verified, and a citation on these rows would point at files that do not contain them:

| Protein | Exposed peri. positions | Δ per reversal | Δ composite per reversal | ΔpH_min per reversal | D/E available | Targets passing §5 | **ΔpH_min at max** |
|---|---|---|---|---|---|---|---|
| NuoN | 43 | 4.65 pts | 0.120 SD | −0.012 | 6 | 5 | **−0.06** |
| NuoL | 85 | 2.35 pts | 0.105 SD | −0.011 | 11 | 9 | **−0.10** |
| YidC | 122 | 1.64 pts | 0.033 SD | −0.003 | 15 | 11 | **−0.04** |

Three consequences, each of which changes the design:

1. **Single reversals are unmeasurable** (−0.003 to −0.012 pH units). They can only ever serve as the bottom of a dose series, never as a comparison in their own right.
2. **No single protein is enough.** Even saturating every available position gives only −0.04 to −0.10, at or below the resolution floor.
3. **Composite shifts are additive across references, so multi-protein variants carry the test.** Combining proteins is the only way to reach a comfortably detectable effect:

| Combined variant | Reversals | ΔpH_min, 7-ref basis | ΔpH_min, 32-ref basis |
|---|---|---|---|
| NuoL + NuoN | 14 | −0.16 | **−0.05** |
| NuoL + NuoN + YidC | 25 | −0.19 | −0.06 |

The two bases differ threefold and nothing in the comparative data chooses between them, so the smaller column sets the sensitivity requirement.

**Design target: reliably resolve 0.02 pH units.** The earlier target of 0.05 was set against effect sizes that have since been revised down twice, and on the 32-reference basis the whole predicted effect of the primary strain is about 0.05, so 0.05 resolution would leave it undetectable. This rules out growth curves as the primary readout and makes the pooled competition assay of Section 7 mandatory rather than preferable. **The margin is thinner than when that was written** (correction 20, `scripts/88 wide`): on the rescaled 32-reference basis the primary strain's whole predicted effect is about **0.036**, so 0.02 resolution leaves a factor of **1.8** rather than 2.5. Whether **n ≥ 126** replicates per variant per pH is affordable is a bench and budget question this project cannot answer, and it is flagged here because doubling a replicate requirement can quietly make a committed design unbuildable while every individual number in it stays defensible.

Note that the triple variant adds 11 reversals for only 0.04 pH units beyond the double, because YidC contributes little per reversal (its 122 positions dilute each change, and its across-species SD is the largest in the panel). The **NuoL + NuoN double is therefore the primary strain**: it captures most of the achievable effect for half the construction work, and leaves YidC free for the independent reverse arm.

## 3. What the per-position scan changed about the design

I had intended to copy natural "acidophile-signature" residues — pick positions where acid-tolerant species carry K/R and *E. coli* carries D/E. The scan (`37a`) shows that design is not available:

- Of 250 exposed periplasmic positions, only **4** are polymorphic between acidic and basic at ≥15% each. Most genomes carry **neither** D/E nor K/R at a given position (means: 9% acidic, 5% basic).
- **No single position survives multiple-testing correction** (0 of 165 at q<0.1; 17 at p<0.05 against 8 expected by chance).

So the signal is a **diffuse aggregate** over many surface positions, consistent with the half-split result (both halves positive in 85–95% of random position splits). Consequences:

1. The experiment must change the **aggregate charge**, not reproduce a polymorphism at a few key sites.
2. No individual position is an established causal residue. Per-position ranking only orders the list.
3. Single mutants are expected to be near-undetectable. **A dose-response series is the core of the design, not an add-on.**

## 3b. Chassis: *E. coli* is the weaker choice

Added 2026-10-26 from `scripts/47`, then **weakened the same day by `scripts/48`**. Comparing two monoderm phyla with two diderm ones gave a threefold difference (+0.303 against +0.093, one-sided p 0.0012). But a meta-regression over all ten phyla with at least 12 species reduces the monoderm increment to +0.124 ± 0.079 (one-sided p 0.058), because the two *large* diderm phyla are low while small ones are not — Campylobacterota, a diderm, has the largest per-phylum estimate of all (+0.723). The matched cytoplasmic control shows no envelope contrast either way.

*E. coli* is a diderm, and diderms do show the smaller estimate, so a monoderm chassis — *Bacillus subtilis* is the obvious one, already in the panel as the second SecY reference and with comparable genetics — may show a larger effect. But the expected gain is now modest and uncertain rather than threefold, and it would be wrong to rebuild the design around a marginal moderator.

**Recommendation:** keep *E. coli* for the first build, since everything here is specified for it. A *B. subtilis* arm is worth adding if capacity allows, but on current evidence it is a hedge against a weakly supported moderator, not a reason to switch organisms.

## 4. Protein and direction choice, from where *E. coli* sits

*E. coli*'s own exposed periplasmic acidic excess, against the bacterial distribution:

| Protein | *E. coli* | Bacteria (mean ± SD) | Percentile | Headroom |
|---|---|---|---|---|
| NuoL | +8.2 pts | +3.7 ± 3.7 | 81st | Good for **forward** (11 D/E available) |
| NuoN | +7.0 pts | +7.0 ± 6.5 | 45th | Forward (6 D/E), largest Δ per reversal |
| YidC | +0.8 pts | +4.8 ± 8.2 | 32nd | Already acid-tolerant end; better for **reverse** (14 K/R) |

**Decisions.** Forward arm in **NuoL and NuoN**: non-essential aerobically on glucose, so they can be edited on the chromosome cleanly, and *E. coli* has room to move. Reverse arm in **YidC**: it already sits at the acid-tolerant end with 14 basic positions to make acidic, and YidC is essential, so this arm runs by plasmid complementation in a conditional-depletion background (strain chosen from the published *yidC*-depletion literature). Running the forward arm first is the lower-risk path.

## 4b. The target choice does not survive the reference-structure check (added 2026-10-27)

`scripts/54` ranks the 32 wide-panel references by their own association with pH_min, correcting each for its own measurement noise. Two results bear directly on Section 4, and neither was available when it was written.

**A single reference's estimate is mostly a measure of annotation quality** (`scripts/54`). Across the 31 references with an estimable reliability, the raw per-reference coefficient correlates **+0.817** with that reference's reliability (its item-rest correlation against the rest of the panel). So the NuoL +0.076 and NuoN +0.088 figures — both from the original 7-reference panel, whose per-reference run is not stored in this dataset, so neither is verifiable by `scripts/69` here and +0.088's apparent pass is a chance match at three decimals — that Section 4 and the effect-size scaling in Section 2 rest on cannot be read as "these proteins matter most". They are consistent with NuoL and NuoN simply being among the cleanest proteins to score in the original 7-reference panel.

**The family's estimate depends on which reference structure scores it — suggestively, not conclusively.** PF00361 is represented in the rule-selected panel by HyfB (hydrogenase-4 component B, P23482) instead of NuoL or NuoN. `scripts/55` puts all three on the same 2,286 genomes with the same model: NuoL **+0.048** (p 0.024), NuoN **+0.049** (p 0.021), HyfB **−0.007** (p 0.72), with the NuoL−HyfB and NuoN−HyfB differences at p 0.055 and p 0.051. Suggestive, not established — an earlier note here called this "zero against clearly positive", which overstated it and is corrected.

**The demonstrated problem is worse than that one, and more general** (`scripts/55`). The three reference structures barely agree on the *measurement*. Scoring the same protein in the same genome against different reference topologies gives correlations of just **+0.252** (NuoL vs NuoN), **+0.079** (NuoL vs HyfB) and **+0.087** (NuoN vs HyfB); SecY *E. coli* against SecY *B. subtilis* gives **+0.126**. Median across all pairs **+0.126**. The pipeline itself is sound — the same accession scored through both panels returns r = 1.000 and 0.999 — so this is about the structures. Structural annotation does not transfer between reference structures of one family, and a per-genome exposed-periplasmic measure is largely specific to the reference it was scored against. The earlier reassurance on this point (two SecY species, difference p 0.42) compared coefficients, which agree, and not measurements, which do not.

**What this does and does not change.** It does not touch the comparative finding, which is an aggregate over references and is unaffected by which structure represents any one family. It does undermine the specific choice of NuoL and NuoN as the forward arm, which was made on per-protein estimates that are now known to track scoring quality. Three options, in order of cost:

1. ~~Score PF00361 from NuoL and NuoN within the wide panel.~~ **Run** (`scripts/55`); the result is above. NuoL and NuoN do carry a positive coefficient on the common sample, so the forward arm is not disqualified — but the three references' measurements barely correlate, so no per-protein estimate here is reliable enough to choose a target on.
2. **Move the forward arm to a reference that ranks high with reliability above the floor** — the leaders are an aminoglycoside efflux pump (PF00873, reliability 0.449, +0.243), AmtB (PF00909, 0.416, +0.207) and a multidrug ABC permease (PF12698). AmtB is the most tractable of the three: non-essential, well characterised, and a single polypeptide.
3. **Accept the target as a reporter rather than a maximiser**, i.e. keep NuoL/NuoN and treat a null as uninformative about the hypothesis. This is the weakest option and it makes the declared asymmetry in Section 9 worse, not better.

Recommendation, revised after option 1 was run: **NuoL and NuoN remain defensible as the forward arm** — they are positive and significant on the common sample, and they are the two references whose measurements agree with each other best of any pair tested (+0.252, against a median of +0.126; `scripts/55`). But the choice can no longer be justified by their per-protein effect sizes, and the honest statement in Section 9 is now stronger than before: because measurements are largely reference-specific, a null in one protein bounds that protein scored against that structure and very little else. ~~Adding **AmtB** (PF00909, reliability 0.416) as a second forward arm is the cheapest way to stop the whole experiment resting on one family's reference structures.~~ **Withdrawn the same day it was made** (`scripts/71`, correction 13). AmtB has 40 exposed periplasmic positions, only **5** of which carry D or E in the reference; applying this design's own criterion from `scripts/37b` — a basic residue tolerated in ≥2% of bacteria — leaves **2**, and both have a coefficient at or below zero, the opposite of the direction the arm would test. The primary arm reverses 14 positions; an arm of 2 pointing the wrong way is not an arm. The recommendation was made on panel-level reliability and coefficient without checking the one property that decides whether an arm can be built, and the filter that would have caught it needs only the class table — no alignment, no scan, no compute — and ranks AmtB 7th of 32.

**The second forward arm is P24177** (probable aminoglycoside efflux pump, PF00873), scanned rather than recommended on the filter alone (`scripts/71`): 262 exposed periplasmic positions, **45** carrying D/E and **39** eligible — against 27 for YidC, NuoL and NuoN combined — with 21 of 39 positive and the top three at +0.112 (p 0.002), +0.067 (p 0.009) and +0.053 (p 0.020). It also has the panel's highest reliability at 0.449 — **but that must not be cited as independent evidence for it** (`scripts/72`): P24177 is the largest reference in the panel by a factor of three, and on a measurement-averaging model a reference of 262 positions is expected to score **0.479** on reliability from size alone, slightly above what it actually scores. The argument for this arm is its acidic density, which ranks 2nd of 32 at 0.172 against a panel median of 0.097 and is independent of size; every reference with comparable density has 17 to 60 positions, so P24177 is the only one both dense and large enough for a workable arm.

**One cost of it that no comparative data can settle.** P24177's eligible positions sit in a large folded periplasmic porter domain, where NuoL and NuoN present short interhelical loops. Whether a functional domain tolerates 39 charge reversals is a biochemical question, and the interpretability gate should be expected to reject a larger fraction of P24177 variants than of the primary arm's. Budget for that attrition rather than discovering it. On buildability it is a better arm than the primary one.

Two things to carry with that (`scripts/71`). Reliability and D/E count correlate **+0.591** across the panel, so panel merit is a moderate guide to buildability rather than none — the error was skipping a free filter, not using the wrong criterion. And the mean coefficient across P24177's 39 eligible positions is **−0.0023**, essentially zero, so choosing positions by their own coefficient is selecting largely on noise; that does not inflate the predicted effect, because arms are sized from the compartment-level composite and not by summing per-position coefficients, but the choice of positions carries no evidential weight of its own.

## 4c. The chassis sits in the live half of the range (added 2026-10-28)

Every predicted effect in this document is a number of charge reversals multiplied by the **average** slope. The association is not linear — `scripts/51` found it one-sided — so the average slope is the wrong multiplier for one specific strain, and the strain here is *E. coli*. `scripts/63` settles where *E. coli* sits.

A hinge model (`scripts/63`, with the break's stability from `scripts/64`), `y = a + b_low·x + b_extra·max(0, x−c)` with the headline's covariates and the break `c` profiled by PGLS likelihood, beats both the line and the quadratic on AIC (**6973.2** against 6982.0 and 6977.2):

| | Break (SD of composite) | Slope below | Slope above | Change |
|---|---|---|---|---|
| Hinge, 2,768 bacteria | +0.217 point estimate — **not a usable number, see below** | **≈ +0.29/SD** | **≈ +0.04/SD** | −0.246 ± 0.069 |

**The break's location cannot be pinned down, and an earlier version of this section quoted it as though it could.** `scripts/64` tests it: across 200 random half-samples of the 428 families the fitted break runs from −1.32 to +1.22, with a 10th-to-90th-percentile span of +0.11 to +1.22, so the profile interval of [−0.147, +0.548] understates the uncertainty badly. The hinge is preferred on AIC in 5 of 6 specifications but **not** once the growth optimum is controlled (ΔAIC +2.54), so some of the curvature may belong to the optimum rather than to the composite. Nothing in this design should depend on where the break is.

**What is robust is the ordinal claim, and that is all this section needs.** All 32 wide-panel references are *E. coli* proteins, so *E. coli*'s composite is not an alignment estimate but the reference value itself — the cleanest measurement available for any organism. It is **−0.274 per SD, the 37th percentile**, and it falls below the fitted break **in all six specifications, below the whole profile interval in three of them, and in 198 of 200 half-samples (99.0%)**. The chassis is in the part of the range where the comparative relationship exists, and the forward arm — which lowers exposed periplasmic acidic excess — moves it further in rather than out. This was the risk worth checking and the answer is the favourable one: had *E. coli* sat above the break, the first reversals would have been made where the comparative data show no relationship at all.

**The effect sizes in Section 2 are deliberately left where they are.** The local slope is 1.55× the average — the two rescaled figures below are derived here from it, like Section 2's own predictions, so no result file holds them — which which would turn the 32-reference prediction of −0.05 into −0.078 and the 7-reference prediction of −0.16 into −0.248, and the assay resolution target of 0.02 units would relax accordingly. Three reasons not to take it. The hinge was fitted on the same data that suggested the one-sidedness, so the local slope is a selected quantity whose sampling distribution is not the one the script prints. A factor of 1.55 is small beside the uncertainty already carried here, where the 32-reference and 7-reference scalings differ by more than threefold and the comparative data do not choose between them (Section 2). And raising a predicted effect is the direction that flatters the experiment, which in this project has been the reliable sign that a number needs deriving a second way before it is used — four of the ten recorded corrections went that way. **The conclusion to carry forward is qualitative: the chassis is well placed, and the sizing target stands at 0.02 units.**

## 5. Target positions

25 positions pass selection (`scripts/37b`): reference carries D/E; a basic residue is tolerated at that position in ≥2% of bacteria; not a documented loss-of-function or activity-reducing site (UniProt Mutagenesis features; this excluded YidC 334–336); and for YidC restricted to residues 24–264, the region reported removable without lethality, so a charge change there tests electrostatics rather than folding.

Top-ranked (`scripts/37b`; β = pH_min association of the acidic state, genome-wide composition controlled; **none individually significant after correction**):

| Protein | Position | *E. coli* | % acidic | % basic | β | p |
|---|---|---|---|---|---|---|
| NuoL | 72 | D | 22% | 3% | +0.080 | 0.039 |
| NuoN | 62 | D | 16% | 2% | +0.078 | 0.046 |
| NuoN | 320 | E | 23% | 4% | +0.077 | 0.027 |
| YidC | 90 | E | 12% | 6% | +0.068 | 0.149 |
| NuoL | 57 | E | 23% | 4% | +0.068 | 0.053 |
| NuoN | 195 | D | 28% | 2% | +0.061 | 0.068 |
| YidC | 240 | D | 37% | 9% | +0.054 | 0.073 |

Counts by protein: NuoL 9, NuoN 5, YidC 11.

## 6. Variant panel

**Primary strain: the NuoL + NuoN double**, all 14 selected positions reversed (predicted −0.16 pH units). This is the only variant whose predicted effect sits clearly above the assay floor, so it carries the test. It is built and assayed first; if it shows nothing, the dose series below is not worth constructing.

**Dose series (secondary, for monotonicity).** Reversals added in the rank order of Section 5, capped at the positions each protein actually has: NuoL 1, 3, 5, 7, 9; NuoN 1, 3, 5; combined 7 and 14. Predicted ΔpH_min runs from −0.011 at one reversal to −0.16 at fourteen, so **the bottom half of the series is below the detection floor by construction** and is included only to establish the shape of the response, not to be individually significant. Prediction: monotonic decrease, with a slope per reversal matching Section 2 (−0.011 to −0.012 per reversal for the Nuo proteins).

**Controls.** Each is chosen so the hypothesis predicts a specific, different outcome — that is the point.

| Control | Construction | Prediction if hypothesis true |
|---|---|---|
| Charge-neutral | D→N, E→Q at the same positions (isosteric, removes charge, adds none) | **79% of the reversal's shift**, not merely "intermediate" — see 6b |
| Same-charge | D↔E at the same positions | **Null** — rules out "any mutation at these sites" |
| Matched cytoplasmic | D/E→K at exposed *cytoplasmic* positions (17 available) | **Null** — this class is null in the comparative data |
| Matched buried periplasmic | D/E→K at buried periplasmic positions (22 available) | **Null**, and doubles as a folding control |
| Reconstructed wild type | Same editing procedure, no coding change | Null — controls for construction artefacts |
| Silent codon variants | Synonymous changes at target codons | Null — controls for mRNA/translation effects |

**Interpretability gate.** Any variant with reduced abundance, mislocalisation, or reduced function *at neutral pH* is excluded from the primary analysis and reported separately. Without this gate a destabilised protein masquerades as an acid-tolerance effect. Requires, per variant: immunoblot abundance, membrane localisation, and activity at pH 7 (NADH:ubiquinone oxidoreductase in inverted membrane vesicles for Nuo; insertase activity against a model substrate for YidC).

## 6b. The charge-neutral control has a sharp prediction, and the two panels disagree about it (added 2026-10-31)

This section exists because the reason for reversing D/E to K rather than neutralising them was taken from the 7-reference hand-picked panel, and `scripts/36c` had never been run on the 32-reference rule-selected panel that carries the headline. Run on it (`scripts/36c`, `scripts/36c wide`; synthesis correction 14), with both scales shown because the figures this section first quoted were mislabelled (synthesis correction 15): `pgls_scan` standardises the predictor it tests, so `beta_per_pt` is per SD of that component's share and the two components' shares have different SDs. The per-point columns are from `scripts/74`, and the comparison after the table needs those rather than the per-SD ones.

| Component | Per SD of its share: 7 | Per SD: 32 | **Per percentage point: 7** | **Per point: 32** |
|---|---|---|---|---|
| D+E | +0.041, 6/6 | +0.040 (p 4.1e-10), 30/32 | +0.0083 | **+0.0086 (p 3.4e-10)** |
| Lysine | −0.035, 1/6 positive | −0.010, 11/32 positive | −0.0136 (p 0.0056) | **−0.0024 (p 0.10, not significant)** |

The acidic half of the rationale is confirmed and sharpened. **The basic half is not** — per point lysine is not significant on the rule-selected panel at all. A reversal is worth its extra cost over removing the acid alone only if the basic residue contributes independently, and flipping one of N exposed periplasmic positions changes the D+E share by −100/N points and the K share by +100/N, so the predicted shift is (100/N)(β_K − β_D+E) for a reversal against (100/N)(−β_D+E) for a neutralisation — N cancels, and the comparison needs no assay parameter and no absolute scaling. The margin reversal buys falls from **2.62× to 1.27×**.

**This is a prediction, not a redesign.** Two things follow and a third does not:

1. **The charge-neutral control becomes the arbiter.** Its expected position was stated in Section 6 only as "intermediate". The two panels make it sharp and they disagree: neutralisation alone should recover **38%** of the reversal's shift on the hand-picked panel and **79%** on the rule-selected one — a **40-percentage-point** gap in exactly what that control measures. Recorded here before any strain exists, so it cannot be fitted to whatever the control turns out to do.
2. **It is not the primary endpoint's question.** Separating 38% from 79% of a shift itself predicted at −0.05 pH units needs resolution the confirmation assay does not have. This is a pooled-assay prediction, and only if Gate 2's fitness-per-pH-unit slope is ever measured (`docs/PREREGISTRATION.md` Section 4). It does not alter the committed effect size, the primary test, or the required *n*, all of which rest on the acidic term and are unchanged.
3. **It is not an argument for switching to neutralisation.** Whether an uncharged substitution is easier to build, better tolerated, or more cleanly interpretable than K is a biochemical question the comparative data cannot weigh, and a smaller predicted margin is no reason to prefer a harder readout. The substitution stands; only its quoted justification changes.

**A related weakness in Section 5's selection, which this arm escapes by luck.** The acidic signal is carried mostly by **glutamate**, and `37b` admits any position where the reference carries D *or* E, so the criterion cannot distinguish the stronger component from the weaker. The primary arm happens to be **9 glutamate and 5 aspartate** (the aspartate positions being NuoL 72, 160 and 210 and NuoN 62 and 195), which is the favourable direction, but by accident of ranking rather than by design, since 37b ranks on each position's own association and never on its residue.

**What those five positions cost the arm is now measured** (`scripts/75 wide`, with `scripts/75` for the hand-picked comparison), by entering D and E simultaneously rather than one at a time — the two shares are nearly independent across genomes (median r −0.143, max VIF 1.16), so the joint fit splits a real distinction. Per percentage point, aspartate pools at +0.0050 (p 1.7e-4, positive in 23 of 32) against glutamate's +0.0119 (p 2.0e-11, 32 of 32), a difference of +0.0069 at p 0.002. Glutamate therefore supplies **81%** of the arm's expected effect and the five aspartate positions **19%**: the arm is worth about **11.1 glutamate-equivalent positions**, a 21% shortfall against counting all 14 equally. On the hand-picked panel the two residues cannot be told apart (p 0.56) and the shortfall is 9% — the pattern a genuine panel difference should show, and which the earlier cross-model version, below 1 on both panels alike, did not. Note that this also retires the sharper claim that aspartate is "a coin flip": its univariate 18-of-32 was attenuated by omitting glutamate, with which it is slightly negatively correlated.

**A consequence for sizing, which is not applied here.** If a position contributes in proportion to its own residue's coefficient, the primary arm's predicted effect is 0.79 of the figure Section 2 gives, and the pre-registered *n* would rise from 64 to 101 per group (`scripts/75 wide`). That proportionality has not been tested by anything, so `docs/PREREGISTRATION.md` Section 3 stands as committed and any change goes through its Section 10. **Recorded now so it is a decision rather than a discovery at the bench** — and it is a reason to prefer glutamate positions if the arm is ever re-selected, which `37b` currently has no way to express.

## 7. Primary readout: pooled barcoded competition

Individual growth curves cannot resolve a 0.16 pH-unit shift in a growth limit, let alone the 0.02 units the dose series needs. Instead: build the panel as a **barcoded pool**, compete it in parallel at several pH values, and read relative fitness by amplicon sequencing. This resolves fitness differences of ~1%, is internally controlled (all variants share the culture), and removes plate-position and batch effects.

**Sensitivity must be demonstrated, not assumed** — and `scripts/70` makes precise what is missing: the pooled readout is log relative fitness, and converting it to pH_min needs a fitness-per-pH-unit slope that has not been measured for this chassis or any other, so the resolution claim below is an assumption until the pilot in `docs/PREREGISTRATION.md` Gate 2 supplies that slope with a standard error. Because the primary prediction is now −0.16 and the per-reversal prediction −0.011, the assay's resolution has to be established empirically before the variants are interpreted: include a spike-in ladder of strains with independently measured small fitness differences, and report the smallest pH shift the pool can distinguish from zero. If that floor is above ~0.05 units, the dose series is uninterpretable and only the combined variant should be analysed.

- Buffered defined medium across pH 3.8–5.0 in 0.2-unit steps, plus pH 7.0 as reference.
- **Both acid conditions:** strong acid (HCl-adjusted) and a weak organic acid at fixed undissociated concentration. This project's schema tracks acid identity because the undissociated fraction matters, and the hypothesis is about surface charge, not uncoupling — the two conditions separate them.
- Primary statistic: log fitness of each variant relative to reconstructed wild type, at challenge pH versus pH 7.0.
- Serial passage with bottleneck control; ≥6 independent pool replicates; barcode counts ≥10⁴ per variant per timepoint.

**Confirmation.** The combined variant and the top of each dose series re-measured as isolated strains with full pH–growth-rate curves, pH_min taken from a fitted continuous curve rather than a visual threshold. Powering for the revised −0.16 effect: at SD 0.1, α 0.05 and 80% power, detecting 0.16 units needs n≈8 per group and detecting 0.10 units needs n≈17, so n≥20 biological replicates per variant per pH was the working floor. **That floor is superseded** (`scripts/70`, `docs/PREREGISTRATION.md`): it was set against −0.16, and against the conservative −0.05 from the 32-reference scaling it detects only 0.091 pH units, underpowered by a factor of 1.8. The committed requirement is **n ≥ 64 per group** for the primary endpoint, or 119 under correction across the 17 declared arms. The assay SD of 0.1 assumed throughout has no source and is the parameter that decides feasibility — at SD 0.05 the conservative case needs n=17, at SD 0.2 it needs n=253 — so measuring it is a gate before the panel is built — more than double the earlier figure, which was set against the inflated effect sizes. Genotype labels blinded during scoring; plate positions randomised.

## 8. Mechanistic chain (secondary)

If growth-limit effects appear, establish *why* — surface charge should act through envelope-proximal physiology, not transcription:

1. **Cytoplasmic pH homeostasis** under acid challenge — ratiometric pH-sensitive GFP; prediction: reversal variants maintain cytoplasmic pH better at low external pH.
2. **Proton-motive force / membrane potential** — prediction: smaller collapse at low pH.
3. **Protein-level stability at low pH** — abundance and activity at pH 7 versus pH 4; prediction: reversal variants retain function at low pH.
4. **Specificity** — the cytoplasmic-matched control should not move any of these.

## 9. What would refute the hypothesis

Stated in advance:

- **Refuted** if the NuoL + NuoN combined variant shows no low-pH fitness advantage, while its neutral-pH abundance and function are intact and the assay has been *demonstrated* (spike-in ladder, Section 7) to resolve a shift of 0.05 pH units or smaller. Without that demonstrated floor a null is uninformative, because the predicted effect is −0.16.
- **Supported** if the combined variant shifts in the predicted direction, the dose-response is monotonic over the range where the assay is sensitive, the effect is absent in cytoplasmic- and buried-matched controls, intermediate in charge-neutral controls, and reversed in the YidC K/R→D/E arm.
- **Quantitative check, not just a sign test.** The predicted magnitudes are now specific — −0.011 per reversal and −0.16 for the combined variant — so an effect that appears but is an order of magnitude larger than predicted argues for a different mechanism (most likely folding or expression damage that the interpretability gate missed), not for the hypothesis.
- **Honest asymmetry:** a null does **not** fully refute the comparative finding. The signal is an aggregate over many proteins; no single protein may be sufficient. A null bounds the per-protein contribution, which is itself worth knowing, but the aggregate model survives it. This asymmetry is a real weakness of the design and should be declared before data collection, not after.

## 10. Parallel computational test — RUN, and it does not discriminate

Run on 2026-10-25 (`scripts/38_substitution_enrichment.py`). Two findings: the test as originally specified is not possible, and the substitute test turns out to have almost no power — which is itself the most useful thing it established.

**The lab-evolution version is not runnable.** ALEdb exposes no API or bulk download (browser CSV export only) and documents no acid or low-pH experiment. The published *E. coli* acid-adaptation experiments report only a few missense mutations per clone, concentrated in regulators (*rpoC*, *rpoB*, *evgS*, *cadA*) rather than in conserved membrane proteins, so the compartment contrast has no power there.

**Substitute: sister pairs on the bacterial tree** — 844 pairs, 600 differing in pH_min, most within a genus, so differences reflect recent change rather than standing correlation.

- *Compartment level:* the exposed periplasmic composite is higher in the less acid-tolerant sister in 52.0% of pairs (p 0.14), rising to 58.3% of the 133 pairs differing by ≥1 pH unit (p 0.013 uncorrected). **Nothing survives Benjamini-Hochberg across the 18 compartment tests** (0 at q<0.1), and the transmembrane control is nominally significant too, which the hypothesis does not predict.
- *Residue level:* among 5,255 substitutions where the sisters differ and exactly one residue is acidic, the acidic residue sits in the less acid-tolerant member 50.8% of the time (p 0.13), with no compartment specificity — buried periplasmic is nominally the strongest at 51.5%.

**Why this is not evidence against the hypothesis** (`scripts/38`). A slope of pH on composition does not imply a large composition difference for a given pH difference. The expected sister-pair difference is the *reverse* regression coefficient, `slope × E[Δx²]/E[Δy²]`, which is tiny because the predictor explains little variance: only **0.014 SD of composite per pH unit** at the between-species slope. Against that expectation the observations are almost exact:

| Threshold | Pairs | Observed | Expected (between-species) | Deviation |
|---|---|---|---|---|
| any difference | 573 | +0.0136 ± 0.0159 SD | +0.0139 | 0.0 SE |
| ≥0.5 units | 371 | +0.0073 ± 0.0184 SD | +0.0182 | 0.6 SE |
| ≥1.0 units | 133 | +0.0697 ± 0.0324 SD | +0.0269 | 1.3 SE |

So the sister-pair null is exactly what the comparative estimate predicts. The test cannot confirm or refute: it was never powered to.

**What it establishes.** No observational test on recent natural substitutions can settle this — the expected signal is an order of magnitude below the noise — which makes the intervention in Sections 1–9 necessary rather than optional.

**Correction (2026-10-25, after `scripts/39a`).** This section first reported that the calibration "favours the between-species slope over the within-family slope", implying the within-family estimate was inflated. That was wrong, and the error was mine: the two slopes were on different scales (per SD against per unit of a composite whose SD is 0.53), so they were never comparable. On a common scale they agree, and the sister-pair observations are consistent with **all** of them — 0.8, 1.5 and 0.5 SE from the +0.107 expectation at the three thresholds, and 0.5 to 1.2 SE from the +0.088 family-fixed-effects expectation. The calibration therefore discriminates between nothing; its only conclusion is the power one above.

**Section 2 was rewritten on this basis**: reversal counts capped at the positions each protein has, and the primary strain is the NuoL + NuoN double at −0.16 pH units rather than any single-protein variant. Sections 6, 7, 9 and 11 were brought into line — a demonstrated sensitivity floor, n≥20 replicates instead of n≥8, and a dose series whose lower half is acknowledged to sit below the detection limit. Those revisions stand, because they were derived from the +0.101 slope, which is the one the correction confirms.

## 11. Order of work

1. ~~Computational enrichment test~~ — **done (Section 10); it could not discriminate, and the design was revised on its calibration.**
2. **Establish the assay floor first.** Spike-in ladder in the pooled competition, to show the smallest resolvable pH shift. If that floor is worse than ~0.05 units the rest of the plan does not work as written and needs a more sensitive readout before any variant is built. This is now the first bench step, because the revised effect sizes are small enough that an unmeasured floor would make every later result uninterpretable.
3. **Primary strain:** NuoL + NuoN combined variant, all 14 positions, chromosomal, with the interpretability gate and the matched cytoplasmic, buried, charge-neutral and same-charge controls in the same pool.
4. Dose series, for the shape of the response over the sensitive range only.
5. Confirmation as isolated strains, n≥20 per variant per pH.
6. Mechanistic chain, if an effect appears.
7. Reverse arm in YidC by complementation — independent sign test, and the only arm where YidC is worth the construction effort.
