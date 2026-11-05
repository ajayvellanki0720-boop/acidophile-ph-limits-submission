# Published experimental evidence, against this project's two residue predictions

Assembled 2026-11-01. The question asked was the one a comparative study can answer without a bench: has anyone already **experimentally manipulated** periplasmic acidic or histidine residues in bacterial proteins and measured something about acid resistance, and does the result agree with what this project predicts?

**Both predictions have direct experimental support. Both arms' compartment is confirmed. One arm's DIRECTION is contradicted.** The histidine prediction is supported from four systems, three of them periplasmic. The acidic prediction has one directly comparable experiment — an exhaustive mutagenesis of every titratable residue in a bacterial membrane protein's periplasmic domain — and in it periplasmic carboxylates are **required** for function at low pH, where this project finds them depleted in acid-tolerant genomes. The two are reconcilable and the reconciliation is specific; §2.2 makes it rather than glossing it.

Search route: Consensus semantic search, then an exhaustive PubMed query (59 records, all screened). **The first version of this document ran only the semantic search, which is capped at three results per query, and concluded that the acidic arm had no direct experiment at all.** That was wrong, and the single most relevant paper in this document is what the capped search missed. Absence below means "not found by a 59-record PubMed screen plus semantic search", which is a stronger statement than before but still not a systematic review.

---

## 1. The histidine prediction — supported, repeatedly and independently

This project finds exposed **periplasmic histidine** associated with a lower growth limit, with a side-specific component (−0.080003 pH per percentage point, p 3.1e−03) and solvent accessibility doing demonstrable work: histidine at non-exposed positions in the same families carries nothing (`scripts/83 wide`). Side-chain pKa ≈ 6.0 puts histidine's titration inside the physiological range, which is the mechanistic reason to expect it to matter.

| Study | System | What was manipulated | Result | Bearing |
|---|---|---|---|---|
| Sen et al. 2017, *J Bacteriol* ([link](https://consensus.app/papers/details/dbbd722c16c251bbb6b24a7408910cb9/)) | *E. coli* EvgS, an acid-sensing histidine kinase with a large periplasmic domain | Mutagenesis of modelled periplasmic-domain residues | Histidine 226 is required for acid induction, and sits next to a proline at the top of the predicted transmembrane helix | **Direct and on the right side.** A single periplasmic histidine is necessary for an acid response in a membrane protein |
| **Weeks & Sachs 2001, *Mol Microbiol*** (PMID 11442825, [DOI](https://doi.org/10.1046/j.1365-2958.2001.02466.x)) | *H. pylori* UreI urea channel, required for acid resistance | **Every** periplasmic histidine, aspartate and glutamate mutated | His-123 and His-131 in periplasmic loop 2 and His-193 in the C-terminus are needed to activate transport at acidic pH; the pH-independent *S. salivarius* homologue lacks them | **The strongest evidence here, and exhaustive rather than selective.** Three periplasmic histidines are necessary for acid activation, and §2 records that the same scan also implicates periplasmic carboxylates |
| Stoib et al. 2026, *Small* ([doi:10.1002/smll.76130](https://doi.org/10.1002/smll.76130)) — peer-reviewed publication of the preprint cited in earlier drafts; retitled on publication | UreI urea channel, *H. pylori* / *H. hepaticus* / *S. salivarius* | Cytoplasmic E136 mutated across homologues | Confirms E136 as a conserved pH sensor, and states that pH-gating in *H. pylori* UreI had been attributed to **periplasmic histidines** | **Supports the histidine side** as prior consensus for that channel, while adding a cytoplasmic acidic sensor the present work's cytoplasmic control does not see |
| Gallo et al. 2025, *J Bacteriol* ([link](https://consensus.app/papers/details/b413e1f9d13c55ff91a58ec874772c59/)) | *E. coli* AdiY, regulator of the arginine-dependent acid resistance system | His34 and His60 replaced by Ala, Gln or Asp | Each replacement abolishes pH-dependent activation in vivo and in vitro | **Mechanism, not sidedness.** Histidine protonation as a pH sensor is established; these histidines are cytoplasmic, so this supports the residue and not the compartment |

**What this does and does not license.** Three independent systems use histidine protonation to sense pH, and in two of them the relevant histidines are periplasmic. That is strong corroboration that exposed histidine is *mechanistically* the right residue to be looking at. It is **not** corroboration of the comparative claim: these are single residues in dedicated sensor proteins, whereas this project measures a compositional share averaged over 32 housekeeping membrane proteins, and nothing above says that share should predict a growth limit. A dedicated pH sensor and a bulk surface property are different hypotheses that happen to involve the same amino acid.

**It does, however, speak to the open question.** `scripts/85 wide` could not identify what the *shared* (non-side-specific) half of the histidine association is shared with, and eleven candidates failed. The literature above suggests a reading that this dataset cannot test: if histidine serves proton-sensing roles on both faces of the membrane — periplasmic in EvgS and UreI, cytoplasmic in AdiY — then a shared component is what one would expect, and the side-specific excess is the part attributable to the outward-facing surface specifically. That is a hypothesis, offered as one, and the project has no way to discriminate it.

---

## 2. The acidic prediction — there IS a direct experiment, and it points the other way

> **Correction, same day.** The first version of this section stated that no study was found manipulating periplasmic acidic residues in a bacterial membrane protein and measuring acid resistance. **That was wrong.** It rested on a semantic search capped at three results per query; an exhaustive PubMed search over the same question returned 59 records, one of which is the most directly relevant experiment in this entire document. The original hedge — "not found by these queries" — was accurate and did not save the conclusion, which was drawn from a search too thin to support it. The superseded text is not retained because it asserted an absence that does not exist.

### 2.1 The experiment that exists

**Weeks & Sachs 2001, *Molecular Microbiology*** — "Sites of pH regulation of the urea channel of *Helicobacter pylori*" (PubMed PMID 11442825, [DOI](https://doi.org/10.1046/j.1365-2958.2001.02466.x)).

UreI is the urea channel that *H. pylori* requires for acid resistance, and it is an inner-membrane protein with defined periplasmic loops. The authors **mutagenised every histidine, aspartate and glutamate, and the lysine, in HpUreI's periplasmic domain** — which is precisely the exhaustive periplasmic titratable-residue scan this project's prediction concerns. The pH-dependence is real and measured: HpUreI's urea uptake is acid-activated, half-maximal at pH 6.0, while the *Streptococcus salivarius* homologue is pH-independent.

What came back as important for acid activation of transport:

| Residue | Location |
|---|---|
| His-123, His-131 | periplasmic loop 2 |
| **Asp-129, Asp-140, Glu-138** | periplasmic loop 2 |
| Lys-132 | periplasmic loop 2 |
| His-193 | C-terminus |

Their stated conclusion is that protonation of three histidines plus the *presence* of three dicarboxylic amino acids in periplasmic loop 2 is needed to activate the channel at acidic pH. And the comparative half is built in: those charged residues are largely **absent from the pH-independent SsUreI**, with loop-swap chimeras behaving accordingly.

### 2.2 What it does to this project's account

**It confirms the compartment and the residue classes, and contradicts the direction.** This project finds exposed periplasmic acidic residues **depleted** in acid-tolerant bacteria (`scripts/57`), and offers charge burden as the mechanism: carboxylates protonate as external pH falls, so fewer of them means less pH-dependent change in outward-facing surface charge. In UreI, periplasmic carboxylates are **required** — their presence is what makes the protein work at low pH, and the homologue that lacks them cannot respond to acid at all.

**The two are compatible, and the distinction has to be made explicitly rather than glossed.** A protein whose job is to *respond* to acid needs titratable residues on the face that senses it; a proteome that must merely *function* at low pH may do better with fewer. UreI's residues are a dedicated gate in one specialised channel; this project measures bulk composition averaged over 32 housekeeping membrane proteins and relates it to a growth limit. Those are different quantities and a trend in one does not contradict a mechanism in the other.

**But it is a genuine counterexample to the rationale, and it favours the account this project did not adopt.** The charge-burden story says periplasmic carboxylates are a liability at low pH. Here they are an asset. Taken with the serratiopeptidase engineering study below, which adds acidic residues to survive gastric pH on an explicit buffering argument, the weight of direct experiment sits with **buffering and function** rather than with charge burden — which is the account this project left open for histidine and declined for aspartate and glutamate.

**It is also the best available precedent for the pre-registered experiment.** An exhaustive scan of periplasmic titratable residues in a single membrane protein with a pH-dependent functional readout is, methodologically, what [`CHARGE_REVERSAL_DESIGN.md`](CHARGE_REVERSAL_DESIGN.md) proposes. That this has been done once, successfully, in 2001, is worth more to the design than any power calculation in it.

### 2.3 The enzyme-engineering literature, which measures something else

No other study was found manipulating periplasmic acidic residues in a bacterial membrane protein with an acid phenotype. The remaining work is surface-charge engineering of **soluble enzymes** assayed for **in vitro stability**, and it is directionally mixed:

| Study | What was done | Result | Bearing |
|---|---|---|---|
| Dhiman et al. 2022, *J Biomol Struct Dyn* ([link](https://consensus.app/papers/details/62e6ca22790d594e821099f626e64f98/)) | Serratiopeptidase; non-conserved amide residues converted **to** acidic (N412D) to survive gastric pH | Rationale is explicitly that a more acidic surface **buffers** against stomach pH | **Opposite sign to the charge-burden account**, and agrees with UreI in favouring a functional/buffering reading |
| Zheng et al. 2025, *Int J Biol Macromol* ([link](https://consensus.app/papers/details/f461adcc43da57d58a776bb3545a95f3/)) | κ-carrageenase; D49I removes a surface acidic residue | Improved **alkaline** tolerance (pH 10) | Right manipulation, wrong environment |
| Li et al. 2026, *Int J Biol Macromol* ([link](https://consensus.app/papers/details/fd256208f9a6541d8de6d3e99ddda6e9/)) | *A. flavus* uricase; surface charge optimisation | Acid stability improved (67%/64% activity retained at pH 4.5 vs 25% wild type), attributed to charge **redistribution** | Supports that surface charge matters at low pH, with no directional claim |

All three measure retained catalytic activity rather than an organism's growth limit, and a designed substitution in one enzyme is not a compositional trend across a proteome.

## 3. One methodological precedent, which is also the sharpest specificity threat found

**Cabello-Yeves et al. 2019, *Microbiome*** ([link](https://consensus.app/papers/details/6a7d716b476b578eb36acc901a75f470/)) compared marine and freshwater prokaryotes, including phylogenetically matched pairs, and found predicted isoelectric point shifting with habitat — more basic in freshwater — with the difference **sharper in secreted than cytoplasmic than membrane proteins**, and concentrated on protein surfaces.

That is the same compartment ordering this project reports for pH, from a different group, for a **different environmental variable**, six years earlier. It cuts two ways and both should be in the manuscript:

- **As precedent it is favourable.** Compartment-resolved surface-charge analysis is an established approach with a published result of this shape, and the secreted > cytoplasmic > membrane ordering is not an artefact peculiar to this pipeline.
- **As a threat it is serious.** If outward-facing surface charge tracks habitat generally, then part of what this project attributes to pH may be a broader habitat signal that correlates with pH. The project's matched controls rule out *position-class* artefacts; they cannot rule out a second environmental variable acting through the same compartment.

**It is not testable here, and the reason is a gap in the data rather than a choice.** There is no salinity or habitat variable for the bacterial sample: the `medium` column is empty for all 3,002 phenotype rows, as the `oxygen` column already was. A halophile indicator exists only for the archaea, where `scripts/43` already uses it — and there, notably, specificity *fails* inside the halophiles, with transmembrane turning positive alongside the extracellular class. That is the one place in this project where the compartment pattern breaks, and it is a salinity-defined subset. **This is now the strongest named candidate for the unidentified shared component** and it is recorded in [`SYNTHESIS.md`](SYNTHESIS.md) §8 as such.

### Two further pieces of context

- **Reed et al. 2013, *Archaea*** ([link](https://consensus.app/papers/details/a1cf9255cc82551c870aa74c3a6b03c3/)), 282 citations, reviews archaeal extremophile protein adaptation and states that adaptations toward pH are less discernible than the thermophilic, psychrophilic and halophilic cases. That is the gap this project addresses, and it is useful for the introduction: the absence of a known compositional signature for pH is a documented position, not an assumption made here.
- **Ullrich et al. 2024, *Front Microbiol*** ([link](https://consensus.app/papers/details/a5d6a72f29c05e0a9f37e5f7fe760f70/)) finds periplasmic redox proteins of acidophiles having more positive redox potentials than neutrophilic homologues, achieved mainly through hydrophobic shielding of cofactor pockets rather than surface charge. So acidophile **periplasmic** proteins are demonstrably different from neutrophilic ones — by a mechanism this project does not measure and whose relationship to surface composition is unknown.

---

## 4. What this changes

| | Before this search | After |
|---|---|---|
| Histidine arm | mechanism unsupported externally | **four** independent systems use histidine protonation to sense pH, three of them periplasmic, one of them an exhaustive scan |
| Acidic arm | mechanism assumed plausible, unexamined | **one direct, exhaustive periplasmic experiment exists** — and periplasmic carboxylates come out as required for low-pH function, not as a liability |
| The charge-burden mechanism | the project's stated account | **not the account the direct evidence favours.** UreI and the serratiopeptidase engineering both point to buffering or function. The association stands on its own evidence; the mechanism offered for it does not |
| Shared histidine component | eleven candidates excluded, none from outside the dataset | a literature-motivated hypothesis exists (proton sensing on both membrane faces), untestable here |
| Specificity | position-class controls all null | **a second environmental variable with the same compartment ordering is documented**, and cannot be tested for want of a habitat variable |
| Charge-reversal experiment | one of several open items | **has a 2001 precedent that worked** — an exhaustive periplasmic titratable-residue scan with a pH-dependent functional readout, which is what the design proposes |

**The one-line summary for the manuscript:** the compartment is corroborated twice over, and the direction is not — the one protein where every periplasmic titratable residue has been tested needs its carboxylates to work in acid, while acid-tolerant genomes carry fewer of them. That is the most interesting unresolved thing in this project, and it is now a stated tension rather than an untested assumption.

**What I would want next, and cannot do here.** The reconciliation in §2.2 — sensors need titratable residues, housekeeping proteins do better without them — is testable in principle: the panel's 32 references are housekeeping proteins, and dedicated pH sensors like EvgS, CadC and UreI are not in it. If the account is right, a panel built from sensors should show the *opposite* sign. That requires structures and detection for a sensor panel, which this project does not have, and it is the single most informative follow-up the search produced.
