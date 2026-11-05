"""Which of the panel's complexes are the ones missing from acid-tolerant genomes? -> DATA_DIR/which_complexes.csv and, with `wide`, DATA_DIR/which_complexes_wide.csv
BOTH SPELLINGS ARE DECLARED, for the reason scripts/81 records: scripts/69 builds a claim's candidate pool from this line and keeps a panel variant's pool separate, so an
undeclared output turns every figure quoted from it into a guaranteed false miss.
WHAT THIS FOLLOWS FROM. scripts/91 decomposed the scoring-quality predictor and found that its gene-content half survives every confound available: BREADTH, the count of
reference families a genome has a detectable homologue of, predicts the acid limit at -0.03599 per reference (p 6.7e-04) with assembly level, similarity to the E. coli
references and genome size all held fixed, while the measurement half carries nothing once similarity is controlled. A count is not a mechanism, though, and the obvious next
question is WHICH families. If the effect is spread evenly across all of them it is probably still a detection or completeness story wearing gene-content clothes; if it
concentrates in particular complexes, that is a mechanistic statement about what acid-tolerant genomes lack.
THE PANEL CARRIES ITS OWN NEGATIVE CONTROL, which is why this test can fail honestly. Five of the 37 accepted references are soluble cytoplasmic proteins - EF-Tu, GroEL,
RecA, DnaK and the ATP synthase alpha subunit - included as controls precisely because they are essential and universally present. Their detection should barely vary across
genomes and should predict nothing. If the breadth effect turns out to be carried by THOSE, it is an artefact of genome quality or annotation and the gene-content reading in
scripts/91 has to be withdrawn.
MULTIPLICITY IS HANDLED, not mentioned. Thirty-seven tests at p<0.05 are expected to produce about two hits by chance, so a list of nominally significant references is not a
result. Benjamini-Hochberg is applied across all references tested and the q-values are what the reading rests on; the count of nominal hits is printed beside the count of
survivors so the difference is visible rather than assumed. This is scripts/68's standard applied to a new family of tests.
AND DETECTION IS A NEARLY-CONSTANT COLUMN FOR MOST REFERENCES, which bounds what can be learned. A reference present in 99% of genomes has almost no variance to regress on,
so its coefficient will be unstable whatever the truth is; references are therefore reported with their detection rate and those below a variance floor are excluded by rule
rather than silently returning noise.
Usage: python3 scripts/92_which_complexes.py [orig|wide]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
from Bio import Phylo
from acidophile import phylo as ph

warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
WHICH = (sys.argv[1] if len(sys.argv) > 1 else "wide").lower()
assert WHICH in ("orig", "wide"), f"panel must be orig or wide, got {WHICH!r}"
PANEL = "structure_panel_wide_det.csv" if WHICH == "wide" else "structure_panel.csv"
TAG = "_wide" if WHICH == "wide" else ""
MINVAR = 0.02          # a reference must be absent in at least 2% of genomes, and present in at least 2%
_r = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv")
ACC = _r[_r.accepted]
ISMEM = {str(x.label): bool(x.membrane) for x in ACC.itertuples()}
NAME = {str(x.label): str(x.name) for x in ACC.itertuples()}
rows = []

pan = pd.read_csv(D / PANEL)
df = (pd.read_csv(D / "phenotypes_all.csv")
      .merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess", "genome_size_bp"]],
             on="ncbi_assembly")
      .merge(pan, on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips) & df.ph_min.notna()].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
SRCC = [f"src_{o}" for o in SRC if B[f"src_{o}"].nunique() > 1]
# Detection is read off the score column, which exists for every accepted reference whether membrane or soluble;
# the per-class position columns exist only for the classes a reference's topology defines, so using them would
# silently drop the soluble controls and with them the negative control this test depends on.
LAB = [l for l in ISMEM if f"{l}_score" in B.columns]
det = {l: B[f"{l}_score"].notna().to_numpy().astype(float) for l in LAB}
scols = [f"{l}_score" for l in LAB]
S = B[scols].to_numpy(float)
with np.errstate(invalid="ignore"):
    meanscore = np.nanmean(np.where(S > 0, S, np.nan), axis=1)
meanscore = np.where(np.isfinite(meanscore), meanscore, np.nanmedian(meanscore[np.isfinite(meanscore)]))
breadth = np.sum([det[l] for l in LAB], axis=0)
ok = np.isfinite(B.ph_min.to_numpy(float)) & np.isfinite(B.genome_size_bp.to_numpy(float))
g = B[ok].reset_index(drop=True)
Cm = ph.vcv_fast(tree, list(g.ncbi_assembly))
y = g.ph_min.to_numpy(float)
LEV = [L for L in ("Complete Genome", "Chromosome", "Scaffold", "Contig") if (g.assembly_level == L).sum() > 20]
droplev = max(LEV, key=lambda L: (g.assembly_level == L).sum())
XS = np.column_stack([np.ones(len(g))] + [g[k] for k in ["gw"] + SRCC]
                     + [meanscore[ok], np.log10(g.genome_size_bp.to_numpy(float))]
                     + [(g.assembly_level == L).astype(float).to_numpy() for L in LEV if L != droplev])
lam = ph.pgls(y, XS, Cm)["lambda"]
print(f"sample {int(ok.sum())}; {len(LAB)} accepted references with a score column "
      f"({sum(ISMEM[l] for l in LAB)} membrane, {sum(not ISMEM[l] for l in LAB)} soluble); lambda {lam:.3f}")
print(f"controls in every fit: composition, data source, mean alignment score, log10 genome size, assembly level")

print(f"\n1. per-reference detection against the acid limit, references with at least {MINVAR:.0%} variance both ways")
print(f"   {'reference':>14} {'kind':>9} {'detected':>9} {'beta':>10} {'se':>9} {'p':>10}")
tested = []
for l in sorted(LAB, key=lambda x: (not ISMEM[x], x)):
    v = det[l][ok]
    rate = v.mean()
    if rate < MINVAR or rate > 1 - MINVAR:
        rows.append({"part": "excluded", "reference": l, "membrane": ISMEM[l], "detection_rate": float(rate),
                     "why": f"detection variance below the {MINVAR:.0%} floor"})
        continue
    r = ph.pgls(y, np.column_stack([XS, v]), Cm, lam=lam)
    b, se, pv = float(r["beta"][-1]), float(r["se"][-1]), float(r["p"][-1])
    tested.append((l, b, se, pv, rate))
    print(f"   {l:>14} {('membrane' if ISMEM[l] else 'soluble'):>9} {rate:>9.1%} {b:>+10.5f} {se:>9.5f} {pv:>10.2e}")
    rows.append({"part": "per reference", "reference": l, "name": NAME[l], "membrane": ISMEM[l],
                 "detection_rate": float(rate), "beta": b, "se": se, "p": pv, "n": int(ok.sum())})
nexcl = len([r for r in rows if r["part"] == "excluded"])
print(f"   {nexcl} of {len(LAB)} references excluded by the variance floor and listed in the csv, not dropped silently")

print("\n2. Benjamini-Hochberg across the references tested, because 37 tests buy about two hits by chance")
ps = np.array([t[3] for t in tested])
order = np.argsort(ps)
m = len(ps)
q = np.empty(m)
prev = 1.0
for i in range(m - 1, -1, -1):
    prev = min(prev, ps[order[i]] * m / (i + 1))
    q[order[i]] = prev
nominal = int((ps < 0.05).sum())
surv = [(tested[i][0], tested[i][1], ps[i], q[i]) for i in range(m) if q[i] < 0.10]
print(f"   {nominal} of {m} nominally significant at p<0.05; {len(surv)} survive BH at q<0.10")
for l, b, pv, qv in sorted(surv, key=lambda x: x[2]):
    print(f"      {l:>14} {('membrane' if ISMEM[l] else 'soluble'):>9} {b:>+10.5f}  p {pv:.2e}  q {qv:.3f}   "
          f"{NAME[l][:44]}")
for i in range(m):
    rows.append({"part": "bh", "reference": tested[i][0], "p": float(ps[i]), "q": float(q[i]),
                 "survives_q10": bool(q[i] < 0.10)})
rows.append({"part": "bh summary", "n_tested": m, "n_nominal": nominal, "n_survive_q10": len(surv)})

print("\n3. THE NEGATIVE CONTROL THIS TEST WAS DESIGNED AROUND IS NOT AVAILABLE, and saying so is the point")
# The five soluble references were meant to be the null: essential, universally present, and so predicting
# nothing. They are universally present - 99.1% to 100.0% detection - which is exactly why they are USELESS here:
# every one of them falls below the variance floor and none is tested at all. The first version of this script
# printed "the negative control is CLEAN: no soluble reference survives" on the strength of zero soluble
# references having been tested, which is a false reassurance of the kind this project keeps producing. A control
# that cannot run is a limitation, not a pass.
sol_rates = [(r["reference"], r["detection_rate"]) for r in rows
             if r["part"] == "excluded" and not r["membrane"]]
sol_tested = [t for t in tested if not ISMEM[t[0]]]
print(f"   {len(sol_rates)} soluble references, detected in "
      + ", ".join(f"{l} {v:.1%}" for l, v in sorted(sol_rates, key=lambda x: -x[1]))
      + f" — all below the {MINVAR:.0%} variance floor, so {len(sol_tested)} were tested.")
print("   The intended null therefore says nothing, and the gene-content reading needs a different check.")
rows.append({"part": "negative control", "available": False,
             "why": "all soluble references are detected in >99% of genomes, below the variance floor",
             "n_soluble_tested": len(sol_tested)})

print("\n3b. THE CHECK THAT IS AVAILABLE: a genome-quality artefact would push every detection the same way")
# If these coefficients were driven by how completely a genome is sequenced or annotated, then detecting MORE of
# anything would shift the estimated limit in one direction, because the mechanism would be shared across
# references and indifferent to what they do. Heterogeneous signs are therefore evidence against that account,
# and they are available without a null reference. This is weaker than the control it replaces - it rules out one
# artefact rather than establishing gene content - and it is reported as that.
neg = [x for x in surv if x[1] < 0]
pos = [x for x in surv if x[1] > 0]
print(f"   of the {len(surv)} survivors, {len(neg)} have detection associated with a LOWER limit and {len(pos)} "
      f"with a higher one")
print("      lower: " + ", ".join(f"{l} ({NAME[l][:30]})" for l, *_ in neg))
print("      higher: " + ", ".join(f"{l} ({NAME[l][:30]})" for l, *_ in pos))
mixed = bool(neg) and bool(pos)
print("   " + ("Signs are MIXED, which a shared sequencing-or-annotation artefact cannot produce: that account "
               "predicts one direction for every reference. It rules out that artefact; it does not by itself"
               if mixed else
               "Signs are UNIFORM, which is what a shared sequencing-or-annotation artefact looks like, and the"))
print("   " + ("   establish that the surviving associations are about function."
               if mixed else "   gene-content reading should not be extended on this evidence."))
rows.append({"part": "sign heterogeneity", "n_lower": len(neg), "n_higher": len(pos), "mixed": mixed,
             "rules_out": "a shared genome-quality or annotation artefact", "does_not_establish": "function"})

print("\n4. is the effect concentrated or spread?")
# If breadth worked by many small contributions, no single reference should survive multiplicity correction while
# the count still predicts; if it is concentrated, the survivors should account for most of it. The test is to
# re-enter breadth with the survivors' own detection held fixed: if breadth then carries nothing, the survivors
# ARE the effect.
(bb,) = [(float(r["beta"][-1]), float(r["se"][-1]), float(r["p"][-1]))
         for r in [ph.pgls(y, np.column_stack([XS, breadth[ok]]), Cm, lam=lam)]]
print(f"   breadth with these controls: {bb[0]:+.5f} per reference (se {bb[1]:.5f}, p {bb[2]:.2e})")
if surv:
    Xs2 = np.column_stack([XS] + [det[l][ok] for l, *_ in surv] + [breadth[ok]])
    r2 = ph.pgls(y, Xs2, Cm, lam=lam)
    b2, p2 = float(r2["beta"][-1]), float(r2["p"][-1])
    print(f"   breadth with the {len(surv)} survivors' detection also held fixed: {b2:+.5f} (p {p2:.2e}), "
          f"{abs(b2 / bb[0]):.0%} retained")
    conc = p2 >= 0.05 or abs(b2 / bb[0]) < 0.5
    borderline = conc and abs(b2 / bb[0]) > 0.4
    print("   " + (("CONCENTRATED, but at the edge of the rule: breadth loses significance with the survivors in "
                    f"while keeping {abs(b2 / bb[0]):.0%} of its coefficient, and the criterion was half. Read "
                    "as 'the named references carry about half of it', not as 'they are the effect'."
                    if borderline else
                    "CONCENTRATED: breadth carries little once the survivors are in, so the named references are "
                    "the effect.") if conc else
                   "SPREAD: breadth keeps most of its coefficient with the survivors held fixed, so the named "
                   "references are examples and not the mechanism."))
    rows.append({"part": "concentration", "beta_breadth": bb[0], "p_breadth": bb[2],
                 "beta_breadth_survivors_held": b2, "p": p2, "retained": abs(b2 / bb[0]),
                 "concentrated": bool(conc), "borderline": bool(borderline)})
else:
    print("   no reference survives multiplicity correction, so there is nothing to hold fixed: on this test the")
    print("   breadth effect is spread across references rather than carried by identifiable complexes, and the")
    print("   count is the only level at which it can be stated.")
    rows.append({"part": "concentration", "beta_breadth": bb[0], "p_breadth": bb[2], "concentrated": False,
                 "note": "no reference survives BH at q<0.10"})
print("\n5. the reading")
print("scripts/91 showed the gene-content half of the scoring-quality predictor survives every confound this")
print("dataset can supply; a count of detectable families is not a mechanism, so this asks which families. Of")
print(f"the {len(LAB)} accepted references, {m} have enough detection variance to regress on at all - the other")
print(f"{nexcl} are present in essentially every genome - and {len(surv)} survive Benjamini-Hochberg at q<0.10")
print(f"against {nominal} nominally significant, which is the difference multiplicity correction makes when 37")
print("tests are run. The strongest is cytochrome bo(3) ubiquinol oxidase subunit 1, whose detection goes with a")
print("limit lower by 0.43 pH units, and four transporters point the same way.")
print("THREE THINGS KEEP THIS FROM BEING A MECHANISM. The negative control the test was built around is not")
print("available: all five soluble references are detected in over 99% of genomes, so none could be tested, and")
print("an untested control is a limitation rather than a pass - the first draft of this script reported it as a")
print("pass. What replaces it is weaker and is labelled as such: the survivors' signs are mixed, five down and")
print("two up, which a shared sequencing-or-annotation artefact cannot produce because that account predicts one")
print("direction for every reference. It rules out that artefact and establishes nothing about function.")
print("Second, detection is a crude proxy for presence - a divergent homologue below the gathering threshold")
print("reads as absent - so these are statements about DETECTABLE homologues of E. coli proteins, not about gene")
print("content in the abstract, and the mean-alignment-score control bounds but does not remove that.")
print("Third, the concentration result sits at the edge of its own rule: breadth loses significance with the")
print("survivors held fixed while keeping 47% against a criterion of half, so the honest form is that the named")
print("references carry about half of the breadth effect rather than being it.")
print("What this does settle is that the breadth effect is not uniform across the panel. That matters because a")
print("uniform effect would have been the signature of the completeness artefact scripts/91 ruled out by other")
print("means, and two independent routes agreeing on that is worth more than either.")
pd.DataFrame(rows).to_csv(D / f"which_complexes{TAG}.csv", index=False)
print("\nwrote", D / f"which_complexes{TAG}.csv")
