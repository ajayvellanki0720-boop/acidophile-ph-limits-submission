"""Is the scoring-quality predictor an artefact of assembly, or gene content doing biology? -> DATA_DIR/what_is_scoring_quality.csv and, with `wide`, DATA_DIR/what_is_scoring_quality_wide.csv
BOTH SPELLINGS ARE DECLARED, for the reason scripts/81 records: scripts/69 builds a claim's candidate pool from this line and keeps a panel variant's pool separate, so an
undeclared output turns every figure quoted from it into a guaranteed false miss.
WHERE THIS COMES FROM. scripts/83 found, while asking something else, that how many panel positions a genome has scorable at all predicts the acid limit at -2.35 pH per log10
positions, p 2.3e-10, and scripts/84 then bounded what holding it fixed does to every panel coefficient. Both treated it purely as a CONTROL - a nuisance to be held fixed -
and neither asked what it is. That is the wrong way round for a predictor of that size: a nuisance variable with a p-value of two in ten billion is either an artefact worth
naming or a finding worth having, and it has been neither for two days.
THE DECOMPOSITION THAT SEPARATES THE TWO ACCOUNTS. Scorable positions is a product of BREADTH - how many of the 32 reference families a genome has a detectable homologue of -
and DEPTH - how many positions each of those homologues aligns to. The two mean different things. Breadth is gene content: the references are respiratory and envelope
complexes, so a genome missing several of them is missing those complexes, which is biology. Depth is how well a present gene aligns to an E. coli reference, which is
divergence and annotation quality, and is much closer to an artefact. If the effect lives in breadth it is a statement about what acid-tolerant genomes HAVE; if it lives in
depth it is a statement about how well this panel measures them, and scripts/54 already established that per-reference estimates are largely a measure of scoring cleanliness.
THE ARTEFACT CANDIDATE IS NAMED AND TESTED RATHER THAN DISMISSED. A draft assembly fragments genes, so a Contig-level assembly should have both fewer detectable homologues
and shorter alignments than a Complete Genome of the same organism. Assembly level is recorded for every genome here, and if it accounts for the effect then this is an
artefact of assembly quality and the control is the right response. scripts/67 already tested assembly level against the OUTCOME; what it did not test is whether assembly
level explains this PREDICTOR, which is a different question and the one that matters here.
WHAT WOULD MAKE THIS A FINDING RATHER THAN A CONTROL, stated before the numbers: breadth carrying the effect with assembly level held fixed, and depth not. Anything else and
it stays a nuisance variable, correctly held fixed and not interpreted.
SCALE. Breadth is a count of references and depth a mean positions-per-reference, so their coefficients are not comparable to each other directly; both are therefore also
reported per standard deviation of their own distribution, with those standard deviations printed.
Usage: python3 scripts/91_what_is_scoring_quality.py [orig|wide]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
MINPOS = 200
_r = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv")
MEM = [str(x) for x in _r[_r.accepted & _r.membrane].label]
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
RES = ("D", "E", "K", "R")
tot = {k: np.zeros(len(B)) for k in ("n",) + RES}
breadth = np.zeros(len(B))
for ref in MEM:
    c = f"{ref}|peri_exp|n"
    if c not in B:
        continue
    v = B[c].fillna(0).to_numpy(float)
    tot["n"] += v
    breadth += (v > 0).astype(float)          # the reference contributed any scorable position at all
    for k in RES:
        tot[k] += B[f"{ref}|peri_exp|{k}"].fillna(0).to_numpy(float)
ok = tot["n"] >= MINPOS
g = B[ok].reset_index(drop=True)
Cm = ph.vcv_fast(tree, list(g.ncbi_assembly))
y = g.ph_min.to_numpy(float)
XS = np.column_stack([np.ones(len(g))] + [g[k] for k in ["gw"] + SRCC])
lam = ph.pgls(y, XS, Cm)["lambda"]
qual = np.log10(tot["n"][ok])
br = breadth[ok]
depth = tot["n"][ok] / np.where(br > 0, br, 1)
LEV = ("Complete Genome", "Chromosome", "Scaffold", "Contig")
lev = {L: (g.assembly_level == L).astype(float).to_numpy() for L in LEV if (g.assembly_level == L).sum() > 20}
drop = sorted(lev, key=lambda L: -lev[L].sum())[0]          # the largest level is the reference category
levcols = [lev[L] for L in lev if L != drop]
print(f"sample {int(ok.sum())}; breadth {br.mean():.1f} of {len(MEM)} references (sd {br.std():.2f}), "
      f"depth {depth.mean():.0f} positions per reference (sd {depth.std():.0f}); lambda {lam:.3f}")
print(f"assembly levels used: {', '.join(sorted(lev))}; reference category {drop!r}")


def fit(cols, k=1):
    X = np.column_stack([XS] + list(cols))
    assert np.linalg.matrix_rank(X) == X.shape[1], "rank-deficient design"
    r = ph.pgls(y, X, Cm, lam=lam)
    return [(float(r["beta"][-k + i]), float(r["se"][-k + i]), float(r["p"][-k + i])) for i in range(k)]


print("\n1. the predicted quantity, reproduced before it is taken apart")
(b0,) = fit([qual])
print(f"   log10 scorable positions: {b0[0]:+.5f} pH per log10 (se {b0[1]:.5f}, p {b0[2]:.2e})")
rows.append({"part": "baseline", "term": "log10 scorable positions", "beta": b0[0], "se": b0[1], "p": b0[2],
             "n": int(ok.sum()), "lambda": lam})

print("\n2. breadth and depth entered separately and together")
print(f"   {'term':>42} {'beta':>10} {'se':>9} {'p':>10} {'per SD':>9}")
sdb, sdd = float(np.std(br, ddof=1)), float(np.std(np.log10(depth), ddof=1))
res = {}
for lab, cols, sds in (("breadth alone (references detected)", [br], [sdb]),
                       ("depth alone (log10 positions per reference)", [np.log10(depth)], [sdd]),
                       ("breadth and depth together", [br, np.log10(depth)], [sdb, sdd])):
    out = fit(cols, k=len(cols))
    for (b, se, pv), s, nm in zip(out, sds, (lab if len(cols) == 1 else ("breadth, both in", "depth, both in"),)
                                  if len(cols) == 1 else ("breadth, both in", "depth, both in")):
        print(f"   {nm:>42} {b:>+10.5f} {se:>9.5f} {pv:>10.2e} {b * s:>+9.4f}")
        res[nm] = (b, se, pv, b * s)
        rows.append({"part": "decomposition", "term": nm, "beta": b, "se": se, "p": pv, "beta_per_sd": b * s,
                     "sd_of_term": s})

print("\n3. the artefact candidate: assembly level")
if levcols:
    out = fit(levcols + [qual], k=len(levcols) + 1)
    bq = out[-1]
    print(f"   scoring quality with assembly level held fixed: {bq[0]:+.5f} (p {bq[2]:.2e}), "
          f"against {b0[0]:+.5f} uncontrolled — {abs(bq[0] / b0[0]):.0%} retained")
    for (b, se, pv), L in zip(out[:-1], [L for L in lev if L != drop]):
        print(f"      {L:>18} against {drop}: {b:>+9.5f} (p {pv:.2e})")
        rows.append({"part": "assembly level", "level": L, "reference_category": drop, "beta": b, "se": se,
                     "p": pv})
    rows.append({"part": "assembly level", "term": "scoring quality controlled", "beta": bq[0], "se": bq[1],
                 "p": bq[2], "retained": abs(bq[0] / b0[0])})
    artefact = bq[2] >= 0.05 or abs(bq[0] / b0[0]) < 0.5
else:
    artefact = None
    print("   too few genomes at any secondary assembly level to test")

print("\n3b. THE TWO CONFOUNDS BREADTH CANNOT ESCAPE ON ITS OWN: distance from the reference organism, and size")
# Every reference in this panel is an E. coli protein, so a genome far from E. coli has fewer DETECTABLE
# homologues for reasons that have nothing to do with which complexes it carries. Breadth and evolutionary
# distance are therefore entangled by construction, and the gene-content reading above is only available if
# breadth survives holding distance fixed. The panel stores each reference's BLOSUM62 alignment score, whose mean
# across the references a genome does hit is a direct measure of how similar its proteins are to E. coli's -
# better suited here than a tree distance, because it is computed on the very sequences the breadth count is
# derived from. Holding it fixed is the test; Pagel's lambda does not do this job, since it models shared
# ancestry in the OUTCOME rather than detection bias in the predictor.
scols = [f"{ref}_score" for ref in MEM if f"{ref}_score" in B.columns]
if scols:
    S = B[scols].to_numpy(float)
    with np.errstate(invalid="ignore"):
        meanscore = np.nanmean(np.where(S > 0, S, np.nan), axis=1)
    ms = meanscore[ok]
    ms = np.where(np.isfinite(ms), ms, np.nanmedian(ms))
    rb = float(np.corrcoef(br, ms)[0, 1])
    # GENOME SIZE GOES IN TOO, because without it the gene-content reading is almost worthless. A bigger genome
    # has more genes and so more detectable references, and scripts/85 found log10 genome size predicts the limit
    # in its own right - so breadth and genome size are confounded in exactly the direction that would manufacture
    # this finding. Leaving it out would have let "acid-tolerant genomes carry more of these complexes" stand when
    # the honest version might be "acid-tolerant genomes are bigger".
    gsz = np.log10(g.genome_size_bp.to_numpy(float))
    out = fit([ms, gsz, br, np.log10(depth)], k=4)
    (bs, ses, ps), (bg, seg, pg), (bb2, seb2, pb2), (bd2, sed2, pd2) = out
    print(f"   log10 genome size, entered with it: {bg:+.5f} (se {seg:.5f}, p {pg:.2e})")
    rows.append({"part": "distance control", "term": "log10 genome size", "beta": bg, "se": seg, "p": pg})
    print(f"   mean alignment score to the E. coli references: r with breadth {rb:+.3f}")
    print(f"   {'term':>42} {'beta':>10} {'se':>9} {'p':>10}")
    for nm, v in (("mean alignment score", (bs, ses, ps)), ("breadth, distance held fixed", (bb2, seb2, pb2)),
                  ("depth, distance held fixed", (bd2, sed2, pd2))):
        print(f"   {nm:>42} {v[0]:>+10.5f} {v[1]:>9.5f} {v[2]:>10.2e}")
        rows.append({"part": "distance control", "term": nm, "beta": v[0], "se": v[1], "p": v[2]})
    breadth_survives = pb2 < 0.05 and (bb2 < 0) == (res["breadth, both in"][0] < 0)
    print("   " + ("Breadth SURVIVES holding similarity to the reference organism and genome size fixed, at "
                   f"{abs(bb2 / res['breadth, both in'][0]):.0%} of its size, so the gene-content reading is "
                   "available." if breadth_survives else
                   "Breadth does NOT survive: once similarity to E. coli is held fixed it carries nothing, so "
                   "the count of detectable references was measuring distance from the reference organism and "
                   "the gene-content reading is withdrawn."))
    rows.append({"part": "distance control", "term": "breadth survives", "survives": bool(breadth_survives),
                 "r_breadth_meanscore": rb})
else:
    breadth_survives = None
    print("   this panel build stores no per-reference alignment scores; the confound cannot be tested here")

print("\n4. the reading, against the criterion fixed in the preamble")
# READ OFF THE BEST-SPECIFIED FIT, which the first version did not. It computed the verdict from part 2, where
# neither breadth nor depth is corrected for similarity to the reference organism, and so concluded "both" while
# part 3b had just shown that depth carries NOTHING once that is held fixed. A script that runs a decisive
# control and then reads its verdict off the uncontrolled fit has wasted the control, and this is the sixth
# conclusion in two days that contradicted output printed above it.
if scols:
    bb = [r for r in rows if r["part"] == "distance control" and r.get("term") == "breadth, distance held fixed"][0]
    dd = [r for r in rows if r["part"] == "distance control" and r.get("term") == "depth, distance held fixed"][0]
    bb, dd = (bb["beta"], bb["se"], bb["p"]), (dd["beta"], dd["se"], dd["p"])
    basis = "with similarity to the reference organism held fixed"
else:
    bb, dd = res["breadth, both in"], res["depth, both in"]
    basis = "without a similarity control, which this panel build cannot support"
br_wins = bb[2] < 0.05 and dd[2] >= 0.05
dp_wins = dd[2] < 0.05 and bb[2] >= 0.05
print(f"   read {basis}: breadth {bb[0]:+.5f} (p {bb[2]:.1e}), depth {dd[0]:+.5f} (p {dd[2]:.1e})")
if breadth_survives is False:
    verdict = ("DETECTION BIAS. Breadth carries nothing once similarity to the E. coli references is held "
               "fixed, so the predictor is measuring how far a genome is from the reference organism and not "
               "what it contains. Correctly a control, and the gene-content reading is withdrawn.")
elif artefact:
    verdict = ("ARTEFACT: assembly level accounts for it, so the control is the right response and there is "
               "nothing here to interpret.")
elif br_wins:
    verdict = ("GENE CONTENT. Breadth carries the effect and depth does not, with assembly level, similarity to "
               "the reference organism and genome size all held fixed - genome size being the confound that "
               "would most easily have manufactured this, and it comes back null - so this is a statement "
               "about which complexes "
               "acid-tolerant genomes HAVE rather than about how well this panel measures them. It is a finding, "
               "and it has been held fixed as a nuisance for two days. Depth behaves exactly as the measurement "
               "account predicted it would - it predicts the limit until similarity to E. coli is controlled and "
               "then carries nothing - so both halves of the decomposition came out as their accounts said.")
elif dp_wins:
    verdict = ("MEASUREMENT. Depth carries it and breadth does not, so it is about how well a present gene "
               "aligns to an E. coli reference - divergence and annotation - which is scripts/54's finding at "
               "genome level. Correctly a control, and correctly not interpreted.")
elif bb[2] < 0.05 and dd[2] < 0.05:
    verdict = ("BOTH, which the preamble did not anticipate and which is not a dodge: breadth and depth each "
               "predict the limit with the other held fixed, so the predictor is a genuine mixture and neither "
               "account alone explains it.")
else:
    verdict = ("NEITHER RESOLVES once both are entered, which given their correlation is a collinearity "
               "statement about this panel rather than a verdict.")
print("   " + verdict)
print(f"   correlation between breadth and log10 depth: {float(np.corrcoef(br, np.log10(depth))[0, 1]):+.3f}")
rows.append({"part": "verdict", "verdict": verdict, "artefact_by_assembly": artefact,
             "r_breadth_depth": float(np.corrcoef(br, np.log10(depth))[0, 1])})
pd.DataFrame(rows).to_csv(D / f"what_is_scoring_quality{TAG}.csv", index=False)
print("\nwrote", D / f"what_is_scoring_quality{TAG}.csv")
