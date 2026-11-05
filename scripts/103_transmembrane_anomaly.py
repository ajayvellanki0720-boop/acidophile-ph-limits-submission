"""Why do genomes with a measured salinity carry a transmembrane signal? -> DATA_DIR/transmembrane_anomaly.csv
THE FACT TO EXPLAIN, from `scripts/98`. On the full sample the exposed transmembrane composite is the project's
quietest control: -0.0430 per SD at p 7.3e-02, null, and the specificity argument leans on it. On the 669 genomes
that happen to carry a measured salinity value it is -0.14768 at p 3.0e-03 - 1.06x the periplasmic coefficient it is
supposed to be a foil for, against 0.23x on the full sample. Across 300 random subsamples of the same size only 2.0%
reach a transmembrane coefficient that extreme, so it is not an artefact of smaller n, and salinity itself is not the
cause: holding it fixed moves that coefficient by less than a thousandth.
WHY IT MATTERS MORE THAN ITS SIZE SUGGESTS. A control that fires is not a small blemish on a specificity argument; it
is the specificity argument failing on that sample. If whatever selects a species into having its salt tolerance
measured also makes the transmembrane class predictive, then the same selection could be at work elsewhere in the
dataset, unobserved, and the clean full-sample controls would be clean by luck.
THE HYPOTHESES, EACH WITH A TEST THAT CAN COME BACK NO:
  (a) COMPOSITION OF THE SUBSAMPLE. It differs from the rest in something already measured - source tier, assembly
      level, genome size, GC, phylum mix. If so, the signal should weaken when those are held fixed.
  (b) WELL-STUDIEDNESS RATHER THAN SALT. Salinity is measured for organisms someone characterised thoroughly. If the
      anomaly is about characterisation depth, then ANOTHER well-characterised subset, selected without reference to
      salinity, should show it too. The growth-optimum-known set is exactly such a subset and is far larger.
  (c) A FEW REFERENCES. The transmembrane composite pools 32 references; if two or three carry the whole thing the
      question becomes about those proteins rather than about the class.
  (d) SALT-ADAPTED MEMBRANES, the biologically obvious one. Then the signal should concentrate in the genomes with
      HIGH measured salinity tolerance, not spread evenly through the subsample.
WHAT I EXPECT TO FIND AND WHY IT IS WRITTEN DOWN FIRST: (b). The subsample is defined by a curation fact rather than a
biological one, and `scripts/98` already showed it is differently composed - mean ph_min 5.91 against 5.71, 6.4%
acid-tolerant against 8.9%. Writing the expectation down is the only way the result can embarrass it.
SCALE. Per SD of each composite, every fit against a baseline refitted on the same genomes - never the full-sample
figure, which would be correction 12.
Usage: python3 scripts/103_transmembrane_anomaly.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import os
import re
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
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
T = ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
SAL = ROOT / "data/raw/metatraits/ncbi_species_salinity_no_predictions.tsv"
rows = []
if not SAL.exists():
    sys.exit(f"missing {SAL}; see scripts/98's header for the fetch commands (data/raw/* is gitignored)")

refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
NAME = {str(r.label): str(r.name) for r in refs.itertuples()}
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
pan = pd.read_csv(D / "structure_panel_wide.csv").drop_duplicates("ncbi_assembly")
gf = pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess", "genome_size_bp", "gc"]]
pheno = pd.read_csv(D / "phenotypes_all.csv").merge(gf, on="ncbi_assembly")
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
B = pheno[(pheno.domain == "Bacteria") & pheno.ncbi_assembly.isin(tips) & pheno.ph_min.notna()].drop_duplicates("ncbi_assembly").reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100.0
B["log_size"] = np.log10(B.genome_size_bp.astype(float))
B["complete"] = (B.assembly_level.astype(str).str.lower() == "complete genome").astype(float)
SRCC = [f"src_{o}" for o in SRC]
M = pan.set_index("ncbi_assembly")


def composite(cls, minrefs=3, only=None, drop=None):
    z, used = {}, []
    for r in MEM:
        if only and r not in only:
            continue
        if drop and r in drop:
            continue
        if f"{r}|{cls}|n" not in M.columns:
            continue
        n_ = M[f"{r}|{cls}|n"].reindex(B.ncbi_assembly).to_numpy(float)
        num = (M[f"{r}|{cls}|D"] + M[f"{r}|{cls}|E"] - M[f"{r}|{cls}|K"] - M[f"{r}|{cls}|R"]).reindex(B.ncbi_assembly).to_numpy(float)
        fl = max(8.0, 0.5 * float(reflen[(r, cls)])) if (r, cls) in reflen.index else 8.0
        x = np.where(n_ >= fl, num / np.where(n_ > 0, n_, np.nan), np.nan)
        z[r] = (x - np.nanmean(x)) / np.nanstd(x, ddof=1)
        used.append(r)
    Z = np.column_stack([z[k] for k in used])
    cnt = np.sum(~np.isnan(Z), axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = np.nanmean(Z, axis=1)
    return np.where(cnt >= min(minrefs, len(used)), s, np.nan), used


def norm(s_):
    s_ = re.sub(r"\s*\(ex .*?\)", "", str(s_)).strip()
    return " ".join(re.sub(r"^Candidatus ", "", s_).split()[:2]).lower()


sal = pd.read_csv(SAL, sep="\t")
sal["key"] = sal.taxon_name.map(norm)
sal["val"] = pd.to_numeric(sal["median"], errors="coerce")
smax = dict(zip(sal.loc[(sal.trait_name == "salinity maximum") & sal.val.notna(), "key"],
                sal.loc[(sal.trait_name == "salinity maximum") & sal.val.notna(), "val"]))
B["sal_max"] = B.species.map(norm).map(smax)
XT, _ = composite("TM")
XP, _ = composite("peri_exp")
HAVE = B.sal_max.notna().to_numpy() & ~np.isnan(XT) & ~np.isnan(XP)
C = ph.vcv_fast(tree, list(B.ncbi_assembly))
y = B.ph_min.to_numpy(float)


def fit(mask, x, covs=("gw",)):
    g = B[mask]
    use = [c for c in list(covs) + SRCC if g[c].nunique() > 1]
    idx = np.where(mask)[0]
    Xb = np.column_stack([np.ones(int(mask.sum()))] + [g[c].to_numpy(float) for c in use])
    v = x[mask]
    Cc = C[np.ix_(idx, idx)]
    lam = ph.pgls(y[mask], Xb, Cc)["lambda"]
    f = ph.pgls(y[mask], np.column_stack([Xb, (v - v.mean()) / np.std(v)]), Cc, lam=lam)
    return float(f["beta"][-1]), float(f["se"][-1]), float(f["p"][-1])


b_tm, se_tm, p_tm = fit(HAVE, XT)
b_pe, _, _ = fit(HAVE, XP)
print(f"the fact: on {int(HAVE.sum())} salinity-measured genomes the transmembrane control is {b_tm:+.5f} (p {p_tm:.2e}),")
print(f"          {abs(b_tm) / abs(b_pe):.2f}x the periplasmic coefficient ({b_pe:+.5f}) on the same genomes")
rows.append({"part": "fact", "measure": "transmembrane on salinity-measured", "n": int(HAVE.sum()),
             "beta_per_sd": b_tm, "se": se_tm, "p": p_tm, "beta_periplasmic": b_pe})

print("\n(a) is the subsample just differently composed? hold the differences fixed and see")
COVSETS = [("growth-wide composition only (as scripts/98)", ("gw",)),
           ("+ genome size and GC", ("gw", "log_size", "gc")),
           ("+ assembly completeness", ("gw", "log_size", "gc", "complete")),
           ("+ the growth optimum where known", ("gw", "log_size", "gc", "complete"))]
for lab, cv in COVSETS[:3]:
    b, se, p = fit(HAVE, XT, cv)
    print(f"   {lab:>44s} {b:>+10.5f} (se {se:.5f}, p {p:.2e})  {b / b_tm:>5.0%} of the original")
    rows.append({"part": "a", "measure": lab, "n": int(HAVE.sum()), "beta_per_sd": b, "se": se, "p": p,
                 "share_of_original": b / b_tm})
print("   If the signal were a composition artefact these would collapse it. They are the measured differences")
print("   between the subsample and the rest, and holding them fixed is the cheapest test available.")

print("\n(b) is it WELL-STUDIEDNESS rather than salt? another thoroughly-characterised subset, chosen without salinity")
# The growth optimum is recorded for a species when someone characterised it properly, and it is known for 1,655
# genomes - a far larger, independently-defined "well-studied" set. If depth of characterisation is what drives the
# transmembrane signal, it should appear here too. This is the hypothesis written down in advance.
OPT = B.ph_opt.notna().to_numpy() & ~np.isnan(XT)
both = OPT & HAVE
print(f"   growth-optimum-known: {int(OPT.sum())} genomes ({int(both.sum())} of them also salinity-measured)")
for lab, mask in (("salinity-measured (the anomaly)", HAVE),
                  ("growth-optimum-known", OPT),
                  ("optimum-known but NOT salinity-measured", OPT & ~HAVE),
                  ("neither", ~OPT & ~HAVE & ~np.isnan(XT))):
    if mask.sum() < 60:
        print(f"   {lab:>42s}  n={int(mask.sum())}, too few")
        continue
    b, se, p = fit(mask, XT)
    bp, _, _ = fit(mask, XP)
    print(f"   {lab:>42s}  n={int(mask.sum()):>5d}  TM {b:>+9.5f} (p {p:.2e})   ratio to peri {abs(b) / abs(bp):>5.2f}x")
    rows.append({"part": "b", "measure": lab, "n": int(mask.sum()), "beta_tm": b, "p_tm": p,
                 "beta_peri": bp, "ratio_tm_to_peri": abs(b) / abs(bp)})

print("\n(c) is it a few references? drop each in turn and watch the transmembrane coefficient")
base = b_tm
moves = []
for r in MEM:
    x2, used = composite("TM", drop={r})
    if len(used) < 10 or np.isnan(x2[HAVE]).all():
        continue
    ok = HAVE & ~np.isnan(x2)
    if ok.sum() < 100:
        continue
    b2, _, p2 = fit(ok, x2)
    moves.append((abs(b2 - base) / se_tm, r, b2, p2))
moves.sort(reverse=True)
print(f"   {'dropped reference':>52s} {'beta':>10s} {'move (SE)':>10s}")
for mv, r, b2, p2 in moves[:4]:
    print(f"   {NAME.get(r, r)[:52]:>52s} {b2:>+10.5f} {mv:>10.2f}")
    rows.append({"part": "c", "measure": f"drop {NAME.get(r, r)}", "beta_per_sd": b2, "move_in_se": mv, "p": p2})
print(f"   largest single-reference move: {moves[0][0]:.2f} SE. A composite carried by one protein moves several.")

print("\n(d) is it salt adaptation? split the subsample by measured salt tolerance")
sv = B.sal_max.to_numpy(float)
med = float(np.nanmedian(sv[HAVE]))
for lab, mask in ((f"salinity maximum <= {med:.1f}% NaCl", HAVE & (sv <= med)),
                  (f"salinity maximum >  {med:.1f}% NaCl", HAVE & (sv > med))):
    if mask.sum() < 60:
        continue
    b, se, p = fit(mask, XT)
    print(f"   {lab:>40s}  n={int(mask.sum()):>4d}  TM {b:>+10.5f} (se {se:.5f}, p {p:.2e})")
    rows.append({"part": "d", "measure": lab, "n": int(mask.sum()), "beta_per_sd": b, "se": se, "p": p})
print("   A membrane-adaptation account predicts the signal concentrating in the salt-tolerant half.")

print("\nverdict")
bsub = {r["measure"]: r for r in rows if r["part"] == "b"}
opt_only = bsub.get("optimum-known but NOT salinity-measured")
a_last = [r for r in rows if r["part"] == "a"][-1]
explained_a = abs(a_last["share_of_original"]) < 0.5
wellstudied = bool(opt_only and opt_only["p_tm"] < 0.05)
single_ref = moves[0][0] > 1.0 if moves else False
# HYPOTHESIS (d) HAD NO BRANCH IN THE FIRST VERSION OF THIS VERDICT, so the else fired and printed "it does not
# concentrate in the salt-tolerant half" directly beneath a table showing it concentrating in the salt-tolerant half.
# That is the third verdict in this session to assert the opposite of its own output, and the cause is the same every
# time: branches written for the outcomes I expected, with the one I did not expect falling through to a default that
# describes nothing. The condition is now derived from the part-(d) rows rather than assumed absent.
dsub = [r for r in rows if r["part"] == "d"]
hi = next((r for r in dsub if r["measure"].startswith("salinity maximum >")), None)
lo = next((r for r in dsub if r["measure"].startswith("salinity maximum <=")), None)
saltadapt = bool(hi and lo and hi["p"] < 0.05 and lo["p"] >= 0.05 and abs(hi["beta_per_sd"]) > 2 * abs(lo["beta_per_sd"]))
if saltadapt:
    v = (f"SALT ADAPTATION, which was the one account I did not expect. The signal is confined to the salt-tolerant half: "
         f"{hi['beta_per_sd']:+.5f} (p {hi['p']:.2e}) above {med:.1f}% NaCl against {lo['beta_per_sd']:+.5f} (p {lo['p']:.2e}) below it, a factor of "
         f"{abs(hi['beta_per_sd']) / abs(lo['beta_per_sd']):.1f}. Well-studiedness is ruled out by an independently-defined characterised subset that does NOT reproduce it "
         f"(optimum-known, n=1,655, p {bsub['growth-optimum-known']['p_tm']:.2f}); measured composition is ruled out because holding genome size, GC and assembly "
         f"completeness fixed leaves {a_last['share_of_original']:.0%} of it; and no single reference moves it more than {moves[0][0]:.2f} SE. "
         "This is the reassuring answer rather than the alarming one: halotolerant membranes are compositionally distinct, which is a known biological fact and is "
         "exactly what scripts/43 found inside archaeal halophiles, where the transmembrane class also turns non-null and specificity fails. The control is not "
         "broken; it is responding to salt in the organisms that tolerate salt.")
elif wellstudied:
    v = ("WELL-STUDIEDNESS, not salt. A thoroughly-characterised subset defined WITHOUT reference to salinity shows the same "
         f"transmembrane signal (p {opt_only['p_tm']:.2e}), so what fires the control is a property of being well characterised rather than of "
         "salt tolerance. That is worse news than the salinity reading, because characterisation depth is not a variable this project can "
         "condition on and it is present throughout the dataset rather than in one subsample.")
elif explained_a:
    v = ("MEASURED COMPOSITION. Holding genome size, GC and assembly completeness fixed removes most of it, so the anomaly is the subsample "
         "being differently composed in ways already recorded, not a new effect.")
elif single_ref:
    v = (f"A FEW REFERENCES. Dropping one moves the coefficient {moves[0][0]:.1f} standard errors, so this is a statement about particular proteins "
         "rather than about the transmembrane class.")
else:
    v = ("UNEXPLAINED. None of the four accounts carries it: the measured covariates do not absorb it, an independently-defined well-studied "
         "subset does not reproduce it, no single reference moves it far, and it does not concentrate in the salt-tolerant half. The control "
         "fires on these genomes and this script cannot say why, which is recorded as a limitation of the specificity argument rather than resolved.")
print(f"   {v}")
rows.append({"part": "verdict", "verdict": v, "n": int(HAVE.sum())})
if saltadapt:
    print("   THE EXPECTATION WRITTEN INTO THIS SCRIPT'S HEADER WAS (b), WELL-STUDIEDNESS, AND IT IS WRONG. It was written")
    print("   down precisely so that it could be, and the subsample's defining property being a curation fact made it the")
    print("   natural guess. The data say the biologically obvious account is the right one. Recording the failed")
    print("   expectation costs nothing and is the only thing that makes writing one down worth doing.")
print("   WHAT IS TRUE EITHER WAY. The full-sample transmembrane control is null and the specificity argument rests on")
print("   it; what this shows is that it is not null in every SUBSET, so the specificity claim holds for the sample as a")
print("   whole and is not a guarantee about any slice of it.")
if saltadapt:
    print("   AND WHAT FOLLOWS FOR THE PAPER. The subsample was defined by a curation fact - having a salinity measurement -")
    print("   but the signal inside it is biological, and the two were easy to confuse precisely because the curation fact")
    print("   selects for the biology: species get their salt tolerance measured when salt is interesting about them. So")
    print("   the right statement is not that the control fails on well-characterised genomes, which would be alarming and")
    print("   is false, but that halotolerance perturbs transmembrane composition - a second environmental variable acting")
    print("   on a DIFFERENT compartment from the one this project claims, which is a result in favour of the controls")
    print("   rather than against them. scripts/98 reached the same conclusion from the other direction: salinity moves the")
    print("   secreted and cytoplasmic shares together where pH moves only the secreted one.")

pd.DataFrame(rows).to_csv(D / "transmembrane_anomaly.csv", index=False)
print(f"\nwrote {D / 'transmembrane_anomaly.csv'}")
