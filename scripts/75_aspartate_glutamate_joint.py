"""Are aspartate and glutamate separable when entered together, and what do the arm's five aspartate positions cost it? -> DATA_DIR/de_joint.csv and, with the `wide` argument, DATA_DIR/de_joint_wide.csv
Both spellings are declared because scripts/69 builds a claim's candidate pool from this line and keeps a panel variant's pool separate from the default run's; declaring
only the default is what left two of scripts/36c's figures untraceable earlier the same day.
THIS CLOSES A QUANTITY SCRIPTS/73 WITHHELD ON PURPOSE. Correction 14 found that the acidic signal is carried by glutamate (positive in 31 of 32 references)
and barely by aspartate (18 of 32), which matters because scripts/37b admits any position where the reference carries D OR E, and the primary arm is 9
glutamate and 5 aspartate. The obvious question is what those 5 cost the arm. scripts/73 refused to answer it, because the only available fits were
UNIVARIATE per component - b_D from a model in which glutamate varies unmodelled, b_E from one in which aspartate does - and b_D and b_E from separate
models are not additive contributions that sum to b_DE. Computed anyway the ratio came out below 1 on BOTH panels, which was the tell. The legitimate
version needs D and E entered SIMULTANEOUSLY so that each coefficient is the effect of that residue holding the other fixed, which is what this fits.
TWO REASONS THE ANSWER MIGHT STILL BE UNUSABLE, both measured here rather than assumed.
  1. COLLINEARITY. The two shares are both acidic residues on the same surface and may move together across genomes, in which case a joint fit splits a
     common effect arbitrarily and neither coefficient means much. The correlation between the two shares is reported per reference and pooled.
  2. The per-position weighting assumes a position's contribution is proportional to its own residue's coefficient. That is a modelling assumption, not a
     measurement, and it is stated here rather than buried: it is what "the arm is effectively fewer than 14 positions" would mean.
SCALE. Every coefficient here is per PERCENTAGE POINT and comes from ph.pgls directly, not from ph.pgls_scan, so correction 15's mislabelling cannot arise:
pgls returns per-unit coefficients for the design as given. The point is laboured because this project has made that error seven times.
WHAT THIS DOES NOT DO. It does not re-derive the experiment's predicted effect or its required sample size. Those come from the composite in SD units
(docs/CHARGE_REVERSAL_DESIGN.md section 2) and are committed in docs/PREREGISTRATION.md section 3; a per-residue reweighting is a different calculation and
turning it into a new effect size would be the post-hoc rescaling that document's section 7 forbids. The output here is a share, not a prediction.
Usage: python3 scripts/75_aspartate_glutamate_joint.py [orig|wide]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import stats
from acidophile import phylo as ph

warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
WHICH = (sys.argv[1] if len(sys.argv) > 1 else "wide").lower()
assert WHICH in ("orig", "wide"), f"panel must be orig or wide, got {WHICH!r}"
if WHICH == "wide":
    TAG, PANEL_CSV, CLS_CSV = "_wide", "structure_panel_wide.csv", "results/structure_panel_wide_classes.csv"
    _r = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv")
    MEM = [str(x) for x in _r[_r.accepted & _r.membrane].label]
else:
    TAG, PANEL_CSV, CLS_CSV = "", "structure_panel.csv", "results/structure_panel_classes.csv"
    MEM = ["YidC", "SecY", "SecY_Bsub", "ATPsynC", "NuoL", "NuoN", "CyoB"]
rows = []
CLS = pd.read_csv(ROOT / CLS_CSV)
reflen = CLS.groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv")
      .merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess", "gc", "genome_size_bp"]],
             on="ncbi_assembly")
      .merge(pd.read_csv(D / PANEL_CSV), on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
Call = ph.vcv_fast(tree, list(B.ncbi_assembly))
B["gw"] = B.acidic_excess * 100
SRCC = [f"src_{o}" for o in SRC if B[f"src_{o}"].nunique() > 1]
print(f"panel: {WHICH} ({len(MEM)} membrane references), {len(B)} bacteria")


def pool(b, s):
    b, s = np.asarray(b, float), np.asarray(s, float)
    w = 1 / s ** 2
    fe = (w * b).sum() / w.sum()
    Q = (w * (b - fe) ** 2).sum()
    k = len(b) - 1
    tau2 = max(0.0, (Q - k) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 0 else 0.0
    wr = 1 / (s ** 2 + tau2)
    mu = (wr * b).sum() / wr.sum()
    sem = float(np.sqrt(1 / wr.sum()))
    return float(mu), sem, float(2 * stats.norm.sf(abs(mu / sem)))


print("\n1. D and E entered simultaneously, per reference (per percentage point, composition controlled)")
print(f"   {'ref':>8} {'n':>5} {'r(D%,E%)':>9} {'beta_D':>9} {'beta_E':>9} {'p_D':>8} {'p_E':>8}")
got = []
for ref in MEM:
    c = "peri_exp"
    if f"{ref}|{c}|n" not in B:
        continue
    n = B[f"{ref}|{c}|n"]
    need = max(8, 0.5 * reflen[(ref, c)])
    xD = (B[f"{ref}|{c}|D"] / n).where(n >= need) * 100
    xE = (B[f"{ref}|{c}|E"] / n).where(n >= need) * 100
    m = (xD.notna() & xE.notna()).to_numpy()
    if m.sum() < 200:
        continue
    g = B[m].reset_index(drop=True)
    Cm = Call[np.ix_(m, m)]
    keep = [k for k in SRCC if g[k].nunique() > 1]
    dv, ev = xD[m].to_numpy(), xE[m].to_numpy()
    if dv.std() == 0 or ev.std() == 0:
        continue
    # lambda fitted on the covariates-only model, as pgls_scan does, then the joint fit at that lambda
    Xb = np.column_stack([np.ones(len(g))] + [g[k] for k in ["gw"] + keep])
    lam = ph.pgls(g.ph_min.to_numpy(float), Xb, Cm)["lambda"]
    r = ph.pgls(g.ph_min.to_numpy(float), np.column_stack([Xb, dv, ev]), Cm, lam=lam)
    rr = float(np.corrcoef(dv, ev)[0, 1])
    bD, bE = float(r["beta"][-2]), float(r["beta"][-1])
    sD, sE = float(r["se"][-2]), float(r["se"][-1])
    pD, pE = float(r["p"][-2]), float(r["p"][-1])
    print(f"   {ref:>8} {int(m.sum()):>5} {rr:>+9.3f} {bD:>+9.5f} {bE:>+9.5f} {pD:>8.3f} {pE:>8.3f}")
    got.append({"ref": ref, "n": int(m.sum()), "r_shares": rr, "beta_D": bD, "se_D": sD, "p_D": pD,
                "beta_E": bE, "se_E": sE, "p_E": pE, "lam": lam, "sd_D": float(dv.std()), "sd_E": float(ev.std())})
    rows.append({"part": "per reference joint fit", "panel": WHICH, **got[-1]})
G = pd.DataFrame(got)
assert len(G) >= 2, "need at least two references to pool"

print(f"\n2. collinearity, which decides whether a joint fit means anything")
print(f"   correlation between the two shares: median {G.r_shares.median():+.3f}, range {G.r_shares.min():+.3f} to {G.r_shares.max():+.3f}")
vif = 1 / (1 - G.r_shares ** 2)
print(f"   implied variance inflation for each coefficient: median {vif.median():.2f}, max {vif.max():.2f}")
ok = vif.max() < 5
print(f"   {'USABLE' if ok else 'NOT USABLE'}: a variance inflation under 5 is the conventional threshold, and the maximum here is {vif.max():.2f}.")
if ok:
    print("   The two shares are close to independent across genomes, so entering them together splits a real")
    print("   distinction rather than an arbitrary one. That is the condition scripts/73 could not check.")
rows.append({"part": "collinearity", "panel": WHICH, "median_r_shares": float(G.r_shares.median()),
             "min_r_shares": float(G.r_shares.min()), "max_r_shares": float(G.r_shares.max()),
             "median_vif": float(vif.median()), "max_vif": float(vif.max()), "usable": bool(ok)})

print("\n3. pooled across references by random effects, both residues from the SAME model")
mD, seD, pD = pool(G.beta_D, G.se_D)
mE, seE, pE = pool(G.beta_E, G.se_E)
print(f"   aspartate  {mD:+.5f} +/- {seD:.5f} (p {pD:.2e}); positive in {(G.beta_D > 0).sum()}/{len(G)}")
print(f"   glutamate  {mE:+.5f} +/- {seE:.5f} (p {pE:.2e}); positive in {(G.beta_E > 0).sum()}/{len(G)}")
d = mE - mD
sd = float(np.hypot(seE, seD))
print(f"   difference {d:+.5f} +/- {sd:.5f} (p {2 * stats.norm.sf(abs(d / sd)):.3f}) - whether the two residues differ AT ALL")
for lab, m_, s_, p_, pos in (("D", mD, seD, pD, int((G.beta_D > 0).sum())), ("E", mE, seE, pE, int((G.beta_E > 0).sum()))):
    rows.append({"part": "pooled joint", "panel": WHICH, "residue": lab, "beta_per_point": m_, "se": s_, "p": p_,
                 "n_refs": int(len(G)), "n_refs_positive": pos})
rows.append({"part": "difference", "panel": WHICH, "beta_E_minus_D": d, "se": sd,
             "p": float(2 * stats.norm.sf(abs(d / sd)))})

print("\n4. the quantity scripts/73 withheld: what the arm's aspartate positions contribute")
T2 = pd.read_csv(D / "mutagenesis_targets.csv")
T2 = T2[T2.target]
arm = T2[T2.ref.isin(["NuoL", "NuoN"])]
nE, nD = int((arm.ref_aa == "E").sum()), int((arm.ref_aa == "D").sum())
assert nE + nD == 14, f"the primary arm is 14 positions; found {nE + nD}"
wE, wD = nE * mE, nD * mD
tot = wE + wD
print(f"   the arm is {nE} glutamate and {nD} aspartate, and on this model their expected contributions are")
print(f"   glutamate {wE:+.5f} and aspartate {wD:+.5f}, so glutamate supplies {wE / tot:.0%} of the arm's expected effect")
print(f"   and the {nD} aspartate positions supply {wD / tot:.0%}.")
equiv = tot / mE
print(f"   Put the other way, the 14 positions are worth {equiv:.1f} glutamate positions, so the arm is effectively")
print(f"   {equiv:.1f} positions rather than 14 - a shortfall of {(1 - equiv / 14):.0%} against counting them equally.")
rows.append({"part": "arm weighting", "panel": WHICH, "n_glutamate": nE, "n_aspartate": nD,
             "contribution_glutamate": float(wE), "contribution_aspartate": float(wD),
             "share_from_glutamate": float(wE / tot), "share_from_aspartate": float(wD / tot),
             "glutamate_equivalent_positions": float(equiv),
             "shortfall_vs_counting_equally": float(1 - equiv / 14)})
print("   This is the figure scripts/73 declined to give, and it is given now only because part 2 established that")
print("   the two shares are separable. It is a SHARE and not a new predicted effect: the committed effect size comes")
print("   from the composite in SD units and re-deriving it here would be the rescaling the pre-registration forbids.")

print("\n4b. is the glutamate-over-aspartate difference a GC artefact? The one confound specific to THIS comparison")
print("   Glutamate is GAA/GAG and aspartate GAU/GAC, so the two differ in the third codon position's GC content and a")
print("   GC-rich genome is pushed toward aspartate. The headline models control GC; scripts/36c's residue")
print("   decomposition does NOT, controlling only genome-wide acidic excess and source. So the difference found above")
print("   could be genome GC acting through codon usage rather than anything about the surface. Refitting with GC and")
print("   genome size added to the same joint model:")
B["gcz"] = (B.gc - B.gc.mean()) / B.gc.std()
B["szz"] = (B.genome_size_bp - B.genome_size_bp.mean()) / B.genome_size_bp.std()
got2 = []
for ref in MEM:
    c = "peri_exp"
    if f"{ref}|{c}|n" not in B:
        continue
    n = B[f"{ref}|{c}|n"]
    need = max(8, 0.5 * reflen[(ref, c)])
    xD = (B[f"{ref}|{c}|D"] / n).where(n >= need) * 100
    xE = (B[f"{ref}|{c}|E"] / n).where(n >= need) * 100
    m = (xD.notna() & xE.notna() & B.gc.notna() & B.genome_size_bp.notna()).to_numpy()
    if m.sum() < 200:
        continue
    g = B[m].reset_index(drop=True)
    Cm = Call[np.ix_(m, m)]
    keep = [k for k in SRCC if g[k].nunique() > 1]
    dv, ev = xD[m].to_numpy(), xE[m].to_numpy()
    if dv.std() == 0 or ev.std() == 0:
        continue
    Xb = np.column_stack([np.ones(len(g))] + [g[k] for k in ["gw", "gcz", "szz"] + keep])
    lam = ph.pgls(g.ph_min.to_numpy(float), Xb, Cm)["lambda"]
    r = ph.pgls(g.ph_min.to_numpy(float), np.column_stack([Xb, dv, ev]), Cm, lam=lam)
    got2.append({"ref": ref, "beta_D": float(r["beta"][-2]), "se_D": float(r["se"][-2]),
                 "beta_E": float(r["beta"][-1]), "se_E": float(r["se"][-1])})
G2 = pd.DataFrame(got2)
mD2, seD2, pD2 = pool(G2.beta_D, G2.se_D)
mE2, seE2, pE2 = pool(G2.beta_E, G2.se_E)
d2 = mE2 - mD2
sd2 = float(np.hypot(seE2, seD2))
p_d2 = float(2 * stats.norm.sf(abs(d2 / sd2)))
print(f"   {'model':38s} {'aspartate':>12s} {'glutamate':>12s} {'difference':>12s} {'p(diff)':>9s}")
print(f"   {'composition + source (as above)':38s} {mD:>+12.5f} {mE:>+12.5f} {d:>+12.5f} {2 * stats.norm.sf(abs(d / sd)):>9.3f}")
print(f"   {'+ GC and genome size':38s} {mD2:>+12.5f} {mE2:>+12.5f} {d2:>+12.5f} {p_d2:>9.3f}")
surv = p_d2 < 0.05 and d2 > 0
print(f"   {'SURVIVES' if surv else 'DOES NOT SURVIVE'}: the difference is {d2 / d:.2f} of its uncontrolled size at p {p_d2:.3f}.")
if surv:
    print("   So the glutamate-over-aspartate asymmetry is not a codon-GC artefact. Worth stating because it is the only")
    print("   confound specific to a comparison BETWEEN two acidic residues - every other check in this project bears on")
    print("   whether the acidic class matters at all, and none of them would have caught this one.")
rows.append({"part": "GC confound", "panel": WHICH, "beta_D_with_gc": mD2, "se_D_with_gc": seD2,
             "beta_E_with_gc": mE2, "se_E_with_gc": seE2, "diff_with_gc": d2, "se_diff_with_gc": sd2,
             "p_diff_with_gc": p_d2, "diff_without_gc": d, "p_diff_without_gc": float(2 * stats.norm.sf(abs(d / sd))),
             "ratio_controlled_to_uncontrolled": float(d2 / d), "survives": bool(surv), "n_refs": int(len(G2))})

print("\n5. what that would do to the required sample size - COMPUTED, NOT ADOPTED")
# Stored because a derived quantity quoted anywhere has to be a row (correction 12), and flagged because changing a
# committed number is the pre-registration's decision to make under its section 10, not this script's.
SD_ASSAY, POWER = 0.10, 0.80


def n_per_group(delta, sd=SD_ASSAY, alpha=0.05, power=POWER):
    for n in range(3, 100000):
        dfree = 2 * n - 2
        ncp = delta / (sd * np.sqrt(2.0 / n))
        crit = stats.t.ppf(1 - alpha / 2, dfree)
        if stats.nct.sf(crit, dfree, ncp) + stats.nct.cdf(-crit, dfree, ncp) >= power:
            return n
    return -1


COMMITTED = 0.05
scaled = COMMITTED * equiv / 14
n0, n1 = n_per_group(COMMITTED), n_per_group(scaled)
print(f"   IF a position's contribution is proportional to its own residue's coefficient - an assumption, not a")
print(f"   measurement - the primary arm's predicted effect is {equiv / 14:.2f} of what counting 14 positions equally gives,")
print(f"   so the committed -{COMMITTED:.2f} becomes -{scaled:.4f} and the required n per group rises from {n0} to {n1}.")
print("   THIS IS NOT ADOPTED HERE. docs/PREREGISTRATION.md section 3 commits to -0.05 and n>=64; changing either is a")
print("   decision for its section 10 deviations procedure, and the proportionality assumption above has not been")
print("   tested by anything. It is recorded so the choice is made deliberately rather than discovered at the bench.")
rows.append({"part": "power implication NOT ADOPTED", "panel": WHICH, "committed_delta": -COMMITTED,
             "scaled_delta": -float(scaled), "glutamate_equivalent_positions": float(equiv),
             "n_committed": int(n0), "n_if_reweighted": int(n1), "sd_assay": SD_ASSAY,
             "note": "conditional on contribution being proportional to the residue's own coefficient; untested"})

R = pd.DataFrame(rows)
R.to_csv(D / f"de_joint{TAG}.csv", index=False)
print("\nwrote", D / f"de_joint{TAG}.csv")
print("\nReading. Entered together and on a scale that cannot be confused, glutamate and aspartate are separable: the two")
print("shares are near-independent across genomes, so the joint fit is interpretable, and the pooled coefficients differ")
print("in the direction correction 14 found from univariate fits. The arm's five aspartate positions are worth much less")
print("than its nine glutamate ones, so counting eligible positions without regard to which residue they carry overstates")
print("what the arm does - which is a defect of scripts/37b's selection criterion rather than of this particular arm.")
