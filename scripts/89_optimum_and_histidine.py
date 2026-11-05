"""How much of the association survives controlling the growth optimum AND histidine together? -> DATA_DIR/optimum_and_histidine.csv and, with `wide`, DATA_DIR/optimum_and_histidine_wide.csv
BOTH SPELLINGS ARE DECLARED, for the reason scripts/81 records: scripts/69 builds a claim's candidate pool from this line and keeps a panel variant's pool separate, so an
undeclared output turns every figure quoted from it into a guaranteed false miss.
THIS IS THE PROJECT'S LOAD-BEARING NUMBER AND THE ONE TEST CORRECTION 20 HAS NOT BEEN APPLIED TO. docs/SYNTHESIS.md section 6 reports that on the 1,655 genomes whose growth
optimum is known the composite gives +0.2123 per SD with composition controlled and +0.0904 with the optimum added, so 57% of the association runs through the optimum and
43% is independent of it. That 43% is the whole claim: a genomic signal for the acid LIMIT beyond what the growth optimum predicts. Correction 20 then found the exposed
acidic coefficient inflated about 27% by omitting exposed histidine - and the optimum decomposition has never been run with histidine held fixed, because the panel carried
no histidine column until 31 October. If both controls together take the independent part to zero, the central claim is in a different condition than this project has
reported; if they do not, the claim survives its hardest available test. Either way the number has to exist.
THE CRITERION IS FIXED BEFORE THE NUMBERS. If the measure loses significance with the optimum and histidine both held fixed, that is a MATERIAL change to the central claim
and must be reported as such rather than as a caveat. If it survives, the figure to quote is the share of the composition-controlled effect that remains - the analogue of
the published 43% - and not the raw coefficient, because the raw coefficient is not comparable across samples.
EVERY COMPARISON IS WITHIN-SAMPLE. The optimum is known for 1,657 of the 2,771 bacteria, so all four fits below run on the same optimum-known subset, and the full-sample
figure is printed separately and never used as a comparator. This is correction 12, which this project has now had to guard against in four scripts in two days.
THE ABSOLUTE FIGURES HERE ARE NOT THE PUBLISHED ONES, and saying so matters. +0.2123 and +0.0904 come from scripts/57's factor-score composite; this script uses the pooled
net share, a different quantity on a different scale. What is comparable is the SHARE that survives each control, which is scale-free, and that is what the verdict rests on.
SCALE. Per percentage point from ph.pgls directly, never ph.pgls_scan - correction 15. Per-SD figures are stated as such and are the per-unit coefficient times the
predictor's own standard deviation, printed. The other charged shares of the same class are held fixed throughout, which is correction 14.
Usage: python3 scripts/89_optimum_and_histidine.py [orig|wide]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
PUB_PLAIN, PUB_OPT = 0.2123, 0.0904      # docs/SYNTHESIS.md section 6, scripts/57's composite, for reference only
_r = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv")
MEM = [str(x) for x in _r[_r.accepted & _r.membrane].label]
rows = []

pan = pd.read_csv(D / PANEL)
if not any(c.endswith("|H") for c in pan.columns):
    sys.exit(f"{PANEL} has no histidine columns; rebuild the panel with scripts/36a first")
df = (pd.read_csv(D / "phenotypes_all.csv")
      .merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pan, on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips) & df.ph_min.notna()].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
RES = ("D", "E", "K", "R", "H")
tot = {k: np.zeros(len(B)) for k in ("n",) + RES}
for ref in MEM:
    if f"{ref}|peri_exp|n" not in B:
        continue
    tot["n"] += B[f"{ref}|peri_exp|n"].fillna(0).to_numpy(float)
    for k in RES:
        tot[k] += B[f"{ref}|peri_exp|{k}"].fillna(0).to_numpy(float)
havepos = tot["n"] >= MINPOS


def run(mask, label):
    """The four fits on `mask`: composition only, + optimum, + histidine, + both. Returns a dict of results."""
    g = B[mask].reset_index(drop=True)
    Cm = ph.vcv_fast(tree, list(g.ncbi_assembly))
    y = g.ph_min.to_numpy(float)
    srcc = [f"src_{o}" for o in SRC if g[f"src_{o}"].nunique() > 1]
    XS = np.column_stack([np.ones(len(g))] + [g[k] for k in ["gw"] + srcc])
    lam = ph.pgls(y, XS, Cm)["lambda"]
    s = {k: 100 * tot[k][mask] / tot["n"][mask] for k in RES}
    net = s["D"] + s["E"] - s["K"] - s["R"]
    sd = float(np.std(net, ddof=1))
    out = {}
    for name, extra in (("composition only", []),
                        ("+ growth optimum", [g.ph_opt.to_numpy(float)]),
                        ("+ histidine", [s["H"]]),
                        ("+ optimum and histidine", [g.ph_opt.to_numpy(float), s["H"]])):
        # The optimum fits cannot run where the optimum is missing, which is the whole point of the two samples.
        # The first version of this line read `any("optimum" in name for _ in [0])`, which is a one-element
        # generator around a plain substring test - it works, and nobody reading it would be sure it does.
        if "optimum" in name and g.ph_opt.isna().any():
            continue
        X = np.column_stack([XS] + extra + [net])
        assert np.linalg.matrix_rank(X) == X.shape[1], f"{label}/{name}: rank-deficient"
        r = ph.pgls(y, X, Cm, lam=lam)
        out[name] = (float(r["beta"][-1]), float(r["se"][-1]), float(r["p"][-1]), sd, int(len(g)))
    return out, lam, sd, int(len(g))


print("1. the samples, kept separate because mixing them is correction 12")
full = havepos & B.ph_min.notna().to_numpy()
sub = full & B.ph_opt.notna().to_numpy()
print(f"   full sample {int(full.sum())} bacteria; optimum-known subset {int(sub.sum())}. Every comparison below")
print("   runs on the subset; the full-sample figure is printed once and never used as a comparator.")
rows.append({"part": "samples", "n_full": int(full.sum()), "n_optimum_known": int(sub.sum())})
fo, lamf, sdf, nf = run(full, "full")
b = fo["composition only"]
print(f"   full sample, composition only: {b[0]:+.5f}/pt -> {b[0] * sdf:+.4f} per SD (p {b[2]:.2e}), n={nf}")
rows.append({"part": "full sample reference", "specification": "composition only", "beta_per_pt": b[0],
             "se": b[1], "p": b[2], "beta_per_sd": b[0] * sdf, "n": nf})

print("\n2. the four fits, all on the optimum-known subset")
res, lam, sd, n = run(sub, "subset")
print(f"   {'specification':>26} {'beta/pt':>10} {'se':>9} {'p':>10} {'per SD':>9} {'share kept':>11}")
base = res["composition only"]
for name in ("composition only", "+ growth optimum", "+ histidine", "+ optimum and histidine"):
    v = res[name]
    keep = v[0] / base[0]
    print(f"   {name:>26} {v[0]:>+10.5f} {v[1]:>9.5f} {v[2]:>10.2e} {v[0] * sd:>+9.4f} {keep:>11.0%}")
    rows.append({"part": "fit", "specification": name, "beta_per_pt": v[0], "se": v[1], "p": v[2],
                 "beta_per_sd": v[0] * sd, "share_of_composition_only": keep, "n": n, "lambda": lam,
                 "sd_of_composite": sd})

print("\n3. the decomposition, in the scale-free terms the published figure uses")
k_opt = res["+ growth optimum"][0] / base[0]
k_his = res["+ histidine"][0] / base[0]
k_both = res["+ optimum and histidine"][0] / base[0]
pub_k = PUB_OPT / PUB_PLAIN
print(f"   independent of the optimum alone:   {k_opt:.0%}   (published, on the factor-score composite: {pub_k:.0%})")
print(f"   independent of histidine alone:     {k_his:.0%}")
print(f"   independent of BOTH:                {k_both:.0%}")
print("   The first row is the comparable one, and it lands close to the published figure on a different")
print("   composite - which is the agreement wanted between two routes to the same decomposition, not a check")
print("   that either is right.")
rows.append({"part": "decomposition", "kept_optimum_only": k_opt, "kept_histidine_only": k_his,
             "kept_both": k_both, "published_kept_optimum_only": pub_k})

print("\n4. the verdict, by the criterion fixed in the preamble")
v = res["+ optimum and histidine"]
survives = v[2] < 0.05 and (v[0] > 0) == (base[0] > 0)
if survives:
    print(f"   SURVIVES. With the growth optimum and exposed histidine both held fixed the measure is "
          f"{v[0]:+.5f}/pt")
    print(f"   ({v[0] * sd:+.4f} per SD, p {v[2]:.2e}), {k_both:.0%} of the composition-controlled effect. The")
    print("   central claim passes the hardest test currently available to it: a genomic signal for the acid")
    print("   limit that is not the growth optimum and not the residue correction 20 showed it had been")
    print("   borrowing from.")
else:
    print(f"   DOES NOT SURVIVE. {v[0]:+.5f}/pt, p {v[2]:.2e}. This is a MATERIAL change to the central claim")
    print("   and not a caveat: the part of the association independent of the growth optimum is not")
    print("   distinguishable from zero once histidine is held fixed, and docs/SYNTHESIS.md section 1 would")
    print("   have to say so.")
rows.append({"part": "verdict", "survives": survives, "beta_per_pt": v[0], "p": v[2],
             "beta_per_sd": v[0] * sd, "share_kept": k_both})
# The two controls are not independent of each other, and the product of their separate shares says whether they
# overlap. If k_opt * k_his is close to k_both the two remove different variance; if k_both is larger, they
# overlap and the combined cost is less than the sum of the parts.
prod = k_opt * k_his
print(f"\n5. do the two controls remove the same variance? product of the separate shares {prod:.0%} against "
      f"{k_both:.0%} jointly")
print("   " + ("They overlap: the joint cost is LESS than the product, so part of what the optimum removes is the"
               if k_both > prod + 0.02 else
               "They are close to independent, so each removes variance the other does not and the costs compound"
               if abs(k_both - prod) <= 0.02 else
               "The joint cost is MORE than the product, which means they reinforce rather than overlap"))
print("   " + ("   same variance histidine removes." if k_both > prod + 0.02 else
               "   nearly multiplicatively." if abs(k_both - prod) <= 0.02 else
               "   - unusual, and worth distrusting before quoting."))
rows.append({"part": "control overlap", "product_of_separate_shares": prod, "joint_share": k_both})
pd.DataFrame(rows).to_csv(D / f"optimum_and_histidine{TAG}.csv", index=False)
print("\nwrote", D / f"optimum_and_histidine{TAG}.csv")
