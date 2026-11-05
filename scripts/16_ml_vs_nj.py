"""Compare the ML (IQ-TREE) and NJ marker trees: agreement with taxonomy, Robinson-Foulds distance, and the Bacteria Kdp / Trk-Ktr
PGLS on the ML tree and on its 100 ultrafast-bootstrap trees. Usage: python3 scripts/16_ml_vs_nj.py [Bacteria|Archaea]   (the Kdp / Trk-Ktr regressions run for Bacteria only)"""
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
DOM = sys.argv[1] if len(sys.argv) > 1 else "Bacteria"
ML, NJ = ROOT / "results/phylo_ml", ROOT / "results/phylo_full"
ph_all = pd.read_csv(ROOT / "results/full/phenotypes_all.csv").drop_duplicates("ncbi_assembly").set_index("ncbi_assembly")


def splits(tree):
    tips = frozenset(t.name for t in tree.get_terminals()); out = set()
    for c in tree.get_nonterminals():
        s = frozenset(t.name for t in c.get_terminals())
        if 1 < len(s) < len(tips) - 1:
            out.add(min(s, tips - s, key=lambda x: (len(x), sorted(x))))
    return out


def mono(tree, members):
    members = set(members)
    if len(members) < 2:
        return None
    m = tree.common_ancestor([t for t in tree.get_terminals() if t.name in members])
    return {t.name for t in m.get_terminals()} == members


rows = []
trees = {"NJ": Phylo.read(NJ / f"{DOM}.nwk", "newick"), "ML": Phylo.read(ML / f"{DOM}.nwk", "newick")}
common = sorted(set(t.name for t in trees["NJ"].get_terminals()) & set(t.name for t in trees["ML"].get_terminals()))
for name, t in trees.items():
    tips = [x.name for x in t.get_terminals()]
    for rank in ("class", "family", "genus"):
        g = ph_all.loc[tips].groupby(rank).size(); g = g[g >= 3]
        ok = sum(bool(mono(t, ph_all.loc[tips][ph_all.loc[tips][rank] == k].index)) for k in g.index)
        rows.append({"tree": name, "rank": rank, "groups(>=3 tips)": len(g), "monophyletic": ok})
print(pd.DataFrame(rows).pivot(index="rank", columns="tree", values="monophyletic").assign(groups=pd.DataFrame(rows).groupby("rank")["groups(>=3 tips)"].first()).to_string())
# RF on the shared tips
def prune_to(t, keep):
    import copy
    t2 = copy.deepcopy(t)
    for x in [x for x in t2.get_terminals() if x.name not in keep]:
        t2.prune(x)
    return t2
a, b = splits(prune_to(trees["NJ"], set(common))), splits(prune_to(trees["ML"], set(common)))
print(f"\nshared tips {len(common)} | NJ splits {len(a)}, ML splits {len(b)}, shared {len(a & b)} | normalised Robinson-Foulds = {(len(a ^ b)) / (len(a) + len(b)):.3f}")
for lab, fam in ((("Alicyclobacillaceae", "Alicyclobacillaceae"),) if DOM == "Bacteria" else (("Sulfolobaceae", "Sulfolobaceae"),)):
    for name, t in trees.items():
        tips = [x.name for x in t.get_terminals()]; mem = [x for x in tips if ph_all.loc[x, "family"] == fam]
        print(f"{lab} monophyletic in {name}: {mono(t, mem)} ({len(mem)} tips)")

# PGLS on ML tree and UFBoot trees (Kdp / Trk-Ktr are bacterial questions)
if DOM != "Bacteria":
    sys.exit(0)
D = ROOT / "results/full"
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv").drop(columns=["n_genes"]), left_on="ncbi_assembly", right_on="genome_id"))
df = df[(df.domain == "Bacteria") & df.ph_opt.notna() & df.ncbi_assembly.isin(common)].reset_index(drop=True)
for o in ("thermobase", "acidb", "bacdive"):
    df[f"src_{o}"] = (df.origin == o).astype(float)
cov = ["ph_opt", "src_thermobase", "src_acidb", "src_bacdive"]
y = df.ph_min.to_numpy(float); z = lambda v: (np.asarray(v, float) - np.mean(v)) / np.std(v)
Xb = np.column_stack([np.ones(len(df))] + [df[c].to_numpy(float) for c in cov]); k = Xb.shape[1]
ZK, ZT = z(df.hm_rate_kdp), z(df.hm_rate_trk_ktr)


def fit(tree, feats):
    C = ph.vcv_fast(tree, list(df.ncbi_assembly))
    X = np.column_stack([Xb] + feats); r = ph.pgls(y, X, C)
    return r


out = []
for name in ("NJ", "ML"):
    t = trees[name]
    for lab, feats in (("Kdp alone", [ZK]), ("Trk/Ktr alone", [ZT]), ("Kdp + Trk/Ktr", [ZK, ZT])):
        r = fit(t, feats)
        for i, nm in enumerate(["Kdp" if "Kdp" in lab else "Trk/Ktr"] if len(feats) == 1 else ["Kdp", "Trk/Ktr"]):
            out.append({"tree": name, "model": lab, "feature": nm, "beta_per_sd": r["beta"][k + i], "p_two_sided": r["p"][k + i], "lambda": r["lambda"]})
R = pd.DataFrame(out); R.to_csv(D / "ml_vs_nj_pgls.csv", index=False)
print("\n" + R.round(5).to_string(index=False))
boots = list(Phylo.parse(ML / "Bacteria_bootstrap.nwk", "newick"))
keep = [t for t in boots]
B = []
for t in keep:
    tips = {x.name for x in t.get_terminals()}
    if not set(df.ncbi_assembly) <= tips:
        continue
    r = fit(t, [ZK, ZT]); B.append({"kdp_b": r["beta"][k], "kdp_p": r["p"][k], "trk_b": r["beta"][k + 1], "trk_p": r["p"][k + 1], "lambda": r["lambda"]})
B = pd.DataFrame(B); B.to_csv(D / "ml_ufboot_pgls.csv", index=False)
print(f"\n{len(B)} ML ultrafast-bootstrap trees (joint model): Kdp beta {B.kdp_b.quantile(.025):.3f} to {B.kdp_b.quantile(.975):.3f}, p median {B.kdp_p.median():.1e}, max {B.kdp_p.max():.1e}, share p<0.05 {(B.kdp_p < .05).mean():.2f}"
      f" | Trk/Ktr beta {B.trk_b.quantile(.025):.3f} to {B.trk_b.quantile(.975):.3f}, p median {B.trk_p.median():.1e}, max {B.trk_p.max():.1e}, share p<0.05 {(B.trk_p < .05).mean():.2f} | lambda {B['lambda'].min():.2f}-{B['lambda'].max():.2f}")
