"""Cross-checks on the Alicyclobacillaceae Kdp result: (1) PGLS with data-source/tier covariates, (2) phylogenetic independent
contrasts, (3) Fitch parsimony count of independent Kdp transitions and the residual change across them.
Usage: python3 scripts/10_phylo_crosschecks.py"""
import copy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import stats
from acidophile import phylo as ph

D = ROOT / "results/within_lineage_v2"; OUT = ROOT / "results/phylo"
df = (pd.read_csv(D / "lineage_phenotypes.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv"), left_on="ncbi_assembly", right_on="genome_id"))
g = df[df.lineage == "Alicyclobacillaceae"].reset_index(drop=True)
main = Phylo.read(OUT / "Bacteria.nwk", "newick")
boots = list(Phylo.parse(OUT / "Bacteria_bootstrap.nwk", "newick"))
g["z_kdp"] = (g.hm_rate_kdp - g.hm_rate_kdp.mean()) / g.hm_rate_kdp.std()
g["db_tier"] = (g.tier != "verified").astype(float)
g["is_thermobase"] = (g.origin == "thermobase").astype(float)
g["is_acidb"] = (g.origin == "acidb").astype(float)
res = []


def pgls_run(label, cols, tree=main):
    C = ph.vcv_from_tree(tree, list(g.ncbi_assembly))
    X = np.column_stack([np.ones(len(g))] + [g[c].to_numpy(float) for c in cols])
    r = ph.pgls(g.ph_min.to_numpy(float), X, C)
    k = cols.index("z_kdp") + 1
    res.append({"check": label, "n": len(g), "lambda": round(r["lambda"], 2), "beta_kdp_per_sd": round(r["beta"][k], 3),
                "p_one_sided_neg": round(float(stats.t.cdf(r["t"][k], r["n"] - X.shape[1])), 4)})


pgls_run("PGLS ph_opt + Kdp (reference)", ["ph_opt", "z_kdp"])
pgls_run("+ database-tier indicator", ["ph_opt", "db_tier", "z_kdp"])
pgls_run("+ source indicators (ThermoBase, AciDB vs literature)", ["ph_opt", "is_thermobase", "is_acidb", "z_kdp"])
pgls_run("+ database tier + growth temperature", ["ph_opt", "db_tier", "temp_opt_c", "z_kdp"]) if g.temp_opt_c.notna().all() else None
# (2) independent contrasts on the induced subtree
sub = copy.deepcopy(main)
for t in [t for t in sub.get_terminals() if t.name not in set(g.ncbi_assembly)]:
    sub.prune(t)
sub.root_with_outgroup  # noqa (tree is already rooted at the midpoint of the full tree)
vals = {r.ncbi_assembly: r for r in g.itertuples()}


def pic(clade, trait):
    """Felsenstein contrasts; returns (value, extra branch length, contrasts list, variances list)."""
    if clade.is_terminal():
        return getattr(vals[clade.name], trait), clade.branch_length or 0.0, [], []
    kids = [pic(c, trait) for c in clade.clades]
    while len(kids) > 2:   # resolve polytomy arbitrarily by merging the first two with zero length
        a, b = kids[0], kids[1]
        kids = [_merge(a, b)] + kids[2:]
    v, ex, cs, vs = _merge(kids[0], kids[1])
    return v, ex + (clade.branch_length or 0.0), cs, vs


def _merge(a, b):
    (va, ea, ca, sa), (vb, eb, cb, sb) = a, b
    la, lb = ea, eb
    tot = (la + lb) or 1e-9
    contrast = (va - vb) / np.sqrt(tot)
    node = (va * lb + vb * la) / tot
    return node, la * lb / tot, ca + cb + [contrast], sa + sb + [tot]


def contrasts(trait):
    return np.array(pic(sub.root, trait)[2])


Y, O, K = contrasts("ph_min"), contrasts("ph_opt"), contrasts("hm_rate_kdp")
X = np.column_stack([O, K])
beta, *_ = np.linalg.lstsq(X, Y, rcond=None)           # regression through the origin
r = Y - X @ beta; dfree = len(Y) - 2
cov = (r @ r / dfree) * np.linalg.inv(X.T @ X); t = beta[1] / np.sqrt(cov[1, 1])
res.append({"check": "independent contrasts: d(pH_min) ~ d(pH_opt) + d(Kdp rate)", "n": len(Y) + 1, "lambda": None,
            "beta_kdp_per_sd": round(float(beta[1]), 3), "p_one_sided_neg": round(float(stats.t.cdf(t, dfree)), 4)})
# (3) Fitch parsimony on Kdp presence
pres = {r.ncbi_assembly: {int(r.hm_cn_kdp > 0)} for r in g.itertuples()}
changes = 0


def fitch(clade):
    global changes
    if clade.is_terminal():
        return pres[clade.name]
    sets = [fitch(c) for c in clade.clades]
    inter = set.intersection(*sets)
    if inter:
        return inter
    changes += 1
    return set.union(*sets)


fitch(sub.root)
res.append({"check": f"Fitch parsimony: minimum independent Kdp gains/losses on the tree = {changes} (Kdp present in {int((g.hm_cn_kdp > 0).sum())} of {len(g)})",
            "n": len(g), "lambda": None, "beta_kdp_per_sd": None, "p_one_sided_neg": None})
R = pd.DataFrame(res); R.to_csv(OUT / "kdp_phylo_crosschecks.csv", index=False)
pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 120)
print(R.to_string(index=False))

# ---- additional checks (appended) ----
# (4) recurrence: among contrasts where Kdp differs, does the pH_min residual (after pH_opt) move opposite to Kdp?
from scipy.stats import binomtest
bO = float((O @ Y) / (O @ O)); rY = Y - bO * O
moving = K != 0
agree = int(((np.sign(rY) == -np.sign(K)) & moving).sum()); nmove = int(moving.sum())
extra = [{"check": f"contrasts where Kdp differs: pH_min residual moves opposite to Kdp in {agree} of {nmove} (one-sided sign test)", "n": nmove,
          "lambda": None, "beta_kdp_per_sd": None, "p_one_sided_neg": round(float(binomtest(agree, nmove, 0.5, alternative='greater').pvalue), 4)}]
# (5) split: original 13 strains vs the 13 added later, each with its own PGLS on the same tree
old = set(pd.read_csv(ROOT / "results/within_lineage/lineage_phenotypes.csv").ncbi_assembly)
full = g.copy()
for label, mask in (("PGLS, 13 strains from the first run only", full.ncbi_assembly.isin(old)), ("PGLS, 13 strains added later only", ~full.ncbi_assembly.isin(old))):
    g = full[mask].reset_index(drop=True)
    g["z_kdp"] = (g.hm_rate_kdp - g.hm_rate_kdp.mean()) / g.hm_rate_kdp.std()
    C = ph.vcv_from_tree(main, list(g.ncbi_assembly))
    X = np.column_stack([np.ones(len(g)), g.ph_opt, g.z_kdp]); r = ph.pgls(g.ph_min.to_numpy(float), X, C)
    extra.append({"check": label, "n": len(g), "lambda": round(r["lambda"], 2), "beta_kdp_per_sd": round(r["beta"][2], 3),
                  "p_one_sided_neg": round(float(stats.t.cdf(r["t"][2], r["n"] - 3)), 4)})
R = pd.concat([R, pd.DataFrame(extra)], ignore_index=True)
R.to_csv(OUT / "kdp_phylo_crosschecks.csv", index=False)
print("\n" + R.tail(3).to_string(index=False))

# (6) held-out lineages: leave-one-genus-out, does adding Kdp beat pH_opt alone?
g = full.copy()
rng2 = np.random.default_rng(3)
pred0, pred1 = pd.Series(np.nan, index=g.index), pd.Series(np.nan, index=g.index)
for genus in g.genus.unique():
    te = g.genus == genus; tr = ~te
    for cols, store in ((["ph_opt"], pred0), (["ph_opt", "hm_rate_kdp"], pred1)):
        A = np.column_stack([np.ones(tr.sum())] + [g.loc[tr, c] for c in cols]); b, *_ = np.linalg.lstsq(A, g.loc[tr, "ph_min"], rcond=None)
        store[te] = np.column_stack([np.ones(te.sum())] + [g.loc[te, c] for c in cols]) @ b
e0, e1 = (g.ph_min - pred0) ** 2, (g.ph_min - pred1) ** 2
gen_idx = {k: np.where(g.genus.values == k)[0] for k in g.genus.unique()}
keys = list(gen_idx)
d = []
for _ in range(3000):
    pick = np.concatenate([gen_idx[k] for k in rng2.choice(keys, len(keys))])
    d.append(np.sqrt(e1.values[pick].mean()) - np.sqrt(e0.values[pick].mean()))
dd = float(np.sqrt(e1.mean()) - np.sqrt(e0.mean()))
R = pd.concat([R, pd.DataFrame([{"check": f"leave-one-genus-out ({len(keys)} genera): RMSE opt+Kdp minus opt-only = {dd:.3f} (bootstrap over genera 95% CI {np.percentile(d,2.5):.3f} to {np.percentile(d,97.5):.3f}); negative = Kdp helps",
                                  "n": len(g), "lambda": None, "beta_kdp_per_sd": None, "p_one_sided_neg": None}])], ignore_index=True)
R.to_csv(OUT / "kdp_phylo_crosschecks.csv", index=False)
print("\n" + R.tail(1).check.iloc[0])
