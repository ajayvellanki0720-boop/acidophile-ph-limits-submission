"""Does the Kdp (negative) / Trk-Ktr (positive) association recur in independent lineages? PGLS fitted separately inside each bacterial
phylum with enough species (ML tree), then random-effects (DerSimonian-Laird) meta-analysis across phyla, heterogeneity, sign test,
and leave-one-phylum-out meta-analysis. Usage: python3 scripts/17_recurrence_by_phylum.py [PHYLO_DIR] [rank=phylum|class] [min_n=12]"""
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
PD = ROOT / (sys.argv[1] if len(sys.argv) > 1 else "results/phylo_ml")
RANK = sys.argv[2] if len(sys.argv) > 2 else "phylum"
MIN_N = int(sys.argv[3]) if len(sys.argv) > 3 else 12
D = ROOT / "results/full"
tree = Phylo.read(PD / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv").drop(columns=["n_genes"]), left_on="ncbi_assembly", right_on="genome_id"))
df = df[(df.domain == "Bacteria") & df.ph_opt.notna() & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
for o in ("thermobase", "acidb", "bacdive"):
    df[f"src_{o}"] = (df.origin == o).astype(float)
z = lambda v: (np.asarray(v, float) - np.mean(v)) / np.std(v)
rows = []
for phy, g in df.groupby(RANK):
    g = g.reset_index(drop=True)
    if len(g) < MIN_N or g.hm_rate_kdp.nunique() < 3 or g.hm_rate_trk_ktr.nunique() < 3:
        continue
    cov = ["ph_opt"] + [c for c in ("src_thermobase", "src_acidb", "src_bacdive") if g[c].nunique() > 1]
    C = ph.vcv_from_tree(tree, list(g.ncbi_assembly))
    X = np.column_stack([np.ones(len(g))] + [g[c].to_numpy(float) for c in cov] + [z(g.hm_rate_kdp), z(g.hm_rate_trk_ktr)])
    r = ph.pgls(g.ph_min.to_numpy(float), X, C); k = X.shape[1] - 2
    rows.append({RANK: phy, "n": len(g), "n_families": g.family.nunique(), "kdp_prev": round((g.hm_cn_kdp > 0).mean(), 2), "trk_prev": round((g.hm_cn_trk_ktr > 0).mean(), 2),
                 "lambda": round(r["lambda"], 2), "kdp_b": r["beta"][k], "kdp_se": r["se"][k], "kdp_p": r["p"][k], "trk_b": r["beta"][k + 1], "trk_se": r["se"][k + 1], "trk_p": r["p"][k + 1]})
R = pd.DataFrame(rows).sort_values("n", ascending=False)
R.to_csv(D / f"recurrence_by_{RANK}.csv", index=False)
pd.set_option("display.width", 220)
print(R.round(3).to_string(index=False))


def dl(b, se):
    w = 1 / se ** 2; fe = (w * b).sum() / w.sum(); Q = (w * (b - fe) ** 2).sum(); df_ = len(b) - 1
    tau2 = max(0.0, (Q - df_) / (w.sum() - (w ** 2).sum() / w.sum())) if df_ > 0 else 0.0
    wr = 1 / (se ** 2 + tau2); mu = (wr * b).sum() / wr.sum(); sem = np.sqrt(1 / wr.sum())
    return mu, sem, 2 * stats.norm.sf(abs(mu / sem)), (max(0.0, (Q - df_) / Q) if Q > 0 else 0.0), stats.chi2.sf(Q, df_) if df_ > 0 else np.nan


for feat, sign in (("kdp", -1), ("trk", +1)):
    b, se = R[f"{feat}_b"].to_numpy(), R[f"{feat}_se"].to_numpy()
    mu, sem, p, i2, pq = dl(b, se)
    agree = int((np.sign(b) == sign).sum())
    print(f"\n{feat.upper()}: random-effects across {len(R)} {RANK} groups: {mu:+.3f} ± {sem:.3f} (p {p:.1e}); I² {i2:.0%} (Q p {pq:.2f}); expected sign in {agree}/{len(R)} groups (sign test p {stats.binomtest(agree, len(R), .5, alternative='greater').pvalue:.3f})")
    loo = []
    for i in range(len(R)):
        m = np.arange(len(R)) != i
        mu_i, _, p_i, _, _ = dl(b[m], se[m]); loo.append((R[RANK].iloc[i], mu_i, p_i))
    print("  leave-one-phylum-out pooled estimate: " + "; ".join(f"-{n[:12]} {m:+.3f} (p {pp:.0e})" for n, m, pp in loo))
