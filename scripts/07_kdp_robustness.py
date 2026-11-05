"""Robustness battery for the Kdp lead in Alicyclobacillaceae (hypothesis fixed in advance: more Kdp -> lower pH_min than
ph_opt predicts, i.e. negative rank correlation with the residual). Writes results/within_lineage_v2/kdp_robustness.csv."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
from acidophile.within_lineage import _ols_resid, _spearman

D = ROOT / "results/within_lineage_v2"
df = (pd.read_csv(D / "lineage_phenotypes.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv"), left_on="ncbi_assembly", right_on="genome_id"))
g = df[df.lineage == "Alicyclobacillaceae"].reset_index(drop=True)
old = set(pd.read_csv(ROOT / "results/within_lineage/lineage_phenotypes.csv").ncbi_assembly)
rng = np.random.default_rng(7)
N = 20000
rows = []


def resid(sub, covs):
    return _ols_resid(sub.ph_min.to_numpy(float), sub[list(covs)].to_numpy(float))[0]


def test(name, sub, covs=("ph_opt",), strata=None, feat="hm_rate_kdp"):
    sub = sub.dropna(subset=list(covs) + [feat]).reset_index(drop=True)
    r = resid(sub, covs); x = sub[feat].to_numpy(float)
    obs = _spearman(r, x)
    if strata is None:
        null = np.array([_spearman(rng.permutation(r), x) for _ in range(N)])
    else:   # permute residuals only within strata (e.g. genus): clade-level differences cannot create the association
        s = sub[strata].to_numpy(); idx = {k: np.where(s == k)[0] for k in np.unique(s)}
        null = np.empty(N)
        for i in range(N):
            rp = r.copy()
            for ix in idx.values():
                rp[ix] = rng.permutation(r[ix])
            null[i] = _spearman(rp, x)
    rows.append({"analysis": name, "n": len(sub), "rho": round(obs, 3), "p_one_sided_neg": round(float((null <= obs + 1e-12).mean()), 4)})


test("all 26, ph_opt", g)
test("new 13 only, ph_opt", g[~g.ncbi_assembly.isin(old)])
test("earlier 13 only, ph_opt", g[g.ncbi_assembly.isin(old)])
test("all 26, ph_opt + growth temperature", g, covs=("ph_opt", "temp_opt_c"))
test("all 26, ph_opt + genome size", g, covs=("ph_opt", "genome_size_bp"))
test("all 26, permuted within genus", g, strata="genus")
test("all 26, permuted within genus, + temperature", g, covs=("ph_opt", "temp_opt_c"), strata="genus")
test("verified-only 11, ph_opt", g[g.tier == "verified"])
test("database-only 15, ph_opt", g[g.tier != "verified"])
test("excluding assumption-based rows", g[~g.species.isin(["Alicyclobacillus mengziensis", "Alicyclobacillus curvatus"])])
test("complete or chromosome-level genomes only", g[g.assembly_level.isin(["Complete Genome", "Chromosome"])])
# presence/absence version (Kdp >= 1 copy) as a Mann-Whitney-style check via same machinery
g2 = g.assign(kdp_present=(g.hm_cn_kdp > 0).astype(float))
test("presence/absence of Kdp, all 26", g2, feat="kdp_present")
# leave-one-out influence
r = resid(g, ("ph_opt",)); x = g.hm_rate_kdp.to_numpy(float)
loo = [_spearman(np.delete(r, i), np.delete(x, i)) for i in range(len(g))]
rows.append({"analysis": f"leave-one-out range (drop each strain); most influential: {g.species[int(np.argmax(loo))]}", "n": len(g) - 1,
             "rho": f"{min(loo):.2f} to {max(loo):.2f}", "p_one_sided_neg": None})
out = pd.DataFrame(rows)
out.to_csv(D / "kdp_robustness.csv", index=False)
print(out.to_string(index=False))
print("\ngenera:", g.genus.value_counts().to_dict(), "| Kdp copies by genus:", g.groupby("genus").hm_cn_kdp.mean().round(2).to_dict())
print("corr(kdp rate, genome size)=%.2f, corr(kdp rate, temp)=%.2f" % (g.hm_rate_kdp.corr(g.genome_size_bp, method="spearman"), g.hm_rate_kdp.corr(g.temp_opt_c, method="spearman")))
