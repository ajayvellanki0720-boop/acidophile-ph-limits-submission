"""Phylogeny-aware screen restricted to species whose minima come from species descriptions (confidence medium/high; all literature,
so no data-source mixing), on the ML trees. Compares with the all-species estimates and calibrates against tree-simulated null features.
Usage: python3 scripts/20_verified_only_phylo.py [n_sims]"""
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
from acidophile.within_lineage import ALL_FEATURES

warnings.filterwarnings("ignore")
NS = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
D, T = ROOT / os.environ.get("DATA_DIR", "results/full"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_ml")
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv").drop(columns=["n_genes"]), left_on="ncbi_assembly", right_on="genome_id"))
df = df[df.ph_opt.notna()].reset_index(drop=True)
df["verified"] = df.confidence.isin(["medium", "high"])
z = lambda v: (np.asarray(v, float) - np.mean(v)) / np.std(v)
rng = np.random.default_rng(21)
out, scr = [], []
for dom in ("Bacteria", "Archaea"):
    tree = Phylo.read(T / f"{dom}.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
    for label, mask in (("verified only", df.verified), ("all species", np.ones(len(df), bool))):
        g = df[(df.domain == dom) & mask & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
        src = [f"src_{o}" for o in ("thermobase", "acidb", "bacdive") if (g.origin == o.replace("src_", "")).any() and g.origin.nunique() > 1]
        for o in ("thermobase", "acidb", "bacdive"):
            g[f"src_{o}"] = (g.origin == o).astype(float)
        cov = ["ph_opt"] + [c for c in ("src_thermobase", "src_acidb", "src_bacdive") if g[c].nunique() > 1]
        C = ph.vcv_fast(tree, list(g.ncbi_assembly)); y = g.ph_min.to_numpy(float)
        Xb = np.column_stack([np.ones(len(g))] + [g[c].to_numpy(float) for c in cov]); k = Xb.shape[1]
        base = ph.pgls(y, Xb, C)
        if dom == "Bacteria":
            ZK, ZT = z(g.hm_rate_kdp), z(g.hm_rate_trk_ktr)
            for lab, feats, names in (("Kdp alone", [ZK], ["Kdp"]), ("Trk/Ktr alone", [ZT], ["Trk/Ktr"]), ("Kdp + Trk/Ktr", [ZK, ZT], ["Kdp", "Trk/Ktr"])):
                r = ph.pgls(y, np.column_stack([Xb] + feats), C)
                for i, nm in enumerate(names):
                    out.append({"domain": dom, "subset": label, "n": len(g), "model": lab, "feature": nm, "beta_per_sd": r["beta"][k + i], "se": r["se"][k + i], "p": r["p"][k + i], "lambda": r["lambda"]})
            # calibration for the joint model on this subset
            Lc = np.linalg.cholesky(C / np.mean(np.diag(C)) + 1e-8 * np.eye(len(g)))
            for nm, x in (("Kdp", g.hm_rate_kdp.to_numpy(float)), ("Trk/Ktr", g.hm_rate_trk_ktr.to_numpy(float))):
                srt = np.sort(x); ts = []
                for _ in range(NS):
                    xs = srt[np.argsort(np.argsort(Lc @ rng.normal(size=len(g))))]
                    ts.append(ph.pgls(y, np.column_stack([Xb, z(xs)]), C)["t"][-1])
                tobs = ph.pgls(y, np.column_stack([Xb, z(x)]), C)["t"][-1]
                out.append({"domain": dom, "subset": label, "n": len(g), "model": f"calibration vs {NS} tree-simulated features", "feature": nm, "beta_per_sd": np.nan, "se": np.nan,
                            "p": (((np.abs(ts) >= abs(tobs)).sum() + 1) / (NS + 1)), "lambda": np.nan})
        for f in [f for f in ALL_FEATURES + ["gc"] if f in g and g[f].nunique() > 2]:
            r = ph.pgls(y, np.column_stack([Xb, z(g[f])]), C)
            scr.append({"domain": dom, "subset": label, "feature": f, "n": len(g), "beta_per_sd": r["beta"][k], "p_two_sided": r["p"][k]})
O = pd.DataFrame(out); O.to_csv(D / "verified_only_pgls.csv", index=False)
S = pd.DataFrame(scr)
for key, idx in S.groupby(["domain", "subset"]).groups.items():
    p = S.loc[idx, "p_two_sided"].to_numpy(); o = np.argsort(p); q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]; S.loc[idx, "q_BH"] = np.minimum(q, 1)
S.sort_values(["domain", "subset", "p_two_sided"]).to_csv(D / "verified_only_screen.csv", index=False)
pd.set_option("display.width", 220)
print(O.round(5).to_string(index=False))
for (dom, sub), g in S.groupby(["domain", "subset"]):
    print(f"\n{dom} / {sub}: top features"); print(g.sort_values("p_two_sided").head(4).round(3).to_string(index=False))
