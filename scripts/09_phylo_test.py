"""Phylogeny-aware test within each lineage: PGLS (Pagel's lambda by ML) of pH_min on pH_opt (+ growth temperature) + one feature,
over the main marker tree and every bootstrap tree. Kdp is a pre-specified one-sided hypothesis (negative coefficient); the other
features are an exploratory screen with Benjamini-Hochberg correction. Usage: python3 scripts/09_phylo_test.py"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import stats
from acidophile import phylo as ph
from acidophile.within_lineage import ALL_FEATURES

D = ROOT / "results/within_lineage_v2"; OUT = ROOT / "results/phylo"
df = (pd.read_csv(D / "lineage_phenotypes.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv"), left_on="ncbi_assembly", right_on="genome_id"))
df["kdp_present"] = (df.hm_cn_kdp > 0).astype(float)
DOM = {"Alicyclobacillaceae": "Bacteria", "Sulfolobales": "Archaea"}


def trees_for(dom):
    main = Phylo.read(OUT / f"{dom}.nwk", "newick")
    boots = list(Phylo.parse(OUT / f"{dom}_bootstrap.nwk", "newick"))
    return main, boots


def fit(tree, g, feat, covs, lam=None):
    tips = list(g.ncbi_assembly)
    C = ph.vcv_from_tree(tree, tips)
    X = np.column_stack([np.ones(len(g))] + [g[c].to_numpy(float) for c in covs] + ([(g[feat] - g[feat].mean()) / g[feat].std()] if feat else []))
    return ph.pgls(g.ph_min.to_numpy(float), X, C, lam)


def one_sided_neg(r, k):
    t, df_ = r["t"][k], r["n"] - len(r["beta"])
    return float(stats.t.cdf(t, df_))


rows, screen = [], []
for lin, g0 in df.groupby("lineage"):
    main, boots = trees_for(DOM[lin])
    for subset, g in (("all", g0), ("verified_only", g0[g0.tier == "verified"])):
        for covs in (["ph_opt"], ["ph_opt", "temp_opt_c"]):
            g = g.dropna(subset=covs).reset_index(drop=True)
            if len(g) < 8:
                continue
            base = fit(main, g, None, covs)
            base0 = fit(main, g, None, covs, lam=0.0)
            lrt = 2 * (base["loglik"] - base0["loglik"])
            rows.append({"lineage": lin, "subset": subset, "covariates": "+".join(covs), "analysis": "baseline (no feature)", "n": len(g),
                         "lambda": round(base["lambda"], 2), "lrt_lambda_vs_0_p": round(float(stats.chi2.sf(max(lrt, 0), 1)), 3)})
            for feat in ("hm_rate_kdp", "kdp_present"):
                if g[feat].nunique() < 2:
                    continue
                k = 1 + len(covs)
                rm = fit(main, g, feat, covs)
                pb = np.array([one_sided_neg(fit(b, g, feat, covs), k) for b in boots])
                bb = np.array([fit(b, g, feat, covs)["beta"][k] for b in boots])
                rows.append({"lineage": lin, "subset": subset, "covariates": "+".join(covs), "analysis": f"PGLS {feat} (one-sided, negative)", "n": len(g),
                             "lambda": round(rm["lambda"], 2), "beta_per_sd": round(float(rm["beta"][k]), 3), "se": round(float(rm["se"][k]), 3),
                             "p_main_tree": round(one_sided_neg(rm, k), 4), "p_median_bootstrap": round(float(np.median(pb)), 4),
                             "frac_boot_trees_p<0.05": round(float((pb < 0.05).mean()), 2), "beta_boot_2.5-97.5%": f"{np.percentile(bb, 2.5):.2f} to {np.percentile(bb, 97.5):.2f}"})
            if subset == "all" and covs == ["ph_opt"]:      # exploratory screen of all features, two-sided
                feats = [f for f in ALL_FEATURES if f in g and g[f].nunique() > 1]
                for f in feats:
                    rm = fit(main, g, f, covs); k = 1 + len(covs)
                    pb = np.median([fit(b, g, f, covs)["p"][k] for b in boots[:40]])
                    screen.append({"lineage": lin, "feature": f, "beta_per_sd": rm["beta"][k], "p_two_sided_main": rm["p"][k], "p_median_boot40": pb, "lambda": rm["lambda"]})
R = pd.DataFrame(rows); R.to_csv(OUT / "pgls_summary.csv", index=False)
S = pd.DataFrame(screen)
if len(S):
    S["q_BH"] = S.groupby("lineage").p_two_sided_main.transform(lambda p: np.minimum.accumulate((p.sort_values(ascending=False) * len(p) / (len(p) - np.arange(len(p)))).values)[::-1].tolist() and 0) if False else np.nan
    for lin, idx in S.groupby("lineage").groups.items():
        p = S.loc[idx, "p_two_sided_main"].to_numpy(); o = np.argsort(p); q = np.empty_like(p)
        q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
        S.loc[idx, "q_BH"] = np.minimum(q, 1)
    S.sort_values(["lineage", "p_two_sided_main"]).to_csv(OUT / "pgls_feature_screen.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(R.to_string(index=False))
print("\nTop of exploratory screen:"); print(S.sort_values("p_two_sided_main").groupby("lineage").head(4).round(3).to_string(index=False))
