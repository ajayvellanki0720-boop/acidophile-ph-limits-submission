"""Close-relative contrasts for the histidine / Kdp association, by lineage.
 1. Descriptives per lineage (spread of pH_min, histidine, Kdp).
 2. Phylogenetic independent contrasts (PIC) within each lineage: d(pH_min) ~ d(pH_opt) [+ source indicators] + d(feature), regression through the origin; all contrasts
    and cherry-only contrasts (sister tips = closest relatives, least affected by deep structure). Features are standardised within the lineage so betas are pH units per SD.
 3. Power: the carrier phyla are repeatedly subsampled to the size of the other lineages; the share of subsamples that detect the effect (one-sided p<0.05) is the power
    available at that sample size if the carrier-phylum effect were real there.
 4. Interaction: PGLS on all bacteria with feature x (non-carrier phylum) interaction.
Usage: python3 scripts/28_close_relative_contrasts.py [n_subsamples]   (env DATA_DIR=results/broad, PHYLO_DIR=results/phylo_broad)"""
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
NSUB = int(sys.argv[1]) if len(sys.argv) > 1 else 300
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad")
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv").drop(columns=["n_genes"]), left_on="ncbi_assembly", right_on="genome_id"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips_all = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ph_opt.notna() & df.ncbi_assembly.isin(tips_all)].reset_index(drop=True)
for o in ("thermobase", "acidb", "bacdive", "literature_auto"):
    B[f"src_{o}"] = (B.origin == o).astype(float)
z = lambda v: (np.asarray(v, float) - np.mean(v)) / (np.std(v) or 1.0)
FEATS = {"His": "frac_His", "Kdp": "hm_rate_kdp"}
CARRIERS, OTHERS = ["Bacillota", "Pseudomonadota"], ["Actinomycetota", "Bacteroidota"]
rng = np.random.default_rng(77)
SKIPPED = []


def pic_fit(g, feats, cherry_only=False):
    """Regress contrasts of pH_min on contrasts of pH_opt, source indicators that vary, and the standardised features. Returns dict feature -> (beta, se, p), n_contrasts."""
    cols = ["ph_opt"] + [c for c in ("src_thermobase", "src_acidb", "src_bacdive") if g[c].nunique() > 1]     # literature_auto is the reference: no intercept in a contrast regression, so a full set of indicators would sum to a constant (zero contrasts)
    M = np.column_stack([g.ph_min.to_numpy(float)] + [g[c].to_numpy(float) for c in cols] + [z(g[FEATS[f]]) for f in feats])
    c, ch = ph.pic_contrasts(tree, list(g.ncbi_assembly), M)
    if cherry_only:
        c = c[ch]
    if len(c) < len(cols) + len(feats) + 5:
        return None, len(c)
    X = c[:, 1:]
    keep = np.where(X.std(axis=0) > 1e-12)[0]                     # a regressor with all-zero contrasts carries no information in this subset
    if len(keep) < X.shape[1] or np.linalg.matrix_rank(X) < X.shape[1]:
        SKIPPED.append((len(g), feats, cherry_only, "rank-deficient: " + ",".join(str(i) for i in range(X.shape[1]) if i not in keep)))
        return None, len(c)
    b, se, p = ph.regress_through_origin(c[:, 0], X)
    k0 = 1 + len(cols) - 1
    return {f: (b[k0 + i], se[k0 + i], p[k0 + i]) for i, f in enumerate(feats)}, len(c)


rows = []
print("1. lineage descriptives")
desc = []
for ph_name in CARRIERS + OTHERS:
    g = B[B.phylum == ph_name]
    desc.append({"phylum": ph_name, "n": len(g), "families": g.family.nunique(), "pH_min sd": g.ph_min.std(), "pH_opt sd": g.ph_opt.std(), "His sd (pts)": 100 * g.frac_His.std(), "Kdp prevalence": (g.hm_cn_kdp > 0).mean(),
                 "share literature_auto": (g.origin == "literature_auto").mean()})
print(pd.DataFrame(desc).round(3).to_string(index=False))
print("\n2. independent contrasts within lineage (beta = pH_min units per SD of the feature, from contrasts)")
for name in CARRIERS + OTHERS + ["Actinomycetes (class)", "Flavobacteriia (class)", "Cytophagia (class)"]:
    if name.endswith("(class)"):
        g = B[B["class"] == name.split(" ")[0]]
    else:
        g = B[B.phylum == name]
    for cherry in (False, True):
        for feats in (["His"], ["Kdp"], ["His", "Kdp"]):
            res, nc = pic_fit(g.reset_index(drop=True), feats, cherry)
            if res is None:
                continue
            for f, (b, se, p) in res.items():
                rows.append({"lineage": name, "n_species": len(g), "contrasts": "cherries only" if cherry else "all", "n_contrasts": nc, "model": "+".join(feats), "feature": f, "beta_per_sd": b, "se": se, "p_two_sided": p})
R = pd.DataFrame(rows); R.to_csv(D / "close_relative_contrasts.csv", index=False)
if SKIPPED:
    print(f"   [{len(SKIPPED)} fits skipped as rank-deficient in a small lineage or cherry subset; first: {SKIPPED[0]}]")
show = R[(R.model == "His+Kdp")].pivot_table(index=["lineage", "n_species", "contrasts", "n_contrasts"], columns="feature", values=["beta_per_sd", "se", "p_two_sided"]).round(3)
print(show.to_string())

print(f"\n3. power: carrier phyla subsampled to the size of the other lineages ({NSUB} subsamples; one-sided p<0.05 in the negative direction, His+Kdp model)")
obs = R[(R.model == "His+Kdp") & (R.contrasts == "all")].set_index(["lineage", "feature"]).beta_per_sd
pw = []
for carrier in CARRIERS:
    full = B[B.phylum == carrier].reset_index(drop=True)
    for target in OTHERS + ["Actinomycetes (class)"]:
        n = int((B["class"] == "Actinomycetes").sum()) if target.endswith("(class)") else int((B.phylum == target).sum())
        res = {"His": [], "Kdp": []}
        for _ in range(NSUB):
            g = full.sample(n=min(n, len(full)), random_state=int(rng.integers(1e9))).reset_index(drop=True)
            r, _ = pic_fit(g, ["His", "Kdp"])
            if r:
                for f in res:
                    b, se, p = r[f]; res[f].append((b, p / 2 if b < 0 else 1 - p / 2))
        for f in res:
            a = np.array(res[f]); ob = obs.get((target, f), np.nan)
            pw.append({"source_phylum": carrier, "subsampled_to_size_of": target, "n": n, "feature": f, "median_beta": np.median(a[:, 0]), "power_onesided_p05": float((a[:, 1] < 0.05).mean()),
                       "observed_target_beta": ob, "share_of_subsamples_as_weak_as_observed": float((a[:, 0] >= ob).mean()) if ob == ob else np.nan})
P = pd.DataFrame(pw); P.to_csv(D / "close_relative_power.csv", index=False); print(P.round(3).to_string(index=False))

print("\n4. interaction (PGLS, all bacteria): effect in carrier phyla (Bacillota + Pseudomonadota) vs Actinomycetota + Bacteroidota vs the ten other phyla")
C = ph.vcv_fast(tree, list(B.ncbi_assembly))
tgt = B.phylum.isin(OTHERS).astype(float).to_numpy(); oth = (~B.phylum.isin(CARRIERS + OTHERS)).astype(float).to_numpy()
cov = ["ph_opt"] + [c for c in ("src_thermobase", "src_acidb", "src_bacdive", "src_literature_auto") if B[c].nunique() > 1]
ints = []
for f, col in FEATS.items():
    x = z(B[col])
    X = np.column_stack([np.ones(len(B))] + [B[c].to_numpy(float) for c in cov] + [tgt, oth, x, x * tgt, x * oth])
    r = ph.pgls(B.ph_min.to_numpy(float), X, C)
    ints.append({"feature": f, "effect_carrier_phyla": r["beta"][-3], "se": r["se"][-3], "extra_in_Actino+Bacteroidota": r["beta"][-2], "se_a": r["se"][-2], "p_a": r["p"][-2],
                 "extra_in_other_phyla": r["beta"][-1], "se_o": r["se"][-1], "p_o": r["p"][-1], "implied_effect_Actino+Bacteroidota": r["beta"][-3] + r["beta"][-2]})
I = pd.DataFrame(ints); I.to_csv(D / "close_relative_interaction.csv", index=False); print(I.round(4).to_string(index=False))
