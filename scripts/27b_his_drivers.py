"""What drives the histidine-fraction association with pH_min? (bacteria; archaea shown for contrast)
 A. per-amino-acid scan (raw and centred-log-ratio composition): is histidine special?
 B. covariate ladder: GC, genome size, thermal index, composition principal components (histidine excluded), the full composition, proteome properties, growth temperature
 C. optimum pathway: histidine vs pH_opt, and pH_min with / without the optimum
 D. where the histidines sit: clustered (metal-binding spacings) vs isolated, membrane vs soluble
Usage: python3 scripts/27b_his_drivers.py   (env DATA_DIR=results/broad, PHYLO_DIR=results/phylo_broad)"""
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
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad")
AA = "ACDEFGHIKLMNPQRSTVWY"
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv"), on="ncbi_assembly")
      .merge(pd.read_csv(D / "hmm_counts.csv").drop(columns=["n_genes"]), left_on="ncbi_assembly", right_on="genome_id"))
hb = D / "his_breakdown.csv"
if hb.exists():
    df = df.merge(pd.read_csv(hb), on="ncbi_assembly", how="left")
for o in ("thermobase", "acidb", "bacdive", "literature_auto"):
    df[f"src_{o}"] = (df.origin == o).astype(float)
df["log_genome"] = np.log10(df.genome_size_bp)
F = df[[f"frac_{a}" for a in AA]].to_numpy(float) + 1e-6
CLR = np.log(F) - np.log(F).mean(axis=1, keepdims=True)
for i, a in enumerate(AA):
    df[f"clr_{a}"] = CLR[:, i]
z = lambda v: (np.asarray(v, float) - np.nanmean(v)) / np.nanstd(v)
OUT = []


def setup(dom):
    tree = Phylo.read(T / f"{dom}.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
    g = df[(df.domain == dom) & df.ph_opt.notna() & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
    return g, ph.vcv_fast(tree, list(g.ncbi_assembly))


def base_cols(g, with_opt=True):
    return (["ph_opt"] if with_opt else []) + [c for c in ("src_thermobase", "src_acidb", "src_bacdive", "src_literature_auto") if g[c].nunique() > 1]


def scan(g, C, feat, extra=(), y="ph_min", with_opt=True, mask=None):
    m = np.ones(len(g), bool) if mask is None else np.asarray(mask)
    h = g[m].reset_index(drop=True); Cm = C[np.ix_(m, m)]
    cols = base_cols(h, with_opt) + list(extra)
    X = np.column_stack([np.ones(len(h))] + [h[c].to_numpy(float) for c in cols])
    if not isinstance(feat, (list, tuple)):
        feat = [feat]
    b, se, p, lam = ph.pgls_scan(Cm, h[y].to_numpy(float), X, h[list(feat)].to_numpy(float))
    return b, se, p, lam, len(h)


pd.set_option("display.width", 220)
for dom in ("Bacteria", "Archaea"):
    g, C = setup(dom)
    print(f"\n########## {dom}: n={len(g)} ##########")
    # A. per-amino-acid scans
    rawb, rawse, rawp, _, _ = scan(g, C, [f"frac_{a}" for a in AA])
    clrb, clrse, clrp, _, _ = scan(g, C, [f"clr_{a}" for a in AA])
    gcb, _, gcp, _, _ = scan(g, C, [f"frac_{a}" for a in AA], extra=("gc", "log_genome"))
    A = pd.DataFrame({"aa": list(AA), "raw_beta": rawb, "raw_p": rawp, "clr_beta": clrb, "clr_p": clrp, "raw_beta_adj_GC_size": gcb, "raw_p_adj_GC_size": gcp,
                      "corr_with_GC": [np.corrcoef(g[f"frac_{a}"], g.gc)[0, 1] for a in AA], "corr_with_His": [np.corrcoef(g[f"frac_{a}"], g.frac_H)[0, 1] if "frac_H" in g else np.corrcoef(g[f"frac_{a}"], g.frac_His)[0, 1] for a in AA]})
    A.insert(0, "domain", dom); OUT.append(A.assign(table="A_per_amino_acid"))
    print("A. per-amino-acid effect on pH_min (pH_opt + sources as covariates; beta per SD):")
    print(A.sort_values("raw_p").head(8).round(4).to_string(index=False))
    H = A[A.aa == "H"].iloc[0]; print(f"   histidine rank by raw p among 20: {int((A.raw_p < H.raw_p).sum()) + 1}; by CLR p: {int((A.clr_p < H.clr_p).sum()) + 1}")
    print(f"   corr(His, GC) = {np.corrcoef(g.frac_His, g.gc)[0, 1]:+.2f}; corr(His, log genome size) = {np.corrcoef(g.frac_His, g.log_genome)[0, 1]:+.2f}; corr(His, IVYWREL) = {np.corrcoef(g.frac_His, g.ivywrel)[0, 1]:+.2f}")
    # B. covariate ladder for histidine
    comp19 = g[[f"clr_{a}" for a in AA if a != "H"]].to_numpy(float); comp19 = comp19 - comp19.mean(axis=0)
    U, S, Vt = np.linalg.svd(comp19, full_matrices=False); pcs = U * S
    for k in range(10):
        g[f"pc{k + 1}"] = pcs[:, k]
    ladder = [("base (pH_opt + sources)", ()), ("+ GC", ("gc",)), ("+ GC + genome size", ("gc", "log_genome")), ("+ GC + genome size + IVYWREL", ("gc", "log_genome", "ivywrel")),
              ("+ composition PC1-5 (His excluded)", tuple(f"pc{k}" for k in range(1, 6))), ("+ composition PC1-10", tuple(f"pc{k}" for k in range(1, 11))),
              ("+ all 19 other amino acids (CLR)", tuple(f"clr_{a}" for a in AA if a != "H")),
              ("+ mean pI, GRAVY, membrane fraction", ("mean_pI", "mean_gravy", "frac_membrane_proxy")), ("+ GC + size + PC1-10 + proteome properties", ("gc", "log_genome", "mean_pI", "mean_gravy", "frac_membrane_proxy") + tuple(f"pc{k}" for k in range(1, 11)))]
    rows = []
    for lab, ex in ladder:
        b, se, p, lam, n = scan(g, C, "frac_His", extra=ex); rows.append({"domain": dom, "adjustment": lab, "n": n, "beta_per_sd": b[0], "se": se[0], "p": p[0], "lambda": lam})
    tk = g.temp_opt_c.notna().to_numpy()
    if tk.sum() > 100:
        for lab, ex in (("[temperature-known subset] base", ()), ("[temperature-known subset] + growth temperature", ("temp_opt_c",))):
            b, se, p, lam, n = scan(g, C, "frac_His", extra=ex, mask=tk); rows.append({"domain": dom, "adjustment": lab, "n": n, "beta_per_sd": b[0], "se": se[0], "p": p[0], "lambda": lam})
    L = pd.DataFrame(rows); OUT.append(L.assign(table="B_covariate_ladder"))
    print("\nB. histidine fraction effect on pH_min under increasing adjustment:"); print(L.round(4).to_string(index=False))
    # C. optimum pathway
    rows = []
    for lab, y, wo in (("His -> pH_min | pH_opt (direct effect)", "ph_min", True), ("His -> pH_min (no optimum covariate; total effect)", "ph_min", False), ("His -> pH_opt", "ph_opt", False)):
        b, se, p, lam, n = scan(g, C, "frac_His", y=y, with_opt=wo); rows.append({"domain": dom, "model": lab, "n": n, "beta_per_sd": b[0], "se": se[0], "p": p[0], "lambda": lam})
    P = pd.DataFrame(rows); OUT.append(P.assign(table="C_optimum_pathway")); print("\nC. optimum pathway:"); print(P.round(4).to_string(index=False))
    tot, direct = P.beta_per_sd.iloc[1], P.beta_per_sd.iloc[0]; print(f"   share of the total His effect on pH_min that runs through the optimum: {(tot - direct) / tot:.0%}")
    for lab, mk in (("optimum <= 5", g.ph_opt <= 5), ("optimum > 5", g.ph_opt > 5)):
        if mk.sum() > 40:
            b, se, p, lam, n = scan(g, C, "frac_His", mask=mk.to_numpy()); print(f"   His -> pH_min | pH_opt within {lab} (n={n}): {b[0]:+.3f} ± {se[0]:.3f} (p {p[0]:.3f})")
    # D. where the histidines sit
    if "his_clustered_frac" in g and g.his_clustered_frac.notna().any():
        cols = ["his_clustered_frac", "his_isolated_frac", "his_frac_soluble", "his_frac_membrane", "frac_proteins_his_rich"]
        rows = []
        for c in cols:
            ok = g[c].notna().to_numpy()
            for lab, ex in (("base", ()), ("+ GC + genome size", ("gc", "log_genome"))):
                b, se, p, lam, n = scan(g, C, c, extra=ex, mask=ok); rows.append({"domain": dom, "component": c, "adjustment": lab, "n": n, "beta_per_sd": b[0], "se": se[0], "p": p[0]})
        Dd = pd.DataFrame(rows); OUT.append(Dd.assign(table="D_his_location")); print("\nD. components of the histidine fraction (each alone):"); print(Dd.round(4).to_string(index=False))
        h = g.dropna(subset=["his_clustered_frac", "his_isolated_frac"]).reset_index(drop=True)
        X = np.column_stack([np.ones(len(h))] + [h[c].to_numpy(float) for c in base_cols(h)] + [z(h.his_clustered_frac), z(h.his_isolated_frac)])
        r = ph.pgls(h.ph_min.to_numpy(float), X, ph.vcv_fast(Phylo.read(T / f"{dom}.nwk", "newick"), list(h.ncbi_assembly)))
        print(f"   joint model: clustered His {r['beta'][-2]:+.3f} ± {r['se'][-2]:.3f} (p {r['p'][-2]:.3f}); isolated His {r['beta'][-1]:+.3f} ± {r['se'][-1]:.3f} (p {r['p'][-1]:.3f})  [correlation between them: {np.corrcoef(h.his_clustered_frac, h.his_isolated_frac)[0, 1]:+.2f}]")
pd.concat(OUT, ignore_index=True).to_csv(D / "his_drivers.csv", index=False)

# E. (bacteria) is histidine a proxy for gene-system presence, and is it the same signal as Kdp?
g, C = setup("Bacteria")
hm = [c for c in g.columns if c.startswith("hm_rate_") and g[c].std() > 0]
rows = []
for lab, ex in (("base", ()), ("+ all Pfam gene-system rates", tuple(hm)), ("+ Kdp rate", ("hm_rate_kdp",)), ("+ Kdp + Trk/Ktr", ("hm_rate_kdp", "hm_rate_trk_ktr")),
                ("+ all gene-system rates + GC + genome size", tuple(hm) + ("gc", "log_genome"))):
    b, se, p, lam, n = scan(g, C, "frac_His", extra=ex); rows.append({"adjustment_for_His": lab, "n": n, "His_beta": b[0], "His_se": se[0], "His_p": p[0]})
b, se, p, lam, n = scan(g, C, "hm_rate_kdp", extra=("frac_His",)); rows.append({"adjustment_for_His": "[Kdp effect] + His", "n": n, "His_beta": b[0], "His_se": se[0], "His_p": p[0]})
E = pd.DataFrame(rows); E.to_csv(D / "his_vs_gene_systems.csv", index=False)
print("\nE. histidine vs gene-system presence and Kdp (bacteria):"); print(E.round(4).to_string(index=False))
print(f"   corr(His, Kdp rate) = {np.corrcoef(g.frac_His, g.hm_rate_kdp)[0, 1]:+.2f}; Spearman = {pd.Series(g.frac_His).corr(pd.Series(g.hm_rate_kdp), method='spearman'):+.2f}")
