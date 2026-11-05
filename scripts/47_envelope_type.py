"""Is the effect larger where the exposed surface is less buffered? -> DATA_DIR/envelope_type.csv
A directional prediction, stated before the test was run and motivated by the lineage breakdown in 44, where the estimate was larger in Bacillota (+0.303) and Actinomycetota
(+0.211) than in Pseudomonadota (+0.087) and Bacteroidota (+0.144). Those pairs differ in envelope architecture: Bacillota and Actinomycetota are monoderms, whose
extracytoplasmic surface sits behind only a cell wall, while Pseudomonadota and Bacteroidota are diderms, whose periplasm is enclosed by an outer membrane that buffers it
against the external medium. If the association arises because external protons act on exposed extracytoplasmic surfaces, the effect should be LARGER in monoderms. If instead
it is a general compositional correlate of acid tolerance, envelope type should not matter.
Classification is by phylum and is deliberately coarse; the known exceptions (diderm Negativicutes inside Bacillota, the atypical envelopes of Deinococcota) are not separated,
which if anything blunts the contrast rather than creating it. Only the four phyla with enough species are used, so the comparison is 822 monoderms against 1,671 diderms.
Usage: python3 scripts/47_envelope_type.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import os, sys, warnings
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np, pandas as pd
from Bio import Phylo
from scipy import stats
from acidophile import phylo as ph
warnings.filterwarnings("ignore")
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
MONODERM = {"Bacillota", "Actinomycetota"}
DIDERM = {"Pseudomonadota", "Bacteroidota"}
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
SOL = [str(r.label) for r in refs.itertuples() if not r.membrane]
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pd.read_csv(D / "structure_panel_wide.csv"), on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
for o in ("thermobase", "acidb", "bacdive", "literature_auto"):
    B[f"src_{o}"] = (B.origin == o).astype(float)
Call = ph.vcv_fast(tree, list(B.ncbi_assembly))
z = lambda s: (s - s.mean()) / s.std()
B["gw"] = B.acidic_excess * 100
SRCC = [f"src_{o}" for o in ("thermobase", "acidb", "bacdive", "literature_auto")]
def composite(names, cls):
    cols = []
    for r in names:
        if f"{r}|{cls}|n" not in B or (r, cls) not in reflen: continue
        n = B[f"{r}|{cls}|n"]
        cols.append(z(((B[f"{r}|{cls}|D"] + B[f"{r}|{cls}|E"] - B[f"{r}|{cls}|K"] - B[f"{r}|{cls}|R"]) / n).where(n >= max(8, 0.5 * reflen[(r, cls)]))))
    M = pd.concat(cols, axis=1)
    return M.mean(axis=1).where(M.notna().sum(axis=1) >= 3)
def fit(mask, pr, extra=()):
    g = B[mask].reset_index(drop=True); Cm = Call[np.ix_(mask, mask)]
    use = [c for c in ["gw"] + list(extra) + SRCC if c in g and g[c].nunique() > 1]
    X = np.column_stack([np.ones(len(g))] + [g[c] for c in use])
    b, se, p, lam = ph.pgls_scan(Cm, g.ph_min.to_numpy(float), X, pr[mask].to_numpy(float)[:, None])
    return b[0], se[0], p[0], len(g)
PRED = {"exposed periplasmic": composite(MEM, "peri_exp"), "exposed cytoplasmic": composite(MEM, "cyto_exp"), "soluble surface": composite(SOL, "surface")}
B["envelope"] = np.where(B.phylum.isin(MONODERM), "monoderm", np.where(B.phylum.isin(DIDERM), "diderm", "other"))
print(f"monoderms {int((B.envelope == 'monoderm').sum())}, diderms {int((B.envelope == 'diderm').sum())}, other phyla {int((B.envelope == 'other').sum())} (excluded)\n")
rows = []
for name, pr in PRED.items():
    est = {}
    for env in ("monoderm", "diderm"):
        m = (pr.notna() & (B.envelope == env)).to_numpy()
        b, se, p, n = fit(m, pr)
        est[env] = (b, se)
        rows.append({"predictor": name, "group": env, "n": n, "beta_per_sd": b, "se": se, "p": p})
        print(f"  {name:20s} {env:9s} n={n:4d}  {b:+.4f} +/- {se:.4f} (p {p:.4f})")
    d = est["monoderm"][0] - est["diderm"][0]; sd = float(np.hypot(est["monoderm"][1], est["diderm"][1]))
    pv = stats.norm.sf(d / sd)   # one-sided: the prediction is monoderm > diderm
    rows.append({"predictor": name, "group": "monoderm - diderm", "n": np.nan, "beta_per_sd": d, "se": sd, "p": pv})
    print(f"  {name:20s} difference {d:+.4f} +/- {sd:.4f}, one-sided p {pv:.4f} (prediction: monoderm larger)\n")
pd.DataFrame(rows).to_csv(D / "envelope_type.csv", index=False)
