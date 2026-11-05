"""The evidence-tier test, on the wide panel at last -> DATA_DIR/tier_test_wide_panel.csv
The synthesis has carried this caveat since the wide panel was built: "The tier comparison has not been repeated on the 32-protein panel, so the +0.225 should be read as the
hand-picked composite's value in that subset, not as the headline figure." It matters more than most outstanding checks, because the evidence tier is the test that killed every
earlier lead in this project. Histidine, Trk/Ktr and the Kdp ATPase all looked respectable in the pooled data and all failed among the species whose growth limits were read
out of the primary literature by hand. The exposed-periplasmic result passed that test on the hand-picked 7-reference panel, and the claim has rested on that ever since
without the rule-selected panel ever being put through it.
The tiers, from the phenotype table:
  origin        how the growth range was obtained - `literature` (hand-curated from species descriptions) is the good tier, `literature_auto` is the rule-based parser,
                and `thermobase`/`acidb`/`bacdive` are database rows
  ph_min_bound  whether the reported minimum is a measured value or merely the lowest pH tested
  confidence    the curation confidence recorded at the time
Each subset is fitted with genome-wide composition controlled, as the headline is. The source indicators are dropped inside a single-origin subset, where they would be
constant. Both the unweighted mean of the 19 high-coverage references and the factor score are tested, since scripts/53 and 56 showed the two differ on the coefficient, and
the hand-picked composite is carried alongside so the new numbers can be read against the +0.225 the claim has been resting on.
The hand-curated tier is small, so the question is not whether it reaches significance but whether the estimate stays positive and of comparable size. A tier estimate near
zero would be the serious outcome, and it is what happened to the earlier leads.
Usage: python3 scripts/58_tier_test_wide_panel.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
WIDE_MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
HAND_MEM = ["YidC", "SecY", "SecY_Bsub", "ATPsynC", "NuoL", "NuoN", "CyoB"]
rl_w = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
rl_o = pd.read_csv(ROOT / "results/structure_panel_classes.csv").groupby(["ref", "class"]).size()
wide = pd.read_csv(D / "structure_panel_wide.csv").drop_duplicates("ncbi_assembly")
orig = pd.read_csv(D / "structure_panel.csv").drop_duplicates("ncbi_assembly")
pheno = pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = pheno[(pheno.domain == "Bacteria") & pheno.ncbi_assembly.isin(tips)].drop_duplicates("ncbi_assembly").reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
Vcv = ph.vcv_fast(tree, list(B.ncbi_assembly))


def per_ref(tab, reflen, label):
    c = f"{label}|peri_exp|n"
    if c not in tab:
        return None
    m = tab.set_index("ncbi_assembly")
    n_ = m[c].reindex(B.ncbi_assembly)
    floor = max(8, 0.5 * reflen[(label, "peri_exp")]) if (label, "peri_exp") in reflen else 8
    num = (m[f"{label}|peri_exp|D"] + m[f"{label}|peri_exp|E"] - m[f"{label}|peri_exp|K"] - m[f"{label}|peri_exp|R"]).reindex(B.ncbi_assembly)
    return pd.Series((num / n_).where(n_ >= floor).to_numpy(), index=B.index)


def build(tab, reflen, labels, cover=0.0, weights=None):
    cols = {}
    for r in labels:
        x = per_ref(tab, reflen, r)
        if x is None or x.notna().mean() < cover:
            continue
        cols[r] = (x - x.mean()) / x.std()
    Z = pd.DataFrame(cols)
    if weights == "pc":
        cc = Z.notna().all(axis=1)
        w = np.linalg.eigh(np.corrcoef(Z[cc].to_numpy(float), rowvar=False))[1][:, -1]
        if w.sum() < 0:
            w = -w
        W = pd.DataFrame(np.tile(w, (len(Z), 1)), columns=Z.columns, index=Z.index).where(Z.notna())
        s = (Z * W).sum(axis=1, min_count=3) / W.abs().sum(axis=1)
    else:
        s = Z.mean(axis=1)
    return s.where(Z.notna().sum(axis=1) >= 3), len(cols)


PREDS = {}
PREDS["19 high-coverage, unweighted"], k19 = build(wide, rl_w, WIDE_MEM, cover=0.80)
PREDS["19 high-coverage, factor-weighted"], _ = build(wide, rl_w, WIDE_MEM, cover=0.80, weights="pc")
PREDS["all 32, unweighted"], k32 = build(wide, rl_w, WIDE_MEM)
PREDS["hand-picked 7"], k7 = build(orig, rl_o, HAND_MEM)
print(f"composites: 19 high-coverage ({k19} refs), all 32 ({k32}), hand-picked ({k7})\n")

SUBSETS = [
    ("all bacteria", pd.Series(True, index=B.index)),
    ("origin = literature (hand-curated)", B.origin.eq("literature")),
    ("origin = literature_auto (parser)", B.origin.eq("literature_auto")),
    ("origin = database rows only", B.origin.isin(["thermobase", "acidb", "bacdive"])),
    ("bound = measured", B.ph_min_bound.eq("measured")),
    ("bound != measured", ~B.ph_min_bound.eq("measured")),
    ("confidence = high or medium", B.confidence.isin(["high", "medium"])),
    ("confidence = low", B.confidence.eq("low")),
]
rows = []
for label, pr in PREDS.items():
    print(f"{label}")
    for sname, mask in SUBSETS:
        m = (pr.notna() & mask & B.ph_min.notna()).to_numpy()
        if m.sum() < 40:
            print(f"   {sname:36s} n={int(m.sum()):4d}  too few, skipped")
            continue
        g = B[m]
        use = [c for c in ["gw"] + [f"src_{o}" for o in SRC] if g[c].nunique() > 1]
        X = np.column_stack([np.ones(int(m.sum()))] + [g[c].to_numpy(float) for c in use])
        b, se, p, lam = ph.pgls_scan(Vcv[np.ix_(m, m)], g.ph_min.to_numpy(float), X, pr[m].to_numpy(float)[:, None])
        rows.append({"composite": label, "subset": sname, "n": int(m.sum()), "beta_per_sd": float(b[0]), "se": float(se[0]), "p": float(p[0]), "lambda": float(lam)})
        print(f"   {sname:36s} n={int(m.sum()):4d}  {b[0]:+.4f} +/- {se[0]:.4f} (p {p[0]:.4f})")
    print()

R = pd.DataFrame(rows)
R.to_csv(D / "tier_test_wide_panel.csv", index=False)
print("wrote", D / "tier_test_wide_panel.csv")
hand = R[R.subset.str.startswith("origin = literature (")]
allb = R[R.subset == "all bacteria"].set_index("composite")
print("\nthe comparison the synthesis has been waiting for - hand-curated tier against the whole set, same model:")
for x in hand.itertuples():
    a = allb.loc[x.composite]
    z = (x.beta_per_sd - a.beta_per_sd) / np.sqrt(x.se ** 2 + a.se ** 2)
    print(f"   {x.composite:36s} hand-curated {x.beta_per_sd:+.4f} +/- {x.se:.4f} (n={x.n})  vs all {a.beta_per_sd:+.4f} (n={int(a.n)})  "
          f"difference {x.beta_per_sd - a.beta_per_sd:+.4f}, z {z:+.2f}")
print("\nFor reference, the claim has rested on +0.225 from the hand-picked composite in the 87-species hand-curated tier.")
