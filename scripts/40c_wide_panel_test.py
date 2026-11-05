"""Does the exposed-periplasmic result hold on a panel I did not choose? -> DATA_DIR/wide_panel_test.csv
SCALE WARNING, correction 15. The `pooled_beta_per_pt` column this writes is per STANDARD DEVIATION of that class's net-acidic share, not per percentage point:
ph.pgls_scan standardises the predictor it is given and nothing here divides the SD back out. The shares' SDs run about 2.5 to 6.5 points, so the per-point
values are three to six times smaller. Per SD is the right scale for comparing one position class against another, which is all this script's pooled table
does, so the figures are correct and only the name was wrong; scripts/74 carries both scales and the conversion.
Uses the systematically selected references from 40a (32 membrane proteins with at least 4 transmembrane helices and at least 3 periplasmic and 3 cytoplasmic topological domains,
a single Pfam family present in at least 80 % of sampled genomes, and an AlphaFold model passing the pLDDT filters) plus 5 conserved soluble controls. Only three of the 32
membrane families overlap the hand-picked panel of 36a (SecY PF00344, the NADH-quinone / hydrogenase family PF00361, and cytochrome bo3 PF00115), so dropping those three gives a
replication on references chosen by the selection rules alone - mostly transporters and permeases - which is the point of the exercise.
Reported for each structural class: random-effects pooling of the per-reference PGLS estimates, then a composite predictor tested across bacteria and within families by
permutation, in each case with genome-wide acidic excess controlled, and the same again after excluding the three shared families.
Usage: python3 scripts/40c_wide_panel_test.py [n_permutations]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
D, T = ROOT / os.environ.get("DATA_DIR", "results/broad2"), ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
SHARED = {"PF00344", "PF00361", "PF00115"}
MCL = ["peri_exp", "peri_bur", "cyto_exp", "cyto_bur", "TM"]
SCL = ["surface", "buried"]
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
SOL = [str(r.label) for r in refs.itertuples() if not r.membrane]
PF = {str(r.label): str(r.pfam) for r in refs.itertuples()}
NOVEL = [m for m in MEM if PF[m] not in SHARED]
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
df = (pd.read_csv(D / "phenotypes_all.csv").merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess", "gc", "genome_size_bp"]], on="ncbi_assembly")
      .merge(pd.read_csv(D / "structure_panel_wide.csv"), on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick"); tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips)].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
Call = ph.vcv_fast(tree, list(B.ncbi_assembly))
B["gw"] = B.acidic_excess * 100
B["gcz"] = (B.gc - B.gc.mean()) / B.gc.std(); B["sz"] = (B.genome_size_bp - B.genome_size_bp.mean()) / B.genome_size_bp.std()
SRCC = [f"src_{o}" for o in SRC if B[f"src_{o}"].nunique() > 1]
print(f"{len(MEM)} membrane references ({len(NOVEL)} not shared with the 36a panel), {len(SOL)} soluble; {len(B)} bacteria\n")


def excess(ref, cls):
    if f"{ref}|{cls}|n" not in B or (ref, cls) not in reflen:
        return None
    n = B[f"{ref}|{cls}|n"]
    return ((B[f"{ref}|{cls}|D"] + B[f"{ref}|{cls}|E"] - B[f"{ref}|{cls}|K"] - B[f"{ref}|{cls}|R"]) / n).where(n >= max(8, 0.5 * reflen[(ref, cls)])) * 100


def pgls(mask, x, cov):
    g = B[mask].reset_index(drop=True); Cm = Call[np.ix_(mask, mask)]
    use = [c for c in cov if g[c].nunique() > 1]
    X = np.column_stack([np.ones(len(g))] + [g[c] for c in use])
    b, se, p, lam = ph.pgls_scan(Cm, g.ph_min.to_numpy(float), X, np.asarray(x, float)[:, None])
    return b[0], se[0], p[0], len(g)


def pool(bs, ses):
    bs, ses = np.asarray(bs), np.asarray(ses)
    w = 1 / ses ** 2; fe = (w * bs).sum() / w.sum(); Q = (w * (bs - fe) ** 2).sum(); k = len(bs) - 1
    tau2 = max(0.0, (Q - k) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 0 else 0.0
    wr = 1 / (ses ** 2 + tau2); mu = (wr * bs).sum() / wr.sum(); sem = np.sqrt(1 / wr.sum())
    return mu, sem, 2 * stats.norm.sf(abs(mu / sem)), (max(0.0, (Q - k) / Q) if k > 0 and Q > 0 else 0.0)


rows = []
print("Per-reference PGLS, pooled by class (beta per percentage point, genome-wide acidic excess controlled)")
for grp, names in (("membrane, all 32", MEM), ("membrane, novel only", NOVEL), ("soluble", SOL)):
    for cls in (SCL if grp == "soluble" else MCL):
        bs, ses, n_ok = [], [], 0
        for r in names:
            x = excess(r, cls)
            if x is None:
                continue
            m = x.notna().to_numpy()
            if m.sum() < 200:
                continue
            b, se, p, nn = pgls(m, x[m].to_numpy(), ["gw"] + SRCC)
            bs.append(b); ses.append(se); n_ok += 1
        if n_ok < 2:
            continue
        mu, sem, pv, i2 = pool(bs, ses)
        # NAME IS A MISNOMER, kept for the files already written: `pooled_beta_per_pt` is per STANDARD DEVIATION of
        # that class's share, not per percentage point, because ph.pgls_scan standardises the predictor it tests and
        # nothing here divides the SD back out. Correction 15; scripts/74 gives both scales. Per SD is the right
        # scale for the class comparison this script makes, so the figures are fine and only the label was wrong.
        rows.append({"group": grp, "class": cls, "references": n_ok, "pooled_beta_per_pt": mu,
                     "scale": "per SD of the class share, NOT per percentage point",
                     "se": sem, "p": pv, "I2": i2, "positive": int(np.sum(np.array(bs) > 0))})
        print(f"  {grp:22s} {cls:9s} {n_ok:2d} refs: {mu:+.4f} +/- {sem:.4f} (p {pv:.1e}); positive in {int(np.sum(np.array(bs) > 0))}/{n_ok}; I2 {i2:.0%}")


def composite(names, cls):
    cols = [excess(r, cls) for r in names]
    cols = [(c - c.mean()) / c.std() for c in cols if c is not None]
    M = pd.concat(cols, axis=1)
    return M.mean(axis=1).where(M.notna().sum(axis=1) >= 3)


print("\nComposite predictor (beta per SD of the composite)")
rng = np.random.default_rng(20261025)
for grp, names, classes in (("membrane, all 32", MEM, MCL), ("membrane, novel only", NOVEL, MCL), ("soluble", SOL, SCL)):
    for cls in classes:
        pr = composite(names, cls)
        # pgls_scan standardises the tested predictor internally, so its beta is ALREADY per SD of that predictor.
        # (An earlier version multiplied by pr.std() here, which silently reported per UNIT of the composite.)
        for lab, cov, need in (("plain", SRCC, False), ("+ composition", ["gw"] + SRCC, False),
                               ("+ composition, GC, size", ["gw", "gcz", "sz"] + SRCC, False), ("+ optimum", ["ph_opt"] + SRCC, True)):
            m = pr.notna().to_numpy() & (B.ph_opt.notna().to_numpy() if need else True)
            if m.sum() < 200:
                continue
            b, se, p, nn = pgls(m, pr[m].to_numpy(), cov)
            rows.append({"group": grp, "class": cls, "model": f"PGLS {lab}", "n": nn, "beta_per_sd": b, "se": se, "p": p})
            if cls in ("peri_exp", "surface"):
                print(f"  {grp:22s} {cls:9s} {lab:24s} n={nn:4d} {b:+.4f} +/- {se:.4f} (p {p:.1e})")
        # within families, permutation
        d = B.assign(x=pr)[["family", "ph_min", "x", "gw"]].dropna()
        d = d[d.family.isin(d.family.value_counts()[lambda s: s >= 4].index)]
        if len(d) < 200:
            continue
        A = np.column_stack([np.ones(len(d)), d.gw])
        xr = d.x - A @ np.linalg.lstsq(A, d.x, rcond=None)[0]
        yr = d.ph_min - A @ np.linalg.lstsq(A, d.ph_min, rcond=None)[0]
        dd = d.assign(x=xr, y=yr)
        c = dd.groupby("family")[["y", "x"]].transform(lambda s: s - s.mean())
        sdc = c.x.std()
        obs = np.polyfit(c.x, c.y, 1)[0]
        gidx = [np.flatnonzero((dd.family == f).to_numpy()) for f in dd.family.unique()]
        xv, yv = c.x.to_numpy(), c.y.to_numpy(); null = np.empty(NPERM)
        for i in range(NPERM):
            xp = xv.copy()
            for ix in gidx:
                xp[ix] = rng.permutation(xp[ix])
            null[i] = np.polyfit(xp, yv, 1)[0]
        pv = (1 + np.sum(np.abs(null) >= abs(obs))) / (NPERM + 1)
        rows.append({"group": grp, "class": cls, "model": f"within families ({dd.family.nunique()} fam), composition removed", "n": len(dd), "beta_per_sd": obs * sdc, "se": null.std() * sdc, "p": pv})
        if cls in ("peri_exp", "surface"):
            print(f"  {grp:22s} {cls:9s} {'within families':24s} n={len(dd):4d} {obs * sdc:+.4f} per SD (permutation p {pv:.4f}, {dd.family.nunique()} families)")
pd.DataFrame(rows).to_csv(D / "wide_panel_test.csv", index=False)
print("\nwrote", D / "wide_panel_test.csv")
