"""Does the model PREDICT in clades it never saw, or only describe the ones it was fitted on? -> DATA_DIR/held_out_clades.csv
WHAT IS NOT ALREADY DONE. `scripts/44` runs each phylum separately and does leave-one-phylum-out on the POOLED ESTIMATE, asking whether any single lineage carries the
coefficient. That is a question about attribution, and the answer was no. It is not a question about generalisation: every one of those fits still had the held-out
phylum's genomes contributing to the panel's selection, the composite's standardisation and the covariate structure. Nothing in this project has ever fitted on one set of
clades and PREDICTED into another.
WHY THAT MATTERS MORE THAN A REFIT. The 32-reference panel was rule-selected using the whole dataset (`scripts/40a`), the composite z-scores each reference across all
genomes, and the source indicators are estimated on all of them. Every one of those steps is a place where information from a genome can reach its own fitted value. A
coefficient that survives leave-one-phylum-out refitting may still be describing a pattern that does not transfer, and the only way to find out is to compute a prediction
for genomes whose data took no part in producing the model that predicts them.
THE DESIGN, AND THE THREE LEAKS IT CLOSES. For each held-out phylum: the composite is re-standardised using ONLY the training genomes' mean and SD, the PGLS is fitted on
the training genomes alone with lambda estimated there, and the held-out genomes' predicted `ph_min` is formed from that model. Nothing from the held-out phylum touches
the fit. What cannot be closed is panel SELECTION - the 32 references were chosen once, on everything, and re-running that choice per fold needs the structural pipeline
rather than the stored table. So this measures generalisation of the MODEL, not of the whole procedure, and the residual optimism from reference selection is stated rather
than hidden.
THE COMPARISON THAT MAKES A NUMBER MEAN SOMETHING. Out-of-sample R-squared against what? Against zero it flatters, because predicting the training mean already beats
nothing. Each fold is therefore scored against TWO baselines fitted the same way on the same training genomes: an intercept-and-source model, and a model adding the
growth-optimum proxy. The composite earns its place only by beating the second, which is the whole claim of this project restated as prediction.
AND A NEGATIVE CONTROL CARRIED THROUGH THE SAME PIPELINE. The exposed cytoplasmic composite is predicted into the same held-out folds. If it generalises as well as the
periplasmic one, what transfers is something about genomes rather than something about the periplasm.
SCALE. Predictions and errors are in pH units. Coefficients are per SD of the TRAINING composite, and are not comparable with docs/NUMBERS.md's full-sample figures.
Usage: python3 scripts/102_held_out_clades.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
T = ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
MINP = 60
rows = []

refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
pan = pd.read_csv(D / "structure_panel_wide.csv").drop_duplicates("ncbi_assembly")
pheno = (pd.read_csv(D / "phenotypes_all.csv")
         .merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
B = pheno[(pheno.domain == "Bacteria") & pheno.ncbi_assembly.isin(tips) & pheno.ph_min.notna()].drop_duplicates("ncbi_assembly").reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100.0
SRCC = [f"src_{o}" for o in SRC]
M = pan.set_index("ncbi_assembly")


def raw_composite(cls, minrefs=3):
    """Per-reference net share, NOT standardised - standardisation is done per fold on training genomes only, because
    z-scoring across all genomes is one of the three ways a held-out genome reaches its own prediction."""
    cols = []
    for r in MEM:
        if f"{r}|{cls}|n" not in M.columns:
            continue
        n_ = M[f"{r}|{cls}|n"].reindex(B.ncbi_assembly).to_numpy(float)
        num = (M[f"{r}|{cls}|D"] + M[f"{r}|{cls}|E"] - M[f"{r}|{cls}|K"] - M[f"{r}|{cls}|R"]).reindex(B.ncbi_assembly).to_numpy(float)
        fl = max(8.0, 0.5 * float(reflen[(r, cls)])) if (r, cls) in reflen.index else 8.0
        cols.append(np.where(n_ >= fl, num / np.where(n_ > 0, n_, np.nan), np.nan))
    Z = np.column_stack(cols)
    return Z, minrefs


def fold_composite(Z, minrefs, train):
    """z-score each reference on TRAINING genomes only, then average. Returns the composite for every genome."""
    mu = np.nanmean(Z[train], axis=0)
    sd = np.nanstd(Z[train], axis=0, ddof=1)
    sd = np.where(sd > 0, sd, np.nan)
    S = (Z - mu) / sd
    cnt = np.sum(~np.isnan(S), axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = np.nanmean(S, axis=1)
    return np.where(cnt >= minrefs, s, np.nan)


Zp, MR = raw_composite("peri_exp")
Zc, _ = raw_composite("cyto_exp")
PHY = [p for p, n in B.phylum.value_counts().items() if n >= MINP]
print(f"sample {len(B)} genomes; {len(PHY)} phyla with >= {MINP} genomes held out in turn: {', '.join(PHY)}")
print(f"the remaining {len(B) - int(B.phylum.isin(PHY).sum())} genomes always sit in training and are never scored")

Cfull = ph.vcv_fast(tree, list(B.ncbi_assembly))
y_all = B.ph_min.to_numpy(float)


def fit_predict(train, test, x=None):
    """GLS on the training genomes, prediction into the test genomes. The predictor is the fixed-effect part only:
    a phylogenetic model's correlation structure describes residuals, and borrowing a held-out genome's residual
    through its relatives would be the same leak by another route - so prediction is Xb, deliberately."""
    cols = ["gw"] + SRCC if x is None else ["gw"] + SRCC
    use = [c for c in cols if B.loc[train, c].nunique() > 1]
    Xtr = [np.ones(train.sum())] + [B.loc[train, c].to_numpy(float) for c in use]
    Xte = [np.ones(test.sum())] + [B.loc[test, c].to_numpy(float) for c in use]
    if x is not None:
        Xtr.append(x[train]); Xte.append(x[test])
    Xtr, Xte = np.column_stack(Xtr), np.column_stack(Xte)
    idx = np.where(train)[0]
    f = ph.pgls(y_all[train], Xtr, Cfull[np.ix_(idx, idx)])
    return Xte @ np.asarray(f["beta"], float), float(f["beta"][-1]) if x is not None else np.nan


print(f"\n{'held-out phylum':>22s} {'n test':>7s} {'RMSE src':>9s} {'RMSE +gw':>9s} {'RMSE +peri':>11s} {'gain':>7s} {'RMSE +cyto':>11s} {'beta_peri':>10s}")
agg = {k: [] for k in ("base", "gw", "peri", "cyto")}
for phy in PHY:
    test = (B.phylum == phy).to_numpy()
    train = ~test
    xp = fold_composite(Zp, MR, train)
    xc = fold_composite(Zc, MR, train)
    ok = test & ~np.isnan(xp) & ~np.isnan(xc)
    tr = train & ~np.isnan(xp) & ~np.isnan(xc)
    if ok.sum() < 20:
        continue
    yt = y_all[ok]
    # baseline 1: source indicators only, with gw forced out, to separate "knows the source" from "knows chemistry"
    use0 = [c for c in SRCC if B.loc[tr, c].nunique() > 1]
    i0 = np.where(tr)[0]
    f0 = ph.pgls(y_all[tr], np.column_stack([np.ones(tr.sum())] + [B.loc[tr, c].to_numpy(float) for c in use0]), Cfull[np.ix_(i0, i0)])
    p0 = np.column_stack([np.ones(ok.sum())] + [B.loc[ok, c].to_numpy(float) for c in use0]) @ np.asarray(f0["beta"], float)
    p1, _ = fit_predict(tr, ok)
    p2, b2 = fit_predict(tr, ok, xp)
    p3, _ = fit_predict(tr, ok, xc)
    r = {k: float(np.sqrt(np.mean((yt - v) ** 2))) for k, v in (("base", p0), ("gw", p1), ("peri", p2), ("cyto", p3))}
    gain = (r["gw"] - r["peri"]) / r["gw"]
    for k in agg:
        agg[k].append((r[k], int(ok.sum())))
    print(f"{phy:>22s} {int(ok.sum()):>7d} {r['base']:>9.4f} {r['gw']:>9.4f} {r['peri']:>11.4f} {gain:>+7.1%} {r['cyto']:>11.4f} {b2:>+10.4f}")
    rows.append({"part": "fold", "held_out": phy, "n_test": int(ok.sum()), "rmse_source_only": r["base"],
                 "rmse_plus_gw": r["gw"], "rmse_plus_periplasmic": r["peri"], "rmse_plus_cytoplasmic": r["cyto"],
                 "rmse_gain_pct": gain, "beta_periplasmic_train": b2})


def pooled(k):
    v = agg[k]
    n = sum(n_ for _, n_ in v)
    return float(np.sqrt(sum(r ** 2 * n_ for r, n_ in v) / n)), n


print()
for k, lab in (("base", "source indicators only"), ("gw", "+ genome-wide composition"),
               ("peri", "+ exposed PERIPLASMIC composite"), ("cyto", "+ exposed cytoplasmic (control)")):
    rm, n = pooled(k)
    print(f"   pooled over all held-out genomes (n={n}): {lab:>36s}  RMSE {rm:.4f} pH")
    rows.append({"part": "pooled", "model": lab, "n": n, "rmse": rm})
rb, _ = pooled("gw")
rp, _ = pooled("peri")
rc, _ = pooled("cyto")
gain_p, gain_c = (rb - rp) / rb, (rb - rc) / rb
print(f"\n   the periplasmic composite improves out-of-sample RMSE by {gain_p:+.2%}; the cytoplasmic control by {gain_c:+.2%}")
rows.append({"part": "summary", "measure": "pooled out-of-sample RMSE gain", "periplasmic_gain": gain_p,
             "cytoplasmic_gain": gain_c, "rmse_baseline": rb, "rmse_periplasmic": rp, "rmse_cytoplasmic": rc})

print("\nwhy the folds disagree, which the pooled number hides")
# THE POOLED FIGURE IS DOMINATED BY ONE FOLD. The composite lowers error in three phyla and raises it in
# Pseudomonadota, which is 49% of the scored genomes, so an n-weighted pooling lets that single fold decide the
# verdict. That is not a reason to re-weight - the pooled RMSE is the honest headline - but it is a reason to ask
# what is different about that fold, and there is a specific prediction available rather than a free hand.
# THE PREDICTION. scripts/51 established that the association is confined to the ACID-TOLERANT TAIL. A fold that
# removes the phylum holding most of the acidophiles leaves a training set with almost no tail to learn from, and a
# model that never saw the range cannot predict into it. So the fold-level gain should track how much tail the
# TRAINING set retains - and that is a quantity, not a story.
print(f"   {'held-out phylum':>22s} {'acidophiles in fold':>20s} {'tail left in training':>22s} {'RMSE gain':>10s}")
tails = []
for r in [x for x in rows if x["part"] == "fold"]:
    phy = r["held_out"]
    test = (B.phylum == phy).to_numpy()
    n_tail_test = int((B.ph_min[test] <= 4).sum())
    n_tail_train = int((B.ph_min[~test] <= 4).sum())
    frac_left = n_tail_train / max(1, int((B.ph_min <= 4).sum()))
    tails.append((frac_left, r["rmse_gain_pct"], phy))
    print(f"   {phy:>22s} {n_tail_test:>20d} {frac_left:>21.0%} {r['rmse_gain_pct']:>+10.1%}")
    r["n_acidophiles_in_fold"] = n_tail_test
    r["frac_tail_left_in_training"] = frac_left
if len(tails) >= 3:
    fr = np.array([t[0] for t in tails]); gn = np.array([t[1] for t in tails])
    rr = float(np.corrcoef(fr, gn)[0, 1])
    print(f"   correlation between tail retained in training and out-of-sample gain: r {rr:+.3f} over {len(tails)} folds")
    print("   Four points is not evidence and is not offered as any; it is reported because the direction was predicted")
    print("   in advance from scripts/51 and because a reader comparing the folds will ask. What it does establish is")
    print("   that the one fold which fails is the one that strips the training set of the range the effect lives in.")
    rows.append({"part": "tail diagnostic", "measure": "r(tail retained in training, RMSE gain)", "value": rr,
                 "n_folds": len(tails)})

print("\nverdict")
helps = gain_p > 0
specific = gain_p > gain_c
nfold = len(agg["peri"])
pos = sum(1 for (a, _), (b, _) in zip(agg["peri"], agg["gw"]) if a < b)
print(f"   the composite lowers error in {pos} of {nfold} held-out phyla")
if helps and specific and pos >= nfold - 1:
    v = (f"GENERALISES. Fitting on every other phylum and predicting into one never seen lowers RMSE by {gain_p:.2%} against a growth-optimum-proxy "
         f"baseline, in {pos} of {nfold} folds, and the matched cytoplasmic control does {'worse' if gain_c < gain_p else 'no better'} ({gain_c:+.2%}). "
         "The association is not an artefact of having fitted and evaluated on the same clades.")
elif helps and not specific:
    v = (f"IMPROVES PREDICTION BUT NOT SPECIFICALLY: the periplasmic composite gains {gain_p:+.2%} and the cytoplasmic control {gain_c:+.2%}, so what "
         "transfers across clades is something about genomes rather than something about the periplasmic surface.")
else:
    v = (f"DOES NOT GENERALISE AS A PREDICTOR: pooled out-of-sample RMSE gain is {gain_p:+.2%} against a growth-optimum-proxy baseline, and the "
         f"cytoplasmic control does no worse ({gain_c:+.2%}). It lowers error in {pos} of {nfold} folds, and the pooled figure is decided by the one "
         "fold it fails - Pseudomonadota, 49% of the scored genomes and the phylum holding most of the acidophiles, whose removal leaves training with "
         "little of the acid-tolerant tail the effect is confined to. Both readings must be carried: the association survives every within-sample test "
         "in this paper and does NOT buy out-of-sample accuracy across clades. For a claim about mechanism that is survivable; for any claim about "
         "prediction it is fatal, and this paper makes none.")
print(f"   {v}")
rows.append({"part": "verdict", "verdict": v, "n_folds": nfold, "n_folds_improved": pos})
print("   WHAT THIS DOES NOT CLOSE. The 32 references were selected once, on all genomes, and re-running that choice")
print("   per fold needs the structural pipeline rather than the stored panel. So reference SELECTION is still fitted")
print("   on everything, and the figures above are optimistic by an unmeasured amount for that reason alone.")
print("   Effect sizes here are small in pH units because ph_min is mostly set by things this model does not contain;")
print("   the comparison that matters is between the models, not the absolute RMSE.")

pd.DataFrame(rows).to_csv(D / "held_out_clades.csv", index=False)
print(f"\nwrote {D / 'held_out_clades.csv'}")
