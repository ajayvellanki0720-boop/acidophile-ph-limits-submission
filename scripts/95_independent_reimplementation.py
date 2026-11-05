"""Do the headline results survive an implementation that shares no code with the one that produced them? -> DATA_DIR/independent_reimplementation.csv
WHY THIS IS WORTH A SCRIPT. Every number in this project comes from one 300-line module, `src/acidophile/phylo.py`, written once and used by 94 scripts. Twenty corrections
have been found in the ANALYSES; none has been found in the ESTIMATOR, and the reason may be that nothing has ever looked. A bug in `pgls` would not announce itself: it
would produce a plausible coefficient with a plausible standard error, and every robustness check in this project would inherit it, because every robustness check calls it.
The existing tests do not close this. `tests/` checks internal identities - that `pgls_scan` returns per-SD figures, that the two VCV builders agree - and two implementations
that agree can be wrong together whenever they share a primitive. `vcv_from_tree` and `vcv_fast` both get their branch lengths from Bio.Phylo's parse and both define
covariance through the same MRCA notion, so their agreement tests the traversal and not the definition.
WHAT IS INDEPENDENT HERE, stated precisely so the claim is checkable rather than rhetorical:
  (a) THE TREE. A recursive-descent Newick parser written here, not Bio.Phylo. If the published results rest on a mis-parsed branch length, both existing VCV functions
      would carry it identically and no test in the project could see it.
  (b) THE COVARIANCE. Built from shared root-to-tip path length via a longest-common-prefix comparison of ancestor paths, not by walking a tree and grouping sibling
      subtrees, and not by calling any distance function. Same mathematical object, no shared route to it.
  (c) THE ESTIMATOR. Lambda by dense grid search, not bounded Brent. Whitening by symmetric eigendecomposition, not Cholesky. Coefficients by explicit normal equations
      with a matrix inverse, not least squares on a whitened design. The log-determinant from slogdet, not from the sum of a Cholesky factor's log diagonal. Every one of
      these is a route a Cholesky bug could not reach.
  (d) THE PREDICTOR. The composite rebuilt from the stored panel's raw per-reference columns, including its position floor and its minimum-reference rule, without calling
      the helper that builds it.
WHAT IS **NOT** REPRODUCED, which matters more than what is. This reproduces the STATISTICAL MODEL and the PREDICTOR'S ARITHMETIC from stored intermediates. It does not
rebuild the structural panel, the HMM detection, the exposure calls or the phenotype curation - those need the structures and the literature pass, and an error in them
would survive this check untouched. The honest scope is: if the headline is wrong, this says whether the fitting is why. It cannot say the pipeline is right.
AND IT IS NOT A TEST THAT CAN ONLY PASS. The tolerances are declared before the comparison and are asserted, so a disagreement stops the script rather than being described.
Two kinds of disagreement are expected and distinguished: arithmetic that must match to floating-point, and convention differences (lambda search precision, ddof in a
standard deviation) that move a figure in its last digits. A difference that is neither is a finding, and the script says which kind each one is.
Usage: python3 scripts/95_independent_reimplementation.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
T = ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
rows = []
FAIL = []


def check(name, mine, theirs, tol, kind):
    """Record a comparison. `kind` is 'exact' for arithmetic that must agree, 'convention' for a declared tolerance."""
    d = abs(float(mine) - float(theirs))
    ok = d <= tol
    rows.append({"part": "comparison", "quantity": name, "independent": float(mine), "published": float(theirs),
                 "abs_diff": d, "tolerance": tol, "kind": kind, "agrees": bool(ok)})
    print(f"   {'OK ' if ok else 'XX '} {name:52s} mine {float(mine):+.6f}  published {float(theirs):+.6f}  diff {d:.2e} (tol {tol:.0e}, {kind})")
    if not ok:
        FAIL.append(name)
    return ok


# ----------------------------------------------------------------------------------------------------------------------
print("1. an independent Newick parser, because both existing VCV builders get their branch lengths from the same place")


def parse_newick(text):
    """Recursive descent over a Newick string -> (children, branch_length, name) nodes. Written here on purpose.
    Handles quoted labels, internal labels and support values, and a missing branch length (treated as 0)."""
    s = text.strip()
    assert s.endswith(";"), "newick must end with a semicolon"
    s = s[:-1]
    pos = 0

    def read_name():
        nonlocal pos
        if pos < len(s) and s[pos] in "'\"":
            q = s[pos]; pos += 1; start = pos
            while pos < len(s) and s[pos] != q:
                pos += 1
            name = s[start:pos]; pos += 1
            return name
        start = pos
        while pos < len(s) and s[pos] not in "(),:;":
            pos += 1
        return s[start:pos].strip()

    def read_length():
        nonlocal pos
        if pos < len(s) and s[pos] == ":":
            pos += 1; start = pos
            while pos < len(s) and (s[pos].isdigit() or s[pos] in ".eE+-"):
                pos += 1
            return float(s[start:pos])
        return 0.0

    def node():
        nonlocal pos
        kids = []
        if pos < len(s) and s[pos] == "(":
            pos += 1
            while True:
                kids.append(node())
                assert pos < len(s), "unbalanced parentheses in newick"
                if s[pos] == ",":
                    pos += 1
                    continue
                assert s[pos] == ")", f"expected ) or , at {pos}, got {s[pos]!r}"
                pos += 1
                break
        name = read_name()
        bl = read_length()
        return {"kids": kids, "bl": bl, "name": name}

    root = node()
    assert pos == len(s), f"trailing characters in newick at {pos}"
    return root


root = parse_newick((T / "Bacteria.nwk").read_text())


def tip_paths(node):
    """For every tip, the list of (node-id, cumulative depth) from the root down to it. No tree library involved."""
    out = {}
    counter = [0]

    def walk(nd, path, depth):
        nid = counter[0]; counter[0] += 1
        d = depth + nd["bl"]
        p = path + [(nid, d)]
        if not nd["kids"]:
            out[nd["name"]] = p
        for k in nd["kids"]:
            walk(k, p, d)

    walk(node, [], 0.0)
    return out


sys.setrecursionlimit(100000)
paths = tip_paths(root)
print(f"   parsed {len(paths)} tips independently of Bio.Phylo")
# Cross-check the PARSE itself against the library, on the one quantity both can produce without sharing a route:
# the root-to-tip depth of every tip. A mis-read branch length shows up here and nowhere else.
from Bio import Phylo  # noqa: E402  - imported only to cross-check the parse, never to produce a result below
lib_tree = Phylo.read(T / "Bacteria.nwk", "newick")
lib_depth = {t.name: lib_tree.distance(lib_tree.root, t) for t in lib_tree.get_terminals()}
dd = np.array([abs(paths[k][-1][1] - lib_depth[k]) for k in lib_depth])
print(f"   root-to-tip depths agree with Bio.Phylo to {dd.max():.3e} (max over {len(dd)} tips)")
rows.append({"part": "1", "quantity": "max root-to-tip depth difference vs Bio.Phylo", "independent": float(dd.max()),
             "published": 0.0, "abs_diff": float(dd.max()), "tolerance": 1e-9, "kind": "exact", "agrees": bool(dd.max() < 1e-9)})
assert dd.max() < 1e-9, f"the independent parser disagrees with Bio.Phylo on branch lengths by {dd.max():.3e}; everything below is suspect"

# ----------------------------------------------------------------------------------------------------------------------
print("\n2. the data and the predictor, rebuilt from the stored panel without the helper that built them")
pan = pd.read_csv(D / "structure_panel_wide.csv").drop_duplicates("ncbi_assembly")
pheno = (pd.read_csv(D / "phenotypes_all.csv")
         .merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly"))
refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv")
MEM = [str(r.label) for r in refs[refs.accepted].itertuples() if r.membrane]
cls_counts = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
B = pheno[(pheno.domain == "Bacteria") & pheno.ncbi_assembly.isin(paths)].drop_duplicates("ncbi_assembly").reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100.0
SRCC = [f"src_{o}" for o in SRC]
M = pan.set_index("ncbi_assembly")


def composite_mine(cls, labels, minrefs=3):
    """The published composite's definition, implemented here: per reference, net charge share (D+E-K-R)/positions under a
    position floor of max(8, half the reference's class length); z-score each across genomes; average, requiring minrefs."""
    z = {}
    for r in labels:
        cn = f"{r}|{cls}|n"
        if cn not in M.columns:
            continue
        n_ = M[cn].reindex(B.ncbi_assembly).to_numpy(float)
        num = (M[f"{r}|{cls}|D"] + M[f"{r}|{cls}|E"] - M[f"{r}|{cls}|K"] - M[f"{r}|{cls}|R"]).reindex(B.ncbi_assembly).to_numpy(float)
        floor = max(8.0, 0.5 * float(cls_counts[(r, cls)])) if (r, cls) in cls_counts.index else 8.0
        x = np.where(n_ >= floor, num / np.where(n_ > 0, n_, np.nan), np.nan)
        mu, sd = np.nanmean(x), np.nanstd(x, ddof=1)      # ddof=1: pandas Series.std(), which the published helper uses
        z[r] = (x - mu) / sd
    Z = np.column_stack([z[k] for k in z])
    cnt = np.sum(~np.isnan(Z), axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = np.nanmean(Z, axis=1)
    return np.where(cnt >= minrefs, s, np.nan), len(z)


x_peri, k_peri = composite_mine("peri_exp", MEM)
print(f"   rebuilt the periplasmic composite over {k_peri} membrane references; {int(np.sum(~np.isnan(x_peri)))} genomes scorable")

# ----------------------------------------------------------------------------------------------------------------------
print("\n3. the covariance matrix, from shared root-to-tip path length by prefix comparison")


def vcv_mine(tips):
    """C[i,j] = total branch length shared by the root-to-i and root-to-j paths = depth of their last common node.
    Computed by comparing ancestor-id sequences, with no tree walk and no distance call."""
    P = [paths[t] for t in tips]
    md = max(len(p) for p in P)
    ids = np.full((len(P), md), -1, dtype=np.int64)
    dep = np.zeros((len(P), md))
    for i, p in enumerate(P):
        ids[i, :len(p)] = [a for a, _ in p]
        dep[i, :len(p)] = [d for _, d in p]
        dep[i, len(p):] = p[-1][1]
    C = np.empty((len(P), len(P)))
    for i in range(len(P)):
        same = ids == ids[i]                       # where j's ancestor chain still matches i's
        shared = np.cumprod(same, axis=1).sum(axis=1).astype(int)   # length of the common prefix
        C[i] = dep[np.arange(len(P)), np.maximum(shared - 1, 0)]
    return C


keep = ~np.isnan(x_peri) & B.ph_min.notna().to_numpy()
tips_used = list(B.ncbi_assembly[keep])
Cm = vcv_mine(tips_used)
# Cross-check the MATRIX against the library's builder. They must agree to floating point: same object, different route.
sys.path.insert(0, str(ROOT / "src"))
from acidophile import phylo as ph  # noqa: E402  - cross-check only; no result below is computed with it
C_lib = ph.vcv_fast(lib_tree, tips_used)
md = float(np.abs(Cm - C_lib).max())
print(f"   {Cm.shape[0]}x{Cm.shape[0]} matrix; max element-wise difference from ph.vcv_fast = {md:.3e}")
check("covariance matrix vs ph.vcv_fast (max element)", md, 0.0, 1e-9, "exact")

# ----------------------------------------------------------------------------------------------------------------------
print("\n4. the estimator, by a different numerical route end to end")


def pgls_mine(y, X, C, lam=None, grid=4001):
    """Pagel's-lambda GLS by symmetric eigendecomposition, explicit normal equations and a dense lambda grid.
    THE ONE EIGENDECOMPOSITION THAT MAKES A 4001-POINT GRID AFFORDABLE, and it is exact rather than an approximation.
    Pagel's lambda scales the off-diagonals only, so V(l) = l*C + (1-l)*D with D = diag(diag(C)). Factor the diagonal out:
        V(l) = D^(1/2) [ l*R + (1-l)*I ] D^(1/2),   R = D^(-1/2) C D^(-1/2)
    and l*R + (1-l)*I has the SAME eigenvectors as R for every l, with eigenvalues l*mu_i + (1-l). So one eigendecomposition
    of R serves the whole sweep: V^-1 and log|V| are then O(n^2) per lambda instead of O(n^3). The first draft of this
    script eigendecomposed V at each of 4001 grid points, which for n=2,768 does not finish, and the temptation at that
    point is to fall back on a Cholesky sweep - which is the one route this script exists not to share.
    Nothing here touches Cholesky, lstsq, Brent search, or a Cholesky log-diagonal determinant."""
    y = np.asarray(y, float); X = np.asarray(X, float)
    n, k = X.shape
    C = C / (float(np.mean(np.diag(C))) or 1.0)       # the published model's normalisation; part of the MODEL, not the code
    d = np.diag(C).copy()
    si = 1.0 / np.sqrt(d)
    R = C * si[:, None] * si[None, :]
    mu, Q = np.linalg.eigh(R)                          # the single decomposition; symmetric, not Cholesky
    Xs = Q.T @ (X * si[:, None])                       # whitened by D^(-1/2) then rotated into R's eigenbasis
    ys = Q.T @ (y * si)
    logdet_d = float(np.sum(np.log(d)))

    def ll_of(l):
        w = l * mu + (1.0 - l)
        w = np.maximum(w, 1e-12)
        iw = 1.0 / w
        A = Xs.T @ (Xs * iw[:, None])                  # X' V^-1 X in the rotated basis
        beta = np.linalg.solve(A, Xs.T @ (ys * iw))    # explicit normal equations
        r = ys - Xs @ beta
        s2 = float(np.sum(r * r * iw) / (n - k))
        ll = -0.5 * (n * np.log(2 * np.pi * s2) + (logdet_d + float(np.sum(np.log(w)))) + (n - k))
        return ll, beta, s2, A

    if lam is None:
        gs = np.linspace(0.0, 1.0, grid)
        lls = np.array([ll_of(float(gg))[0] for gg in gs])
        lam = float(gs[int(np.argmax(lls))])
    ll, beta, s2, A = ll_of(lam)
    cov = s2 * np.linalg.inv(A)
    se = np.sqrt(np.diag(cov))
    t = beta / se
    return {"beta": beta, "se": se, "p": 2 * stats.t.sf(np.abs(t), n - k), "lambda": lam, "loglik": ll, "n": n}


y = B.ph_min.to_numpy(float)[keep]
xv = x_peri[keep]
g = B[keep]
use = [c for c in ["gw"] + SRCC if g[c].nunique() > 1]
X0 = np.column_stack([np.ones(len(y))] + [g[c].to_numpy(float) for c in use])
# The published headline comes from pgls_scan, whose convention is: lambda fitted ONCE on the covariates-only model, the
# tested predictor standardised internally so the returned coefficient is already per SD. Both are reproduced explicitly.
lam0 = pgls_mine(y, X0, Cm)["lambda"]
sd_pop = float(np.std(xv))                 # ddof=0, which is what scripts/57 divides by to get its per-unit figure
z = (xv - xv.mean()) / np.std(xv)
fit = pgls_mine(y, np.column_stack([X0, z]), Cm, lam=lam0)
b_sd, se_sd, p_sd = float(fit["beta"][-1]), float(fit["se"][-1]), float(fit["p"][-1])
print(f"   n {len(y)}, lambda {lam0:.4f} (grid of 4001), per SD {b_sd:+.6f} (se {se_sd:.6f}, p {p_sd:.3e}), per unit {b_sd / sd_pop:+.6f}")

print("\n   against docs/NUMBERS.md's headline row:")
pub = pd.read_csv(D / "headline_numbers.csv")
H = pub[pub.figure == "exposed periplasmic, 32-reference rule-selected panel"].iloc[0]
check("headline n", len(y), H.n, 0, "exact")
check("headline references used", k_peri, H.references_used, 0, "exact")
check("headline SD of predictor", sd_pop, H.sd_of_predictor, 5e-4, "convention (ddof)")
check("headline lambda", lam0, H.lam, 2e-3, "convention (grid vs Brent)")
check("headline beta per SD", b_sd, H.beta_per_sd, 5e-4, "convention")
check("headline se per SD", se_sd, H.se_per_sd, 5e-4, "convention")
check("headline beta per unit", b_sd / sd_pop, H.beta_per_unit, 2e-3, "convention")
# p-values span orders of magnitude, so they are compared on the log scale, which is how they are read.
check("headline log10(p)", np.log10(p_sd), np.log10(float(H.p)), 0.05, "convention")

print("\n   the per-SD/per-unit identity, checked on my own numbers rather than assumed:")
ident = abs((b_sd / sd_pop) * sd_pop - b_sd)
check("beta_per_unit * SD == beta_per_sd", ident, 0.0, 1e-12, "exact")

# ----------------------------------------------------------------------------------------------------------------------
print("\n5. the same, for the rows the claim actually rests on")
SPECS = [("exposed periplasmic, 32-reference panel, + growth optimum", "peri_exp", MEM, ["gw", "ph_opt"], ()),
         ("exposed periplasmic, 32 refs, + composition, on the optimum-known sample", "peri_exp", MEM, ["gw"], ("ph_opt",)),
         ("exposed periplasmic, 32 references, NO covariates", "peri_exp", MEM, [], ()),
         ("exposed cytoplasmic control, 32 references", "cyto_exp", MEM, ["gw"], ()),
         ("transmembrane control, 32 references", "TM", MEM, ["gw"], ())]
for label, cls, labels, covs, require in SPECS:
    xc, kc = composite_mine(cls, labels)
    m = (~np.isnan(xc)) & B.ph_min.notna().to_numpy()
    for c in list(covs) + list(require):
        m &= B[c].notna().to_numpy()
    gg = B[m]
    uu = [c for c in list(covs) + SRCC if gg[c].nunique() > 1]
    Xb = np.column_stack([np.ones(int(m.sum()))] + [gg[c].to_numpy(float) for c in uu])
    Cc = vcv_mine(list(gg.ncbi_assembly))
    yy = gg.ph_min.to_numpy(float)
    lm = pgls_mine(yy, Xb, Cc)["lambda"]
    vv = xc[m]
    f = pgls_mine(yy, np.column_stack([Xb, (vv - vv.mean()) / np.std(vv)]), Cc, lam=lm)
    P = pub[pub.figure == label].iloc[0]
    print(f"\n   {label}")
    check(f"{label[:40]} | n", int(m.sum()), P.n, 0, "exact")
    check(f"{label[:40]} | beta per SD", float(f["beta"][-1]), P.beta_per_sd, 5e-4, "convention")
    check(f"{label[:40]} | log10(p)", np.log10(float(f["p"][-1])), np.log10(float(P.p)), 0.05, "convention")

# ----------------------------------------------------------------------------------------------------------------------
print("\n6. the two sided results, which are the project's second and third claims")
# These come from scripts/82's pooled parametrisation, which is a different predictor construction from the composite
# above: counts summed ACROSS references, then divided, rather than z-scored per reference and averaged. Reproducing it
# is a separate exercise and is done here with the same independent estimator.
pan_det = pd.read_csv(D / "structure_panel_wide_det.csv").drop_duplicates("ncbi_assembly")
Md = pan_det.set_index("ncbi_assembly")
Bd = pheno[(pheno.domain == "Bacteria") & pheno.ncbi_assembly.isin(paths) & pheno.ph_min.notna()].drop_duplicates("ncbi_assembly").reset_index(drop=True)
for o in SRC:
    Bd[f"src_{o}"] = (Bd.origin == o).astype(float)
Bd["gw"] = Bd.acidic_excess * 100.0


def pooled_mine(cls, keys=("D", "E", "H", "K", "R")):
    tot = {k: np.zeros(len(Bd)) for k in ("n",) + tuple(keys)}
    for r in MEM:
        if f"{r}|{cls}|n" not in Md.columns or f"{r}|{cls}|H" not in Md.columns:
            continue
        tot["n"] += np.nan_to_num(Md[f"{r}|{cls}|n"].reindex(Bd.ncbi_assembly).to_numpy(float))
        for k in keys:
            tot[k] += np.nan_to_num(Md[f"{r}|{cls}|{k}"].reindex(Bd.ncbi_assembly).to_numpy(float))
    return tot


tpp, tcc = pooled_mine("peri_exp"), pooled_mine("cyto_exp")
ok2 = (tpp["n"] >= 200) & (tcc["n"] >= 200)
Bd2 = Bd[ok2].reset_index(drop=True)
Pp = {k: 100.0 * tpp[k][ok2] / tpp["n"][ok2] for k in ("D", "E", "H", "K", "R")}
Cc2 = {k: 100.0 * tcc[k][ok2] / tcc["n"][ok2] for k in ("D", "E", "H", "K", "R")}
C2 = vcv_mine(list(Bd2.ncbi_assembly))
y2 = Bd2.ph_min.to_numpy(float)
uu2 = [c for c in ["gw"] + SRCC if Bd2[c].nunique() > 1]
Xb2 = np.column_stack([np.ones(len(y2))] + [Bd2[c].to_numpy(float) for c in uu2])
lam2 = pgls_mine(y2, Xb2, C2)["lambda"]
print(f"   n {len(y2)}, lambda {lam2:.4f}")
sid = pd.read_csv(D / "histidine_sidedness_wide.csv")
for lab, pp, cc, others, pubname in (("acidic D+E", Pp["D"] + Pp["E"], Cc2["D"] + Cc2["E"], (), "acidic D+E"),
                                     ("histidine", Pp["H"], Cc2["H"], ("D", "E", "K"), "histidine")):
    Xs = np.column_stack([Xb2] + [Pp[k] for k in others] + [pp + cc, pp - cc])
    f = pgls_mine(y2, Xs, C2, lam=lam2)
    for j, term in ((-2, "sum"), (-1, "difference")):
        r = sid[(sid.part == "contrast") & (sid.measure == pubname) & (sid.term == term)]
        if not len(r):
            continue
        check(f"{lab} {term} (per pt)", float(f["beta"][j]), float(r.iloc[0].beta_per_pt), 5e-4, "convention")

# ----------------------------------------------------------------------------------------------------------------------
print("\n7. verdict")
nex = sum(1 for r in rows if r.get("kind") == "exact")
ncv = sum(1 for r in rows if r.get("kind", "").startswith("convention"))
if FAIL:
    print(f"   {len(FAIL)} of {nex + ncv} comparisons DISAGREE beyond their declared tolerance:")
    for f_ in FAIL:
        print(f"      {f_}")
    print("   That is a finding about the estimator or the predictor's arithmetic, not a convention difference, and the")
    print("   figures it touches should not be quoted until it is resolved.")
else:
    print(f"   all {nex + ncv} comparisons agree: {nex} to floating point, {ncv} within a declared convention tolerance.")
    print("   The headline, its optimum-controlled form, its two matched controls and both sided results are reproduced by")
    print("   an implementation that shares no Newick parser, no covariance construction, no matrix factorisation and no")
    print("   lambda search with the one that produced them. Twenty corrections were found in the analyses and none in the")
    print("   estimator; this is the first evidence that the second fact is not merely an absence of looking.")
print("   WHAT THIS DOES NOT SAY. The structural panel, the HMM detection, the exposure calls and the phenotype curation")
print("   are inputs here, not reproductions. An error in any of them is invisible to this script and would be inherited")
print("   by every figure it just confirmed. The scope is the fitting and the predictor's arithmetic, and nothing wider.")

pd.DataFrame(rows).to_csv(D / "independent_reimplementation.csv", index=False)
print(f"\nwrote {D / 'independent_reimplementation.csv'}")
assert not FAIL, f"{len(FAIL)} comparison(s) disagree beyond tolerance: {FAIL}"
