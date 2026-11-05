"""Marker-gene tree and phylogeny-aware regression.

Tree: 29 ribosomal-protein Pfam families (gathering cutoffs), best hit per genome per family, aligned to the Pfam profile
(pyhmmer hmmalign, match columns only), concatenated, columns with < 50 % occupancy and genomes with < 50 % of families
removed, gamma-corrected Poisson distances on shared sites, neighbour-joining, midpoint root. This is NOT a maximum-likelihood
tree (no ML software is installed and none was downloaded); topology uncertainty is carried by column-bootstrap trees, and
every phylogenetic test is repeated over those trees.
Regression: PGLS with Pagel's lambda fitted by maximum likelihood."""
from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import pandas as pd
import pyhmmer
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from scipy import stats
from scipy.optimize import minimize_scalar

MARKERS = ("PF00181 PF03947 PF00297 PF00573 PF00673 PF00347 PF00238 PF00252 PF00861 PF00237 PF00276 PF00189 PF00333 "
           "PF03719 PF00410 PF00338 PF00366 PF00203 PF00828 PF00572 PF00177 PF00164 PF00411 PF00380 PF00318 PF00298 "
           "PF00687 PF00466 PF00416").split()
_s = lambda x: x.decode() if isinstance(x, bytes) else str(x)


def best_marker_hits(proteome: Path, hmms) -> dict:
    """family -> (protein name, sequence text) with the highest-scoring gathering-cutoff hit."""
    alpha = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(proteome, digital=True, alphabet=alpha) as sf:
        seqs = sf.read_block()
    by_name = {_s(s.name): s for s in seqs}
    best = {}
    for top in pyhmmer.hmmsearch(hmms, seqs, bit_cutoffs="gathering", cpus=4):
        fam = _s(top.query.accession).split(".")[0]
        for h in top:
            if h.included and (fam not in best or h.score > best[fam][0]):
                best[fam] = (h.score, _s(h.name))
    return {fam: (n, by_name[n].textize().sequence) for fam, (_, n) in best.items()}


def align_family(hmm, named_seqs: dict) -> dict:
    """genome -> aligned match-column string (upper case / '-') for one family."""
    alpha = pyhmmer.easel.Alphabet.amino()
    block = pyhmmer.easel.TextSequenceBlock([pyhmmer.easel.TextSequence(name=g.encode(), sequence=s) for g, s in named_seqs.items()]).digitize(alpha)
    msa = pyhmmer.hmmalign(hmm, block, trim=True)
    out = {}
    for name, seq in zip(msa.names, msa.alignment):
        keep = [c.upper() if c.isalpha() and c.isupper() else "-" for c in seq if c not in ".abcdefghijklmnopqrstuvwxyz"]
        out[_s(name)] = "".join(keep)
    return out


def concatenate(per_family: dict, genomes: list, min_family_frac=0.5, min_col_occ=0.5):
    """per_family: fam -> {genome: aligned str}. Returns (genome x column character matrix as np.array, genome list)."""
    lens = {f: len(next(iter(a.values()))) for f, a in per_family.items() if a}
    cols = []
    for g in genomes:
        cols.append("".join(per_family[f].get(g, "-" * lens[f]) for f in sorted(lens)))
    M = np.array([list(c) for c in cols])
    frac = np.array([sum(g in per_family[f] for f in lens) / len(lens) for g in genomes])
    M = M[frac >= min_family_frac]
    kept = [g for g, f in zip(genomes, frac) if f >= min_family_frac]
    occ = (M != "-").mean(axis=0)
    return M[:, occ >= min_col_occ], kept


def distance_matrix(M: np.ndarray, alpha=1.0):
    """Gamma-corrected Poisson distance on shared non-gap sites."""
    n = M.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            ok = (M[i] != "-") & (M[j] != "-")
            if ok.sum() < 30:
                D[i, j] = D[j, i] = 5.0
                continue
            p = min((M[i][ok] != M[j][ok]).mean(), 0.93)
            D[i, j] = D[j, i] = alpha * ((1 - p) ** (-1 / alpha) - 1)
    return D


def nj_tree(D: np.ndarray, names: list):
    lower = [list(D[i, :i + 1]) for i in range(len(names))]
    tree = DistanceTreeConstructor().nj(DistanceMatrix(names, lower))
    for c in tree.find_clades():
        if c.branch_length is not None and c.branch_length < 0:
            c.branch_length = 0.0
    tree.root_at_midpoint()
    return tree


def bootstrap_trees(M, names, B=100, seed=0):
    rng = np.random.default_rng(seed)
    trees = []
    for _ in range(B):
        cols = rng.integers(0, M.shape[1], M.shape[1])
        trees.append(nj_tree_fast(distance_matrix_fast(M[:, cols]), names))
    return trees


def to_newick(tree) -> str:
    buf = io.StringIO(); Phylo.write(tree, buf, "newick"); return buf.getvalue().strip()


def vcv_from_tree(tree, tips: list) -> np.ndarray:
    """Brownian covariance for the induced subtree of `tips`, rooted at their MRCA."""
    clades = {c.name: c for c in tree.get_terminals()}
    mrca = tree.common_ancestor([clades[t] for t in tips]) if len(tips) > 1 else clades[tips[0]]
    d = np.array([tree.distance(mrca, clades[t]) for t in tips])
    n = len(tips)
    C = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            C[i, j] = C[j, i] = d[i] if i == j else (d[i] + d[j] - tree.distance(clades[tips[i]], clades[tips[j]])) / 2
    return C


def pgls(y, X, C, lam=None):
    """GLS with Pagel's lambda (off-diagonals of C scaled). Returns dict with beta, se, t, p, lambda, loglik. X includes intercept."""
    y = np.asarray(y, float); X = np.asarray(X, float); n, k = X.shape
    scale = np.mean(np.diag(C)) or 1.0
    C = C / scale

    def fit(l):
        V = C * l + np.diag(np.diag(C)) * (1 - l)
        L = np.linalg.cholesky(V + 1e-10 * np.eye(n))
        Xs, ys = np.linalg.solve(L, X), np.linalg.solve(L, y)
        beta, *_ = np.linalg.lstsq(Xs, ys, rcond=None)
        r = ys - Xs @ beta
        s2 = float(r @ r / (n - k))
        ll = -0.5 * (n * np.log(2 * np.pi * s2) + 2 * np.log(np.diag(L)).sum() + (n - k))
        return beta, s2, ll, Xs

    if lam is None:
        res = minimize_scalar(lambda l: -fit(l)[2], bounds=(0.0, 1.0), method="bounded")
        lam = float(res.x)
    beta, s2, ll, Xs = fit(lam)
    cov = s2 * np.linalg.inv(Xs.T @ Xs)
    se = np.sqrt(np.diag(cov)); t = beta / se
    p = 2 * stats.t.sf(np.abs(t), n - k)
    return {"beta": beta, "se": se, "t": t, "p": p, "lambda": lam, "loglik": ll, "sigma2": s2, "n": n}


def vcv_fast(tree, tips: list) -> np.ndarray:
    """Same matrix as vcv_from_tree for tips whose MRCA is the tree root, in O(n^2) array operations (single traversal)."""
    index = {t: i for i, t in enumerate(tips)}
    n = len(tips)
    C = np.zeros((n, n))

    def walk(clade, depth):
        if clade.is_terminal():
            i = index.get(clade.name)
            if i is not None:
                C[i, i] = depth
            return [] if i is None else [i]
        groups = [walk(c, depth + (c.branch_length or 0.0)) for c in clade.clades]
        groups = [g for g in groups if g]
        for a in range(len(groups)):
            for b in range(a + 1, len(groups)):
                C[np.ix_(groups[a], groups[b])] = depth
                C[np.ix_(groups[b], groups[a])] = depth
        return [i for g in groups for i in g]

    walk(tree.root, 0.0)
    return C


def distance_matrix_fast(M: np.ndarray, alpha=1.0, min_shared=30):
    """Vectorised equivalent of distance_matrix (gamma-corrected Poisson distance on shared non-gap sites)."""
    codes = np.unique(M, return_inverse=True)[1].reshape(M.shape).astype(np.int16)
    gap = np.searchsorted(np.unique(M), "-")
    present = (codes != gap)
    n = M.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        ok = present & present[i]
        shared = ok.sum(axis=1)
        diff = ((codes != codes[i]) & ok).sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            p = np.minimum(diff / np.where(shared > 0, shared, 1), 0.93)
        d = alpha * ((1 - p) ** (-1 / alpha) - 1)
        d[shared < min_shared] = 5.0
        D[i] = d
    D = (D + D.T) / 2
    np.fill_diagonal(D, 0.0)
    return D


def nj_tree_fast(D: np.ndarray, names: list):
    """Neighbour-joining (Saitou & Nei) with numpy; returns a midpoint-rooted Biopython tree. Equivalent topology to nj_tree."""
    from Bio.Phylo.BaseTree import Clade, Tree
    n = len(names)
    D = D.astype(float).copy()
    clades = [Clade(branch_length=0.0, name=nm) for nm in names]
    active = list(range(n))
    while len(active) > 3:
        sub = D[np.ix_(active, active)]
        m = len(active)
        r = sub.sum(axis=1)
        Q = (m - 2) * sub - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = np.unravel_index(np.argmin(Q), Q.shape)
        a, b = active[i], active[j]
        dab = D[a, b]
        la = 0.5 * dab + (r[i] - r[j]) / (2 * (m - 2))
        lb = dab - la
        clades[a].branch_length, clades[b].branch_length = max(la, 0.0), max(lb, 0.0)
        u = Clade(branch_length=0.0, clades=[clades[a], clades[b]])
        newd = 0.5 * (D[a, :] + D[b, :] - dab)
        D[a, :] = newd; D[:, a] = newd; D[a, a] = 0.0       # reuse slot `a` for the new node
        clades[a] = u
        active.pop(j)
    a, b, c = active
    dab, dac, dbc = D[a, b], D[a, c], D[b, c]
    la, lb, lc = 0.5 * (dab + dac - dbc), 0.5 * (dab + dbc - dac), 0.5 * (dac + dbc - dab)
    for k, l in ((a, la), (b, lb), (c, lc)):
        clades[k].branch_length = max(l, 0.0)
    tree = Tree(root=Clade(clades=[clades[a], clades[b], clades[c]]), rooted=False)
    tree.root_at_midpoint()
    return tree


def pgls_scan(C, y, Xbase, F, lam=None):
    """Closed-form PGLS scan: lambda fitted once on the covariates-only model (or given), then each column of F (standardised inside) is tested as an
    added regressor via Frisch-Waugh on the whitened design. Returns arrays (beta, se, p) and lambda. Matches pgls() up to the fixed-lambda approximation.

    SCALE: F is standardised inside this function, so the returned beta is ALREADY per standard deviation of the column as passed in.
    Do NOT multiply it by that column's standard deviation to "convert to per SD" - that converts it BACK to per unit. This mistake was
    made four times in scripts/ before being caught, each time understating a reported effect by a factor of the predictor's SD."""
    from scipy.linalg import solve_triangular
    y = np.asarray(y, float); Xbase = np.asarray(Xbase, float); F = np.asarray(F, float)
    if lam is None:
        lam = pgls(y, Xbase, C)["lambda"]
    Cn = C / np.mean(np.diag(C))
    V = Cn * lam + np.diag(np.diag(Cn)) * (1 - lam)
    L = np.linalg.cholesky(V + 1e-10 * np.eye(len(y)))
    Xs = solve_triangular(L, Xbase, lower=True); ys = solve_triangular(L, y, lower=True)
    Q, _ = np.linalg.qr(Xs)
    yt = ys - Q @ (Q.T @ ys)
    sd = F.std(axis=0); ok = sd > 0
    Fz = (F - F.mean(axis=0)) / np.where(ok, sd, 1)
    Fs = solve_triangular(L, Fz, lower=True); Ft = Fs - Q @ (Q.T @ Fs)
    xx = (Ft ** 2).sum(axis=0); xx[~ok] = np.nan
    beta = (Ft * yt[:, None]).sum(axis=0) / xx
    dfree = len(y) - Xs.shape[1] - 1
    s2 = ((yt ** 2).sum() - beta ** 2 * xx) / dfree
    se = np.sqrt(s2 / xx)
    p = 2 * stats.t.sf(np.abs(beta / se), dfree)
    return beta, se, p, lam


def pic_contrasts(tree, tips: list, traits: np.ndarray, eps=1e-8):
    """Felsenstein phylogenetic independent contrasts for several traits at once on the subtree induced by `tips`.
    traits: (n_tips, k) aligned with `tips`. Returns (contrasts (n-1, k), is_cherry (n-1,) bool): a cherry contrasts two sister tips,
    i.e. the closest relatives. Polytomies are resolved arbitrarily with zero-length internal branches."""
    import copy
    import sys
    sys.setrecursionlimit(max(10000, sys.getrecursionlimit()))
    keep = set(tips)
    t2 = copy.deepcopy(tree)
    for leaf in [x for x in t2.get_terminals() if x.name not in keep]:
        t2.prune(leaf)
    idx = {name: i for i, name in enumerate(tips)}
    traits = np.asarray(traits, float)
    cons, cher = [], []

    def merge(a, b):
        (xa, ea, ta), (xb, eb, tb) = a, b
        tot = ea + eb
        cons.append((xa - xb) / np.sqrt(tot)); cher.append(ta and tb)
        return (xa * eb + xb * ea) / tot, ea * eb / tot, False

    def walk(c):
        if c.is_terminal():
            return traits[idx[c.name]], max(c.branch_length or 0.0, eps), True
        kids = [walk(k) for k in c.clades]
        while len(kids) > 1:
            kids = [merge(kids[0], kids[1])] + kids[2:]
        x, e, t = kids[0]
        return x, e + max(c.branch_length or 0.0, 0.0), t

    walk(t2.root)
    return np.array(cons), np.array(cher)


def regress_through_origin(y, X):
    """OLS without intercept (the correct regression for contrasts). Returns beta, se, p."""
    y = np.asarray(y, float); X = np.asarray(X, float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta; dfree = len(y) - X.shape[1]
    cov = (r @ r / dfree) * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    return beta, se, 2 * stats.t.sf(np.abs(beta / se), dfree)
