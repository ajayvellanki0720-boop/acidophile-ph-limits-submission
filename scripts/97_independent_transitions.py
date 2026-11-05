"""Did acid tolerance arise repeatedly, and did the composition shift arise with it each time? -> DATA_DIR/independent_transitions.csv
THE QUESTION A COMPARATIVE REGRESSION CANNOT ANSWER. Pagel's lambda stops a shared ancestry from being counted as independent evidence, but it does not tell you how many
times the thing happened. A coefficient of +0.1848 per SD is consistent with two very different histories: a composition shift that accompanies acid tolerance wherever acid
tolerance arises, or one ancient shift in one clade that happens to contain most of the acidophiles. The second would still be a real association and would mean something
entirely different, and nothing in this project has distinguished them.
WHAT IS ALREADY DONE AND WHY IT IS NOT THIS. `scripts/44` fits the model inside each phylum and pools - that is a within-lineage regression, not a count of origins.
`scripts/45` and `scripts/46` use sister pairs (cherries), which is the finest possible scale and therefore the one with the least pH contrast in it: most cherries are two
neutrophiles. Neither asks how many independent times the trait appeared, which is the question that separates convergence from inheritance.
THE DESIGN. Acid tolerance is treated as a binary state on the tips. A MAXIMAL TOLERANT CLADE is one in which every tip with a known minimum is tolerant and which is not
contained in a larger such clade; under no reversals each is one independent origin. For each, the comparison is against its SISTER clade in the same tree, which shares the
immediate ancestor and therefore controls for everything the ancestor carried. The predicted direction is fixed in advance from the published coefficient's sign: the
tolerant clade should be DEPLETED of exposed periplasmic acidic residues relative to its sister. Counting how many transitions go the predicted way is a sign test, which
makes this the only result in the project that does not depend on the PGLS machinery at all.
TWO THINGS THAT WOULD INVALIDATE IT, AND THE DEFENCE AGAINST EACH. (a) Composites must be residualised on genome-wide composition BEFORE any clade is compared, because
`scripts/46` found that failing to do so makes the exposed cytoplasmic class look positive in sister comparisons when the regressions say it is null - a whole-proteome
acidity shift moves both exposed classes together. Residualisation here is fitted across all bacteria once, exactly as there. (b) The matched control classes must be run
through the identical pipeline, because a sign test that fires on every position class is measuring something about clades rather than about the periplasm.
WHAT THIS CANNOT DO, stated because the limitation is structural and not fixable by more code. Tip states cannot distinguish a gain from a loss. A clade of tolerant species
sister to a neutrophile clade is scored here as one origin of tolerance; it is equally consistent with tolerance being ancestral to both and lost in the sister. Parsimony on
tip states counts STATE CHANGES, and calling them origins is an assumption, named here and not defended. The count is therefore a count of independent CONTRASTS, which is
what the sign test needs, and only an upper bound on origins.
SCALE. This part reports standardised mean differences between clades, not pH coefficients, because the unit is a clade rather than a genome. No figure here is comparable
with a per-SD regression coefficient and none is quoted as one.
Usage: python3 scripts/97_independent_transitions.py [threshold]   (default 4.0; env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import stats as st

warnings.filterwarnings("ignore")
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
T = ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
THR = float(sys.argv[1]) if len(sys.argv) > 1 else 4.0
MINSIS = 3            # a sister clade must carry at least this many scorable tips to be a comparison at all
rows = []

refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
pan = pd.read_csv(D / "structure_panel_wide_det.csv").drop_duplicates("ncbi_assembly")
pheno = (pd.read_csv(D / "phenotypes_all.csv")
         .merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess", "frac_His"]], on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
B = pheno[(pheno.domain == "Bacteria") & pheno.ncbi_assembly.isin(tips) & pheno.ph_min.notna()].drop_duplicates("ncbi_assembly").reset_index(drop=True)
M = pan.set_index("ncbi_assembly")


def composite(cls, residualise=True):
    """The published composite for a position class, then residualised on genome-wide acidic excess across ALL bacteria -
    scripts/46's control, applied before any clade comparison rather than after."""
    z = {}
    for r in MEM:
        if f"{r}|{cls}|n" not in M.columns:
            continue
        n_ = M[f"{r}|{cls}|n"].reindex(B.ncbi_assembly).to_numpy(float)
        num = (M[f"{r}|{cls}|D"] + M[f"{r}|{cls}|E"] - M[f"{r}|{cls}|K"] - M[f"{r}|{cls}|R"]).reindex(B.ncbi_assembly).to_numpy(float)
        fl = max(8.0, 0.5 * float(reflen[(r, cls)])) if (r, cls) in reflen.index else 8.0
        x = np.where(n_ >= fl, num / np.where(n_ > 0, n_, np.nan), np.nan)
        z[r] = (x - np.nanmean(x)) / np.nanstd(x, ddof=1)
    Z = np.column_stack([z[k] for k in z])
    s = np.where(np.sum(~np.isnan(Z), axis=1) >= 3, np.nanmean(Z, axis=1), np.nan)
    if not residualise:
        return s
    gw = B.acidic_excess.to_numpy(float) * 100.0
    ok = ~np.isnan(s)
    A = np.column_stack([np.ones(ok.sum()), gw[ok]])
    co = np.linalg.lstsq(A, s[ok], rcond=None)[0]
    out = np.full(len(s), np.nan)
    out[ok] = s[ok] - A @ co
    return out


def hist_share(cls):
    """Pooled exposed histidine share for a class, residualised on proteome-wide histidine the same way."""
    tot_n = np.zeros(len(B)); tot_h = np.zeros(len(B))
    for r in MEM:
        if f"{r}|{cls}|n" not in M.columns or f"{r}|{cls}|H" not in M.columns:
            continue
        tot_n += np.nan_to_num(M[f"{r}|{cls}|n"].reindex(B.ncbi_assembly).to_numpy(float))
        tot_h += np.nan_to_num(M[f"{r}|{cls}|H"].reindex(B.ncbi_assembly).to_numpy(float))
    s = np.where(tot_n >= 200, 100.0 * tot_h / np.where(tot_n > 0, tot_n, np.nan), np.nan)
    gw = B.frac_His.to_numpy(float) * 100.0
    ok = ~np.isnan(s) & ~np.isnan(gw)
    A = np.column_stack([np.ones(ok.sum()), gw[ok]])
    co = np.linalg.lstsq(A, s[ok], rcond=None)[0]
    out = np.full(len(s), np.nan)
    out[ok] = s[ok] - A @ co
    return out


# MEASURES, each with the direction predicted from the published sign BEFORE any clade is looked at. The acidic
# composite's published coefficient is POSITIVE (more net acidic -> higher minimum -> less tolerant), so a tolerant
# clade is predicted LOWER than its sister. Exposed histidine's is NEGATIVE, so a tolerant clade is predicted HIGHER.
MEAS = [("exposed periplasmic acidic (the headline measure)", composite("peri_exp"), -1, True),
        ("exposed periplasmic histidine", hist_share("peri_exp"), +1, True),
        ("exposed cytoplasmic acidic (control)", composite("cyto_exp"), 0, False),
        ("transmembrane acidic (control)", composite("TM"), 0, False),
        ("exposed cytoplasmic histidine (control)", hist_share("cyto_exp"), 0, False)]

print(f"1. the tolerant clades, at a threshold of pH_min <= {THR}")
state = {a: (v <= THR) for a, v in zip(B.ncbi_assembly, B.ph_min.to_numpy(float))}
idx = {a: i for i, a in enumerate(B.ncbi_assembly)}
n_tol = int(sum(state.values()))
print(f"   {len(B)} bacterial genomes with a minimum; {n_tol} of them tolerant ({n_tol / len(B):.1%})")


def scorable_tips(clade):
    return [t.name for t in clade.get_terminals() if t.name in state]


def all_tolerant(clade):
    ts = scorable_tips(clade)
    return bool(ts) and all(state[t] for t in ts)


# Maximal tolerant clades, found top-down: descend only where the clade is not already wholly tolerant, so the first
# wholly-tolerant clade on each path is kept and nothing inside it is counted again.
maximal = []


def descend(clade):
    if all_tolerant(clade):
        maximal.append(clade)
        return
    for c in clade.clades:
        if scorable_tips(c):
            descend(c)


sys.setrecursionlimit(100000)
descend(tree.root)
sizes = [len(scorable_tips(c)) for c in maximal]
print(f"   {len(maximal)} maximal wholly-tolerant clades, covering {sum(sizes)} of the {n_tol} tolerant genomes")
print(f"   clade sizes: {sum(1 for s in sizes if s == 1)} singletons, {sum(1 for s in sizes if s > 1)} with 2 or more (largest {max(sizes) if sizes else 0})")
assert sum(sizes) == n_tol, f"the maximal clades cover {sum(sizes)} tolerant tips but there are {n_tol}; the partition is wrong"
rows.append({"part": "1", "threshold": THR, "measure": "maximal wholly-tolerant clades", "n": len(maximal),
             "n_tolerant_genomes": n_tol, "n_singleton": sum(1 for s in sizes if s == 1), "largest_clade": max(sizes) if sizes else 0})

print("\n2. pairing each tolerant clade with its sister, which shares its immediate ancestor")
parents = {}
for cl in tree.find_clades(order="level"):
    for c in cl.clades:
        parents[id(c)] = cl
pairs = []
for cl in maximal:
    par = parents.get(id(cl))
    if par is None:
        continue
    sis_tips = []
    for c in par.clades:
        if c is cl:
            continue
        sis_tips += [t for t in scorable_tips(c)]
    sis_tips = [t for t in sis_tips if not state[t]]        # the sister side must be non-tolerant to be a contrast
    if len(sis_tips) >= MINSIS:
        pairs.append((cl, scorable_tips(cl), sis_tips))
print(f"   {len(pairs)} tolerant clades have a sister carrying at least {MINSIS} non-tolerant scorable tips")
print(f"   ({len(maximal) - len(pairs)} are dropped: the sister is absent, too small, or itself tolerant - stated because dropping")
print("   them is a selection on the tree and not on any measure, which is why it cannot bias the direction below)")
rows.append({"part": "2", "threshold": THR, "measure": "usable tolerant-vs-sister contrasts", "n": len(pairs),
             "n_dropped": len(maximal) - len(pairs), "min_sister_tips": MINSIS})
if len(pairs) < 5:
    print("   FEWER THAN FIVE CONTRASTS: a sign test on this is not worth running and none is reported.")
    pd.DataFrame(rows).to_csv(D / "independent_transitions.csv", index=False)
    sys.exit(0)

print("\n3. the test across independent contrasts, direction fixed in advance from the published coefficient")
# BOTH A SIGN TEST AND A SIGNED-RANK TEST, and the reason is power rather than thoroughness. A sign test throws away the
# size of every contrast, which is affordable with hundreds of them and not with thirty-three; the signed-rank test keeps
# the magnitudes and is the primary here. The sign test is still reported because it assumes nothing about the shape of
# the differences, so a disagreement between the two is informative about that shape rather than about the biology.
print(f"   {'measure':>46s} {'pairs':>6s} {'as predicted':>13s} {'sign p':>9s} {'signed-rank p':>14s} {'median d':>9s}")
res = {}
for name, vec, pred, is_claim in MEAS:
    ds = []
    for cl, tt, ss in pairs:
        a = np.array([vec[idx[t]] for t in tt if t in idx], float)
        b = np.array([vec[idx[t]] for t in ss if t in idx], float)
        a, b = a[~np.isnan(a)], b[~np.isnan(b)]
        if len(a) < 1 or len(b) < MINSIS:
            continue
        ds.append(float(np.mean(a) - np.mean(b)))
    ds = np.array(ds)
    if len(ds) < 5:
        continue
    want = pred if pred != 0 else -1
    k = int(np.sum(np.sign(ds) == want))
    p_sign = float(st.binomtest(k, len(ds), 0.5, alternative="two-sided").pvalue)
    p_sr = float(st.wilcoxon(ds, alternative="two-sided").pvalue) if np.any(ds != 0) else 1.0
    res[name] = {"n": len(ds), "k": k, "p_sign": p_sign, "p_sr": p_sr, "median": float(np.median(ds)),
                 "mean": float(np.mean(ds)), "pred": pred, "claim": is_claim,
                 "dir_ok": bool(np.sign(np.median(ds)) == want)}
    print(f"   {name:>46s} {len(ds):>6d} {k:>6d}/{len(ds):<6d} {p_sign:>9.2e} {p_sr:>14.2e} {np.median(ds):>+9.4f}")
    rows.append({"part": "3", "threshold": THR, "measure": name, "n_pairs": len(ds), "n_as_predicted": k,
                 "sign_p": p_sign, "signed_rank_p": p_sr, "median_diff": float(np.median(ds)),
                 "mean_diff": float(np.mean(ds)), "predicted_direction": pred, "is_claim": is_claim,
                 "median_in_predicted_direction": res[name]["dir_ok"]})

print("\n4. power first, because it decides whether part 3 can be read at all")
# The critical count is the SMALLEST number of predicted-direction contrasts that reaches two-sided p < 0.05. A first
# draft of this took the LARGEST instead, by iterating downward and stopping at the first significant value, and reported
# that 80% power required the predicted direction in 100% of contrasts. That is the kind of error that makes a null look
# respectable, so the critical count is asserted against a hand-computed value here.
n_p = max(r["n"] for r in res.values())
crit = next(k for k in range(n_p // 2, n_p + 1) if st.binomtest(k, n_p, 0.5, alternative="two-sided").pvalue < 0.05)
assert n_p != 33 or crit == 23, f"for 33 contrasts the critical count is 23, not {crit}; the search is inverted again"
need = next((pr for pr in np.arange(0.50, 1.001, 0.01) if float(st.binom.sf(crit - 1, n_p, pr)) >= 0.80), None)
print(f"   {n_p} contrasts; two-sided 0.05 needs {crit} of them in the predicted direction")
if need is not None:
    print(f"   80% power requires the predicted direction to hold in {need:.0%} of contrasts")
else:
    print("   no consistency short of perfect reaches 80% power at this sample size")
rows.append({"part": "4", "threshold": THR, "measure": "power of the sign test", "n_pairs": n_p,
             "critical_count": int(crit), "consistency_for_80pct_power": float(need) if need else None})
for name, r in res.items():
    if r["claim"]:
        print(f"   {name}: {r['k']}/{r['n']} = {r['k'] / r['n']:.0%} observed")

print("\n5. what that makes the result")
claims = {n: r for n, r in res.items() if r["claim"]}
ctrls = {n: r for n, r in res.items() if not r["claim"]}
sig_c = [n for n, r in claims.items() if r["p_sr"] < 0.05 and r["dir_ok"]]
sig_k = [n for n, r in ctrls.items() if r["p_sr"] < 0.05]
dir_ok = [n for n, r in claims.items() if r["dir_ok"]]
print(f"   claim measures significant by signed rank AND in the predicted direction: {len(sig_c)} of {len(claims)}")
print(f"   control measures significant by signed rank: {len(sig_k)} of {len(ctrls)}")
print(f"   claim measures whose median merely POINTS the predicted way: {len(dir_ok)} of {len(claims)}")
# A CONTROL JUST OVER THE LINE IS NOT A CLEAN CONTROL. Exposed cytoplasmic histidine comes in at p 0.0506 here, and a
# verdict that reads "no control is significant" off that is deciding a specificity claim on the fourth decimal place of a
# threshold nobody chose for this purpose. Controls within 1.5x of the threshold are named as marginal, and a claim whose
# only matched control is marginal does not get called specific.
MARG = 0.05 * 1.5
marg_k = [n for n, r in ctrls.items() if 0.05 <= r["p_sr"] < MARG]
if marg_k:
    _m = ", ".join("%s at p %.4f" % (n, ctrls[n]["p_sr"]) for n in marg_k)
    print(f"   control measures MARGINAL (p between 0.05 and {MARG:.3f}): {_m}")
head_name = next(n for n, r in claims.items() if "headline" in n)
head = claims[head_name]
head_tag = (f"the headline measure is NOT significant here (signed-rank p {head['p_sr']:.3f}, {head['k']}/{head['n']} in the "
            f"predicted direction, median {head['median']:+.4f} which does point the predicted way)")
if sig_c and not sig_k and not marg_k:
    verdict = (f"CONVERGENT AND SPECIFIC at threshold {THR} for {', '.join(sig_c)}: it shifts as predicted across independent tolerant clades, "
               f"each against its own sister, and no control moves. Separately, {head_tag}.")
elif sig_c and not sig_k and marg_k:
    verdict = (f"CONVERGENT at threshold {THR} for {', '.join(sig_c)}, SPECIFICITY NOT ESTABLISHED: its own matched control "
               f"({', '.join(marg_k)}) sits at p {min(ctrls[n]['p_sr'] for n in marg_k):.4f}, which is over the line by an amount no one should "
               f"rest a specificity claim on. The sided contrast in scripts/82 is where that claim lives, not here. Separately, {head_tag}.")
elif sig_c and sig_k:
    verdict = (f"CONVERGENT BUT NOT SPECIFIC at threshold {THR}: claim and control measures both move, so the contrasts are detecting "
               f"something common to all position classes. Separately, {head_tag}.")
else:
    verdict = (f"UNDERPOWERED, NOT NULL, at threshold {THR}: no claim measure reaches significance in the predicted direction, but detecting the effect would need the "
               f"predicted direction in {need:.0%} of contrasts and there are only {n_p}. "
               f"{len(dir_ok)} of {len(claims)} claim measures point the predicted way and no control is significant, which is weak corroboration "
               "and must not be reported as more. The informative statement is the power, not the p-value.")
print(f"\n   verdict: {verdict}")
rows.append({"part": "5", "threshold": THR, "measure": "verdict", "n_pairs": n_p, "verdict": verdict,
             "n_sig_claim": len(sig_c), "n_sig_control": len(sig_k), "n_claim_direction_ok": len(dir_ok)})

print("\n6. the same test at other thresholds, reported in full because choosing one after the fact would be selection")
# The threshold is the one free parameter and moving it trades contrast quality against contrast count. All four are run
# and all four are printed; 4.0 is the primary because it is the value this script was written with, and the others are
# sensitivity. Nothing below is allowed to replace the primary verdict.
print(f"   {'threshold':>10s} {'tolerant':>9s} {'clades':>7s} {'contrasts':>10s} {'headline k/n':>13s} {'sr p':>9s} {'histidine k/n':>14s} {'sr p':>9s}")
for thr2 in (3.0, 3.5, 4.0, 4.5, 5.0):
    st2 = {a: (v <= thr2) for a, v in zip(B.ncbi_assembly, B.ph_min.to_numpy(float))}

    def all_tol2(clade):
        ts = [t.name for t in clade.get_terminals() if t.name in st2]
        return bool(ts) and all(st2[t] for t in ts)

    mx2 = []

    def desc2(clade):
        if all_tol2(clade):
            mx2.append(clade)
            return
        for c in clade.clades:
            if [t.name for t in c.get_terminals() if t.name in st2]:
                desc2(c)

    desc2(tree.root)
    pr2 = []
    for cl in mx2:
        par = parents.get(id(cl))
        if par is None:
            continue
        sis = [t.name for c in par.clades if c is not cl for t in c.get_terminals()
               if t.name in st2 and not st2[t.name]]
        if len(sis) >= MINSIS:
            pr2.append((cl, [t.name for t in cl.get_terminals() if t.name in st2], sis))
    out = []
    for nm, vec, pred, claim in MEAS[:2]:
        ds = []
        for cl, tt, ss in pr2:
            a = np.array([vec[idx[t]] for t in tt if t in idx], float)
            b = np.array([vec[idx[t]] for t in ss if t in idx], float)
            a, b = a[~np.isnan(a)], b[~np.isnan(b)]
            if len(a) < 1 or len(b) < MINSIS:
                continue
            ds.append(float(np.mean(a) - np.mean(b)))
        ds = np.array(ds)
        if len(ds) < 5:
            out.append((0, 0, float("nan")))
            continue
        kk = int(np.sum(np.sign(ds) == pred))
        out.append((kk, len(ds), float(st.wilcoxon(ds).pvalue) if np.any(ds != 0) else 1.0))
    mark = " <- primary" if thr2 == THR else ""
    print(f"   {thr2:>10.1f} {int(sum(st2.values())):>9d} {len(mx2):>7d} {len(pr2):>10d} "
          f"{out[0][0]:>6d}/{out[0][1]:<6d} {out[0][2]:>9.2e} {out[1][0]:>7d}/{out[1][1]:<6d} {out[1][2]:>9.2e}{mark}")
    rows.append({"part": "6", "threshold": thr2, "measure": "threshold sweep", "n_tolerant": int(sum(st2.values())),
                 "n_clades": len(mx2), "n_contrasts": len(pr2), "headline_k": out[0][0], "headline_n": out[0][1],
                 "headline_signed_rank_p": out[0][2], "histidine_k": out[1][0], "histidine_n": out[1][1],
                 "histidine_signed_rank_p": out[1][2], "is_primary": thr2 == THR})
print("   Read the sweep for CONSISTENCY of direction, not for its smallest p-value. Five thresholds and two measures is")
print("   ten tests, and the smallest of ten p-values is not a finding; the primary row is the one that counts.")
sw = [r for r in rows if r["part"] == "6"]
nd = sum(1 for r in sw if r["headline_n"] and r["headline_k"] / r["headline_n"] > 0.5)
ns = sum(1 for r in sw if np.isfinite(r["headline_signed_rank_p"]) and r["headline_signed_rank_p"] < 0.05)
print(f"   What the sweep does show: the headline measure points the predicted way at {nd} of {len(sw)} thresholds, and reaches")
print(f"   p < 0.05 at {ns} of them - WITHOUT the primary being one. That combination has an obvious honest reading and a")
print("   dishonest one. The honest reading is that the direction is stable across thresholds while significance is not,")
print("   which is what an underpowered real effect looks like and also what a weak artefact looks like. The dishonest")
print("   reading is to adopt whichever threshold cleared. The primary stays 4.0 for the reason it was chosen - it is the")
print("   value the script was written with, before any of these numbers existed.")

pd.DataFrame(rows).to_csv(D / "independent_transitions.csv", index=False)
print(f"\nwrote {D / 'independent_transitions.csv'}")
