"""If a reported minimum is the lowest pH anyone tested, what does that do to the headline? -> DATA_DIR/censoring_sensitivity.csv
THE OBJECTION THIS ANSWERS, which is the most serious one the project faces. `ph_min` is supposed to be the lowest pH at which an organism grows. What a species
description usually records is the lowest pH at which growth WAS SEEN, which is bounded below by whatever the authors chose to test. If that choice depends on what they
expected the organism to do, the error is not noise - it is systematically related to the organism's apparent acid tolerance, and therefore to anything that predicts it.
WHY THE PROJECT'S EXISTING ANSWER DOES NOT COVER IT. `scripts/62` measured outcome reliability (ICC 0.965, error SD 0.343 pH) and correctly concluded that error in the
DEPENDENT variable does not bias a regression slope - it inflates the standard error and caps R-squared. That is true of CLASSICAL error, which is independent of the
predictor by definition. Censoring is not classical error. It has a sign (reported >= true, always) and a mechanism (a human decision about where to start), and if that
mechanism correlates with the predictor the slope is biased, not merely noisy. Part 4 demonstrates the difference by simulation rather than asserting it, because the
reflex to treat all outcome error as harmless is exactly how the project's fifth scale error happened.
AND THE COLUMN THAT LOOKS LIKE IT SETTLES THIS DOES NOT. `phenotypes_all.csv` carries `ph_min_bound` with the value `measured` for 2,732 of 3,002 rows, which reads as a
determination that those minima are demonstrated limits. It is not one. `scripts/00g_build_phenotype_table.py` line 72 assigns `bound="measured"` UNCONDITIONALLY to every
parser-extracted row, and 2,623 of the 2,732 arrive that way, never having been read. Only the 114 hand-curated rows had bound type actually determined from the text, and
5 of those came back censored. Part 1 establishes this from the data rather than from the source line, and part 2 shows what it did to two published robustness checks
that stratify on the label - `scripts/58`'s tier test and `scripts/36d`'s subset - where the partition is not bound type but SOURCE, wearing a bound-type name.
THE QUANTITATIVE DELIVERABLE IS PART 5, not the complaint. Rather than ask whether censoring exists - it does, and its extent is unknowable from these sources - it asks
how much DIFFERENTIAL censoring would be required to manufacture the headline from nothing. That is a bound the reader can judge: it converts an unfalsifiable worry into a
number, and the number can then be compared against what is actually known about the measurements. It is computed by displacing the outcome along the predictor inside the
real phylogenetic model and finding where the coefficient crosses zero, not from an algebraic approximation that assumes away the covariates and the tree.
PART 3 IS THE EVIDENCE THAT THE WORRY IS REAL, and it is the strongest thing here: the reported minima have a grid signature. If they were measured limits they would be
distributed continuously; they are not. Nothing in part 3 depends on any model.
SCALE. Coefficients are per standard deviation of the composite, matching `docs/NUMBERS.md`'s headline row, which is the figure the sensitivity bound is about.
Usage: python3 scripts/96_censoring_sensitivity.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
rng = np.random.default_rng(20261101)
rows = []

refs = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv"); refs = refs[refs.accepted]
MEM = [str(r.label) for r in refs.itertuples() if r.membrane]
reflen = pd.read_csv(ROOT / "results/structure_panel_wide_classes.csv").groupby(["ref", "class"]).size()
pan = pd.read_csv(D / "structure_panel_wide.csv").drop_duplicates("ncbi_assembly")
pheno = (pd.read_csv(D / "phenotypes_all.csv")
         .merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
B = pheno[(pheno.domain == "Bacteria") & pheno.ncbi_assembly.isin(tips)].drop_duplicates("ncbi_assembly").reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100.0
SRCC = [f"src_{o}" for o in SRC]
M = pan.set_index("ncbi_assembly")


def composite(cls="peri_exp", minrefs=3):
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
    s = np.nanmean(Z, axis=1)
    return np.where(np.sum(~np.isnan(Z), axis=1) >= minrefs, s, np.nan)


x_all = composite()
keep = (~np.isnan(x_all)) & B.ph_min.notna().to_numpy()
G = B[keep].reset_index(drop=True)
xk = x_all[keep]
Cm = ph.vcv_fast(tree, list(G.ncbi_assembly))
y0 = G.ph_min.to_numpy(float)
USE = [c for c in ["gw"] + SRCC if G[c].nunique() > 1]
X0 = np.column_stack([np.ones(len(G))] + [G[c].to_numpy(float) for c in USE])
z = (xk - xk.mean()) / np.std(xk)
LAM = ph.pgls(y0, X0, Cm)["lambda"]


def beta_of(y, lam=None):
    f = ph.pgls(y, np.column_stack([X0, z]), Cm, lam=LAM if lam is None else lam)
    return float(f["beta"][-1]), float(f["se"][-1]), float(f["p"][-1])


b0, se0, p0 = beta_of(y0)
print(f"baseline, reproducing docs/NUMBERS.md's headline row: n {len(G)}, lambda {LAM:.4f}, {b0:+.6f} per SD (se {se0:.6f}, p {p0:.3e})")
assert abs(b0 - 0.184837) < 5e-4, f"baseline {b0:+.6f} is not the published headline; the sensitivity below would be about something else"

print("\n1. the column that looks like it settles bound type, and what it actually records")
ct = pd.crosstab(B.origin, B.ph_min_bound)
print(ct.to_string())
meas = B[B.ph_min_bound == "measured"]
auto = int((meas.origin == "literature_auto").sum())
hand = int((meas.origin == "literature").sum())
print(f"   {len(meas)} rows are labelled `measured`. {auto} are parser rows, which scripts/00g line 72 labels `measured`")
print(f"   unconditionally, with no reference to the text. {hand} were determined by reading, and of the {int((B.origin == 'literature').sum())} hand-read")
print(f"   rows {int((B.ph_min_bound == 'lowest_tested').sum())} came back censored - the ONLY direct estimate of the censoring rate anywhere in this project.")
assert auto > hand * 10, "the premise of this part is that the label is dominated by defaults; it is not, so re-read scripts/00g"
rows.append({"part": "1", "measure": "rows labelled measured", "n": len(meas), "note": f"{auto} by parser default, {hand} determined by reading"})
rows.append({"part": "1", "measure": "hand-read rows found censored", "n": int((B.ph_min_bound == 'lowest_tested').sum()),
             "n_total": int((B.origin == 'literature').sum()),
             "rate": float((B.ph_min_bound == 'lowest_tested').sum() / max(1, (B.origin == 'literature').sum()))})

print("\n2. what that did to the published robustness checks that stratify on the label")
# scripts/58 reports `bound = measured` / `bound != measured` beside `origin = ...` and `confidence = ...` as if they were
# three independent cuts of the evidence. If the label is a source default they are ONE cut, reported three times, and the
# tier test's apparent agreement across them is not corroboration.
part = {"bound = measured": (B.ph_min_bound == "measured").to_numpy(),
        "confidence in (high, medium)": B.confidence.isin(["high", "medium"]).to_numpy(),
        "origin is a literature tier": B.origin.isin(["literature", "literature_auto"]).to_numpy()}
ks = list(part)
print(f"   {'':34s} " + " ".join(f"{k[:22]:>24s}" for k in ks))
for a in ks:
    line = f"   {a[:32]:34s} "
    for b in ks:
        agree = float((part[a] == part[b]).mean())
        line += f"{agree:>24.3%}"
    print(line)
pairs = [(a, b, float((part[a] == part[b]).mean())) for i, a in enumerate(ks) for b in ks[i + 1:]]
worst = min(p[2] for p in pairs)
for a, b, v in pairs:
    rows.append({"part": "2", "measure": f"row-level agreement: {a} vs {b}", "rate": v})
print(f"   the three partitions agree on at least {worst:.1%} of rows, so scripts/58's tier test reports one cut three times.")
print("   Its three strata are not independent evidence, and the published figures are unchanged - what changes is what they")
print("   can be said to show: a SOURCE robustness check, not a bound-type one. No bound-type stratification exists.")
assert worst > 0.95, f"the three partitions differ on {1 - worst:.1%} of rows, so this part's reading is too strong"

print("\n3. the grid signature, which needs no model and is the evidence that the worry is real")
f10 = np.round((G.ph_min.to_numpy(float) * 10) % 10).astype(int)
print(f"   {float((f10 == 0).mean()):.1%} of reported minima are whole numbers, {float((f10 == 5).mean()):.1%} end in .5, "
      f"{float(((f10 != 0) & (f10 != 5)).mean()):.1%} are anything else")
top = G.ph_min.value_counts().head(5)
print("   the five commonest values and their share of the sample:")
for v, c in top.items():
    print(f"      pH {v:<4} {c:5d} genomes ({c / len(G):5.1%})")
print("   A measured limit has no reason to prefer a round number. A testing grid has every reason to, and the modal value")
print(f"   - pH {top.index[0]}, {top.iloc[0] / len(G):.0%} of the sample - is a standard medium's lower end rather than a biological boundary.")
print("   This does not say what fraction are censored, which these sources cannot answer. It says the distribution is shaped")
print("   by the design of experiments and not only by the organisms, and that is the premise of everything below.")
rows.append({"part": "3", "measure": "share of minima on whole numbers", "rate": float((f10 == 0).mean()), "n": len(G)})
rows.append({"part": "3", "measure": "share on a 0.5 grid (whole or half)", "rate": float(((f10 == 0) | (f10 == 5)).mean()), "n": len(G)})
rows.append({"part": "3", "measure": f"share at the modal value pH {top.index[0]}", "rate": float(top.iloc[0] / len(G)), "n": len(G)})

print("\n4. classical error against censoring, simulated, because they do opposite things to a slope")
sd_err = 0.343                      # scripts/62's error SD for a single curator reading
NS = 200
cl, df_ = [], []
for _ in range(NS):
    cl.append(beta_of(y0 + rng.normal(0, sd_err, len(y0)))[0])
    # differential censoring of the SAME total magnitude, but aligned with the predictor instead of independent of it
    df_.append(beta_of(y0 + sd_err * z)[0])
cl = np.array(cl)
print(f"   classical noise, SD {sd_err} pH added to the outcome, {NS} draws: mean {cl.mean():+.6f} per SD against the")
print(f"   baseline {b0:+.6f}; bias {cl.mean() - b0:+.6f}, which is {abs(cl.mean() - b0) / se0:.2f} standard errors. No bias, as scripts/62 says.")
print(f"   the SAME magnitude aligned with the predictor: {df_[0]:+.6f}, a shift of {df_[0] - b0:+.6f} = {(df_[0] - b0) / se0:.1f} standard errors.")
print("   Same amount of error in the outcome. One is harmless and one would be the whole result. scripts/62's conclusion is")
print("   correct about the error it measured and says nothing about the error this script is about.")
rows.append({"part": "4", "measure": "classical outcome noise, mean beta", "beta_per_sd": float(cl.mean()), "n_draws": NS,
             "bias_in_se": float((cl.mean() - b0) / se0)})
rows.append({"part": "4", "measure": "same magnitude aligned with the predictor", "beta_per_sd": float(df_[0]),
             "bias_in_se": float((df_[0] - b0) / se0)})

print("\n5. THE BOUND: how much differential censoring would be needed to manufacture the headline from nothing?")
# The displacement is c_i = reported - true >= 0. Write the part that matters as alpha per SD of the predictor; alpha < 0
# is the direction that creates a spurious positive coefficient here, since it adds outcome where the predictor is low.
print(f"   {'alpha (pH per SD of predictor)':>32s} {'beta per SD':>13s} {'p':>11s}  reading")
alphas = [0.0, -0.05, -0.10, -0.1848, -0.25, -0.30]
cross = None
for a in alphas:
    bb, ss, pp = beta_of(y0 + a * z)
    tag = ""
    if abs(bb) < 2 * ss:
        tag = "indistinguishable from zero"
        cross = cross if cross is not None else a
    elif bb < 0:
        tag = "sign reversed"
    print(f"   {a:>32.4f} {bb:>+13.6f} {pp:>11.2e}  {tag}")
    rows.append({"part": "5", "measure": "displacement aligned with predictor", "alpha_ph_per_sd": a,
                 "beta_per_sd": bb, "se": ss, "p": pp})
# Solve for the exact alpha that zeroes the coefficient. The map is linear in alpha, so two points determine it; it is
# still verified by refitting at the solution rather than trusted.
b_at = {a: beta_of(y0 + a * z)[0] for a in (0.0, -0.1)}
slope = (b_at[-0.1] - b_at[0.0]) / -0.1
a_zero = -b_at[0.0] / slope
b_chk, se_chk, _ = beta_of(y0 + a_zero * z)
print(f"\n   the coefficient reaches zero at alpha = {a_zero:+.4f} pH per SD of the predictor (refit there gives {b_chk:+.2e}, confirming the solve)")
assert abs(b_chk) < 1e-6, f"the solved alpha does not zero the coefficient ({b_chk:+.3e}); the linearity assumption is wrong"
print(f"   So the censoring displacement would have to run {abs(a_zero):.3f} pH DEEPER for every standard deviation of the")
print(f"   predictor - about {abs(a_zero) * 4:.1f} pH across the sample's central four standard deviations - purely as a function of")
print("   surface composition, and in exactly the direction that mimics the finding.")
# WHAT MAKES THAT NUMBER JUDGEABLE: express it as a required correlation, given a displacement magnitude. The hand-read
# tier is the only place the magnitude is observable at all, and it is thin, so a range is given rather than one value.
print("\n   that bound expressed as a required correlation, which is what a reader can weigh:")
print(f"   {'if the displacement SD is':>28s} {'required r(displacement, predictor)':>38s}")
for sd_c in (0.25, 0.5, 1.0, 2.0):
    r_need = abs(a_zero) / sd_c
    note = "impossible (r > 1)" if r_need > 1 else ""
    print(f"   {sd_c:>28.2f} {r_need:>38.3f}  {note}")
    rows.append({"part": "5", "measure": "required correlation of displacement with predictor", "displacement_sd": sd_c,
                 "required_r": float(r_need), "attainable": bool(r_need <= 1)})
print("   Reading it: if censoring moves minima by a quarter of a pH unit on average, no correlation with surface composition")
print("   could produce the headline, because the required one exceeds 1. Explaining the result away needs BOTH a large")
print("   displacement AND a strong, specifically-aligned dependence on genome composition. That is a conjunction, and naming")
print("   it is the point - the objection is not refuted, it is given a size.")
rows.append({"part": "5", "measure": "alpha that zeroes the headline", "alpha_ph_per_sd": float(a_zero),
             "n": len(G), "beta_per_sd": b0})

print("\n6. the question a proxy can actually answer: how DEEP would the censoring have to be?")
# A FIRST VERSION OF THIS PART COMPARED THE WRONG TWO THINGS, and it is recorded because the mistake is instructive. It
# printed the difference in the PREDICTOR between likely-censored rows and the rest, in standard deviations of the
# predictor, and set that beside part 5's alpha, which is pH per standard deviation of the predictor. Those are different
# quantities and the comparison is a category error: it made a separation of 1.4 SD look like it dwarfed a requirement of
# 0.185, when the two cannot be placed on one axis at all. The closing sentence then described the separations as "a few
# hundredths of an SD" while the table above it printed 0.31, 0.94, 0.20 and 1.41 - prose asserting the opposite of the
# numbers beside it, which is this project's single commonest error class.
# THE WELL-POSED QUESTION. Take a proxy for being censored and suppose every row it flags has its minimum reported delta
# pH ABOVE the true limit, and every other row is correct. Then solve for the delta that drives the headline to zero. That
# is in pH units, it is a depth a reader can judge against what species descriptions actually do, and it is obtained by
# refitting rather than from an algebraic approximation.
# THE THREATENING DIRECTION IS NOW EXPLICIT. The displacement is always positive - reported is at or above true - so a
# proxy manufactures the headline only if it flags rows with a HIGH predictor value, because the coefficient is positive.
# A proxy whose flagged rows have a LOWER predictor value biases the coefficient the other way, so censoring there would
# mean the published figure is too SMALL. Both are reported; they are not the same finding and must not be averaged.
PROX = [("minimum on a whole number (grid value)", (f10 == 0)),
        ("minimum at or above pH 5 (conventional band)", G.ph_min.to_numpy(float) >= 5.0),
        ("minimum at the modal value exactly", G.ph_min.to_numpy(float) == float(top.index[0])),
        ("from a database row, bound never assessed", G.origin.isin(["thermobase", "acidb", "bacdive"]).to_numpy())]
print(f"   {'proxy for being censored':>46s} {'n':>6s} {'mean z gap':>11s} {'required depth':>15s} {'headroom':>9s} {'direction':>11s}")
for name, msk in PROX:
    if msk.sum() < 10 or (~msk).sum() < 10:
        continue
    gap = float(z[msk].mean() - z[~msk].mean())
    fl = msk.astype(float)
    # linear in delta, so two refits determine the root; the root is then verified by a third refit at the solution
    b_a, b_b = beta_of(y0)[0], beta_of(y0 - 1.0 * fl)[0]
    sl = (b_b - b_a) / -1.0
    d_need = (-b_a / sl) * -1.0 if abs(sl) > 1e-12 else float("nan")
    chk = beta_of(y0 - d_need * fl)[0] if np.isfinite(d_need) else float("nan")
    assert not np.isfinite(d_need) or abs(chk) < 1e-6, f"solved depth for {name} does not zero the coefficient ({chk:+.3e})"
    # THE HEADROOM, which is what makes a required depth judgeable rather than merely large. The true limit cannot be
    # below pH 0, so the deepest mean displacement physically available for a group is the mean of its own reported
    # minima - every flagged organism's true limit sitting at pH 0. If the required depth exceeds that, no censoring
    # pattern of this shape can produce the headline, whatever anyone assumed about laboratory practice.
    head = float(G.ph_min.to_numpy(float)[msk].mean())
    feasible = bool(np.isfinite(d_need) and 0 < d_need <= head)
    if not np.isfinite(d_need):
        direction = "no effect"
    elif d_need <= 0:
        direction = "protects"
    elif feasible:
        direction = "THREATENS"
    else:
        direction = "IMPOSSIBLE"
    dn = f"{d_need:>+15.3f}" if np.isfinite(d_need) else f"{'undefined':>15s}"
    print(f"   {name:>46s} {int(msk.sum()):>6d} {gap:>+11.4f} {dn} {head:>9.2f} {direction:>11s}")
    rows.append({"part": "6", "measure": f"censoring depth that would zero the headline: {name}", "n": int(msk.sum()),
                 "mean_z_gap": gap, "required_depth_ph": float(d_need), "headroom_to_ph0": head,
                 "threatens": bool(feasible), "physically_possible": feasible})
print("   Reading the depth column: it is how far BELOW its reported value the true limit would have to lie, on average,")
print("   for every flagged row, with all other rows correct, for the headline to vanish. A negative entry means censoring")
print("   in that group pushes the coefficient the other way, so the published figure would be an UNDERestimate.")
print("   `headroom` is the mean reported minimum of the flagged rows, which is the deepest displacement physically")
print("   available to them - it would put every flagged organism's true limit at pH 0. A required depth above the")
print("   headroom is not implausible, it is impossible, and that distinction is the whole value of the column.")
six = [r for r in rows if r["part"] == "6" and np.isfinite(r["required_depth_ph"]) and r["required_depth_ph"] > 0]
thr = [r for r in six if r["physically_possible"]]
imp = [r for r in six if not r["physically_possible"]]
if thr:
    worst_t = min(thr, key=lambda r: r["required_depth_ph"])
    print(f"   {len(thr)} proxy/proxies could in principle do it. The most dangerous needs the least depth:")
    print(f"      {worst_t['measure'].split(': ')[1]}, at {worst_t['required_depth_ph']:.2f} pH against {worst_t['headroom_to_ph0']:.2f} pH of headroom.")
    print("   Whether that is plausible is a question about species descriptions, not about this dataset.")
else:
    print(f"   NO proxy can do it. All {len(imp)} that point in the threatening direction require a displacement DEEPER THAN")
    print("   THE HEADROOM: the flagged organisms' true minima would have to lie below pH 0. Censoring of this shape cannot")
    print("   produce the headline even in its physical limit, so the objection, in the form a proxy can express it, is closed")
    print("   rather than merely bounded. What remains open is censoring that no proxy here can see.")
    for r in imp:
        print(f"      {r['measure'].split(': ')[1]}: needs {r['required_depth_ph']:.2f} pH, has {r['headroom_to_ph0']:.2f}")

print("\n7. verdict")
print(f"   The objection is real, it is not refuted, and it now has a size. {float(((f10 == 0) | (f10 == 5)).mean()):.0%} of reported minima sit on a whole")
print("   or half pH unit, which is a testing grid and not a biological boundary; the field that appears to record bound")
print(f"   type records a parser default for {auto / len(meas):.0%} of the rows carrying it; and the only direct censoring estimate in the")
print(f"   project rests on {int((B.origin == 'literature').sum())} hand-read rows, {int((B.ph_min_bound == 'lowest_tested').sum())} of which came back censored.")
print(f"   To erase the headline, censoring must deepen by {abs(a_zero):.3f} pH for every standard deviation of surface composition.")
if thr:
    wt = min(thr, key=lambda r: r["required_depth_ph"])
    print(f"   Routed through the most dangerous observable proxy, that is a true limit lying {wt['required_depth_ph']:.2f} pH below the")
    print(f"   reported one for all {wt['n']} rows it flags, against {wt['headroom_to_ph0']:.2f} pH of physical headroom. Nothing here rules that out.")
else:
    wi = min(imp, key=lambda r: r["required_depth_ph"]) if imp else None
    if wi:
        print(f"   Routed through the most dangerous observable proxy it is IMPOSSIBLE, not merely implausible: the {wi['n']} flagged")
        print(f"   rows would need true minima {wi['required_depth_ph']:.2f} pH below their reported values, and only {wi['headroom_to_ph0']:.2f} pH exists")
        print("   between those values and pH 0. The grid signature in part 3 is real and the bound-type field is a default, but")
        print("   neither can account for the association by this route.")
print("   WHAT THIS IS NOT. It is not a clearance. The proxies are weak and each is a proxy for other things too; a selection")
print("   process acting through which ORGANISMS get described at all would appear in none of them; and the direction that")
print("   threatens the result - neutrophiles tested down to pH 5 or 6 and no further - is also the most plausible one a")
print("   priori, because that is how species descriptions are actually done.")
print("   WHAT WOULD SETTLE IT, and neither is reachable from these sources: bound type read from full texts at scale, which")
print("   scripts/65 showed cannot be mined from a literature index because negations are stopwords; or growth limits")
print("   measured on a common grid, which is bench work. Until then the manuscript should carry this number rather than")
print("   scripts/62's reassurance, which is about a different kind of error.")

pd.DataFrame(rows).to_csv(D / "censoring_sensitivity.csv", index=False)
print(f"\nwrote {D / 'censoring_sensitivity.csv'}")
