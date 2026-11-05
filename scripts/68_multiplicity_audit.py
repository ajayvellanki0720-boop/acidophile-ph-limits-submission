"""How many hypotheses has this project tested, and does the headline survive correcting for them? -> DATA_DIR/multiplicity_audit.csv
This is the one criticism that no amount of per-test rigour answers, and the synthesis has never put a number on it. The project screened genome-wide
features, then amino-acid composition, then position classes within protein structures, then individual positions, and it reports a result from the
third of those. A reader is entitled to ask how large the search was and whether the surviving p-value means anything once that is accounted for.
The audit needs a unit, and the obvious one is wrong. There are several thousand rows carrying a p- or q-value across the project's result files (the script counts
them rather than quoting a figure, because an earlier draft of this docstring quoted a number from a different glob than the script used), but the overwhelming majority are
the SAME predictor refitted under a different model, subset, tree or panel - robustness checks, which do not multiply the error rate, because they are
deliberate re-examinations of one hypothesis rather than independent chances to find something. The unit used here is therefore a distinct
(feature or position class) tested against pH_min, counted once no matter how often it was refitted, enumerated from the screen files that record the
project's actual search:
  - the genome-wide feature and system screen, 18 features
  - the amino-acid composition screen, 20 residues
  - the structural position classes, which is where the result lives
Both a principled count and a deliberately absurd one are reported, because the claim should not depend on which a reader prefers.
THE MORE INTERESTING RESULT IS THE SECOND HALF. This project has three documented false leads - histidine, Trk/Ktr and the Kdp ATPase - each of which
looked respectable and was withdrawn. If multiplicity were the project's problem, correcting for it should have caught them. It would not have: Kdp
reached p 1.2e-08 in the same screen, which clears Bonferroni at any count considered here by four orders of magnitude, and it was still wrong. What
killed it was the within-family test and the evidence tier, not its marginal p-value. So this audit establishes two things at once: that the headline
survives correction comfortably, and that surviving correction is nearly worthless as evidence in this dataset. The binding constraints are the
structural tests, which is where the project's confidence should be located and where this document should direct a sceptical reader.
WHAT THIS DOES NOT BOUND. Choices that were never recorded as a test are not countable: the position-count floor, which covariates enter a model, how a
position class is defined, which panel is primary. Those are the forking paths that matter most and no audit of recorded p-values can see them. The
count below is a lower bound on the search, the tests are not independent of each other, and the honest reading is that this addresses one specific
form of multiplicity rather than the general concern.
Usage: python3 scripts/68_multiplicity_audit.py   (env DATA_DIR=results/broad2)"""
import glob
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
HEADLINE_P = None
rows = []

print("1. the size of the search, counted by distinct hypothesis rather than by refit")
screen = pd.read_csv(D / "broadened_screen_total_effect.csv")
feat = screen[["feature", "p"]].drop_duplicates("feature")
print(f"   genome-wide feature and system screen: {len(feat)} distinct features")
aa = pd.read_csv(ROOT / "results/broad/his_drivers.csv")
aa = aa[aa.domain == "Bacteria"][["aa", "raw_p"]].drop_duplicates("aa").rename(columns={"raw_p": "p"})
print(f"   amino-acid composition screen:         {len(aa)} distinct residues")
panel = pd.read_csv(D / "structure_panel_test.csv")
classes = sorted(panel["class"].dropna().unique())
print(f"   structural position classes:           {len(classes)} distinct classes {classes}")
# the class-level hypothesis is the composite over a panel, which is what the headline is; take each class's best recorded p as its hypothesis
cls_p = panel.groupby("class").p.min().rename("p").reset_index().rename(columns={"class": "feature"})
H = pd.concat([feat.assign(layer="genome-wide feature"),
               aa.rename(columns={"aa": "feature"}).assign(layer="amino-acid composition"),
               cls_p.assign(layer="structural position class")], ignore_index=True)
N = len(H)
print(f"   -> {N} distinct hypotheses in the recorded search")
allp, nfiles = [], 0
for f in sorted(glob.glob(str(ROOT / "results/**/*.csv"), recursive=True)):
    try:
        d = pd.read_csv(f)
    except Exception:
        continue
    hit = False
    for c in d.columns:
        if c.lower() in ("p", "pval", "p_value", "p_raw", "q", "qval", "q_bh"):
            allp += list(pd.to_numeric(d[c], errors="coerce").dropna())
            hit = True
    nfiles += int(hit)
N_absurd = len(allp)
print(f"   for comparison, every p- or q-value row in every result file ({nfiles} files), counting each refit as a")
print(f"   separate hypothesis, which is indefensible but is the most conservative count available: {N_absurd}")
rows.append({"part": "search size", "layer": "genome-wide feature", "n_hypotheses": len(feat)})
rows.append({"part": "search size", "layer": "amino-acid composition", "n_hypotheses": len(aa)})
rows.append({"part": "search size", "layer": "structural position class", "n_hypotheses": len(cls_p)})
rows.append({"part": "search size", "layer": "principled total", "n_hypotheses": N})
rows.append({"part": "search size", "layer": "every p-value row (absurdly conservative)", "n_hypotheses": N_absurd})

print("\n2. the headline against those counts")
hn = pd.read_csv(D / "headline_numbers.csv")
row = hn[hn.figure.str.contains("32-reference rule-selected panel", regex=False)]
HEADLINE_P = float(row.p.iloc[0])
HEADLINE_B = float(row.beta_per_sd.iloc[0])
print(f"   headline: {HEADLINE_B:+.4f} per SD, p = {HEADLINE_P:.2e}")
for label, n in (("principled count", N), ("every p-value row", N_absurd)):
    thr = 0.05 / n
    print(f"   Bonferroni at {label:18s} (n={n:5d}): threshold {thr:.2e}; headline clears it by "
          f"{np.log10(thr / HEADLINE_P):.1f} orders of magnitude")
    rows.append({"part": "headline", "basis": label, "n_hypotheses": n, "bonferroni_threshold": thr,
                 "headline_p": HEADLINE_P, "orders_of_margin": float(np.log10(thr / HEADLINE_P))})
# Benjamini-Hochberg over the principled set
pv = np.sort(H.p.dropna().to_numpy())
k = np.arange(1, len(pv) + 1)
passing = pv <= 0.05 * k / len(pv)
n_bh = int(np.max(np.where(passing)[0]) + 1) if passing.any() else 0
print(f"   Benjamini-Hochberg at 0.05 over the {len(pv)} recorded hypotheses: {n_bh} clear it")
rows.append({"part": "headline", "basis": "BH over principled set", "n_hypotheses": len(pv), "n_passing_bh": n_bh})

print("\n3. the part that matters: correction would not have caught this project's false leads")
FAILED = {"hm_rate_kdp": "Kdp ATPase", "frac_His": "histidine fraction", "hm_rate_trk_ktr": "Trk/Ktr potassium uptake"}
thr = 0.05 / N
caught = []
for key, name in FAILED.items():
    r = screen[screen.feature == key]
    if r.empty:
        print(f"   {name:26s} not in the screen file")
        continue
    p = float(r.p.iloc[0])
    verdict = "WOULD have survived Bonferroni" if p < thr else "would have been caught by Bonferroni"
    print(f"   {name:26s} screen p {p:.2e}  ({verdict} at n={N})")
    caught.append((name, p < thr))
    rows.append({"part": "false leads", "lead": name, "screen_p": p, "bonferroni_threshold": thr,
                 "survives_bonferroni": bool(p < thr)})
print(f"   Kdp cleared the threshold by {np.log10(thr / float(screen[screen.feature == 'hm_rate_kdp'].p.iloc[0])):.1f} orders of magnitude")
nc = sum(1 for _, surv in caught if not surv)
print(f"   and was withdrawn anyway. So correction would have caught {nc} of the {len(caught)} false leads - the two weak ones - and missed")
print("   the one that mattered. Multiplicity correction is not what distinguishes the surviving lead from the failed ones.")

print("\n4. what did distinguish them")
print("   The failed leads were killed by the within-family test and the evidence tier, both of which the exposed-periplasmic")
print("   composite passes. Those are the binding constraints and are where a sceptical reader should look:")
try:
    tier = pd.read_csv(D / "tier_test_wide_panel.csv")
    t = tier[(tier.composite == "all 32, unweighted") & (tier.subset.str.startswith("origin = literature ("))]
    if not t.empty:
        print(f"     hand-curated tier, all 32 references: {float(t.beta_per_sd.iloc[0]):+.4f} +/- {float(t.se.iloc[0]):.4f} (p {float(t.p.iloc[0]):.4f}, n={int(t.n.iloc[0])})")
        rows.append({"part": "binding constraint", "test": "hand-curated evidence tier", "beta_per_sd": float(t.beta_per_sd.iloc[0]),
                     "p": float(t.p.iloc[0]), "n": int(t.n.iloc[0])})
except Exception as e:
    print("     (tier file unavailable:", e, ")")
try:
    rb = pd.read_csv(D / "exposure_robustness.csv")
    wf = rb[rb.test.astype(str).str.contains("family", case=False, na=False)]
    for x in wf.head(3).itertuples():
        print(f"     {str(x.test)[:44]:44s} {x.beta_per_sd:+.4f} (p {x.p:.2e}, n={int(x.n)})")
except Exception:
    pass

R = pd.DataFrame(rows)
R.to_csv(D / "multiplicity_audit.csv", index=False)
print("\nwrote", D / "multiplicity_audit.csv")
print("\nReading. The headline clears Bonferroni by several orders of magnitude on any count a reader might prefer, so multiplicity is")
print("not this project's problem - and the same arithmetic shows why that reassurance is nearly worthless, because the Kdp lead cleared")
print("it too and was wrong. Cite the within-family test and the evidence tier instead; they are what the failed leads failed.")
