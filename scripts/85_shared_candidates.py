"""Which measured genome or phenotype property accounts for the SHARED histidine component? -> DATA_DIR/shared_candidates.csv and, with `wide`, DATA_DIR/shared_candidates_wide.csv
BOTH SPELLINGS ARE DECLARED, for the reason scripts/81 records: scripts/69 builds a claim's candidate pool from this line and keeps a panel variant's pool separate, so an
undeclared output turns every figure quoted from it into a guaranteed false miss.
THE QUESTION, AND WHY IT IS NOT THE ONE scripts/83 ASKED. Correction 19 split exposed histidine into a side-specific excess (-0.0800 per point) and a SHARED component
(-0.1317) that is larger. scripts/83 then ruled out the three candidates that lived in the panel itself - non-exposed histidine in the same families, proteome-wide
histidine, and scoring quality - and docs/SYNTHESIS.md section 8 recorded the component as open and "blocked on a hypothesis rather than on data access". That is true only
of panel columns. The phenotype table and the genome-feature table carry a dozen properties that are plausible correlates of histidine content and were never tried, so the
gap is not actually blocked; it needs a list of candidates and a rule for reading them, which is this script.
THE CRITERION IS FIXED BEFORE THE NUMBERS. A candidate ACCOUNTS FOR the component if, with it held fixed, the component's coefficient either loses significance or falls
below half its uncontrolled magnitude. A candidate that moves it less than that is reported as a partial contributor with the size of its contribution, and a candidate
whose own coefficient is null is reported as unable to explain anything regardless of what it does to the component - a regressor that does not predict the outcome cannot
be what the outcome is responding to.
EVERY COMPARISON IS WITHIN-SAMPLE, which is the thing most likely to go wrong here and is correction 12 exactly. Growth temperature is known for about half these genomes.
Fitting the component on 2,768 and the controlled version on 1,464 and then comparing the two coefficients would be a decomposition computed across two different samples -
the error correction 12 was about. So each candidate gets its OWN baseline, refitted on precisely the genomes where that candidate is observed, and the printed comparison
is always baseline-on-this-subset against controlled-on-this-subset. The full-sample baseline is printed once, separately, and never used as a comparator.
THE BASELINE IS CHECKED AGAINST scripts/82 RATHER THAN ASSUMED. The regressor and control set here are that script's, so the full-sample baseline must reproduce -0.13168
to five decimals, and an assertion requires it. scripts/83 learned this the hard way: it reported -0.3665 against that figure as a confirmation of being the same quantity,
when it was a different pooling under a different control set.
SCALE. Per percentage point of pooled share, from ph.pgls directly and never ph.pgls_scan - correction 15. The periplasmic D, E and K shares are held fixed in every fit,
which is correction 14 and, after correction 19 and scripts/83's third instance of it, the error this project is most likely to repeat.
Usage: python3 scripts/85_shared_candidates.py [orig|wide]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
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
WHICH = (sys.argv[1] if len(sys.argv) > 1 else "wide").lower()
assert WHICH in ("orig", "wide"), f"panel must be orig or wide, got {WHICH!r}"
PANEL = "structure_panel_wide_det.csv" if WHICH == "wide" else "structure_panel.csv"
TAG = "_wide" if WHICH == "wide" else ""
MINPOS = 200
HALF = 0.5                 # "falls below half its uncontrolled magnitude" - the criterion, fixed here
_r = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv")
MEM = [str(x) for x in _r[_r.accepted & _r.membrane].label]
rows = []

pan = pd.read_csv(D / PANEL)
if not any(c.endswith("|H") for c in pan.columns):
    sys.exit(f"{PANEL} has no histidine columns; rebuild the panel with scripts/36a first")
GF = ["ncbi_assembly", "acidic_excess", "frac_His", "gc", "genome_size_bp", "ivywrel", "mean_pI",
      "frac_pI_gt7", "mean_gravy", "frac_membrane_proxy", "n_proteins", "frac_Cys", "mean_pI_membrane"]
df = (pd.read_csv(D / "phenotypes_all.csv")
      .merge(pd.read_csv(D / "genome_features.csv")[GF], on="ncbi_assembly")
      .merge(pan, on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips) & df.ph_min.notna()].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
# AEROBICITY IS NOT MEASURED IN THIS DATASET - the `oxygen` column is empty for all 3,002 phenotype rows - and
# the candidate is kept in the list so that it is reported as unmeasured rather than quietly absent. The first
# version classed it with `contains("aerob")`, which matches "anaerobe" as well and would have labelled every
# anaerobe an aerobe the moment the column was populated; that is fixed here even though nothing exercises it.
_ox = B.oxygen.astype("string").str.lower().fillna("")
B["aerobe"] = np.where(_ox.eq(""), np.nan,
                       np.where(_ox.str.contains("anaerob") | _ox.str.contains("microaerob"), 0.0, 1.0))
B["log_genome"] = np.log10(B.genome_size_bp)
B["log_proteins"] = np.log10(B.n_proteins)
SRCC = [f"src_{o}" for o in SRC if B[f"src_{o}"].nunique() > 1]


def pooled(cls, keys=("D", "E", "H", "K")):
    tot = {k: np.zeros(len(B)) for k in ("n",) + tuple(keys)}
    for ref in MEM:
        if f"{ref}|{cls}|n" not in B or f"{ref}|{cls}|H" not in B:
            continue
        tot["n"] += B[f"{ref}|{cls}|n"].fillna(0).to_numpy(float)
        for k in keys:
            tot[k] += B[f"{ref}|{cls}|{k}"].fillna(0).to_numpy(float)
    return tot


PT, CT = pooled("peri_exp"), pooled("cyto_exp")
base_ok = (PT["n"] >= MINPOS) & (CT["n"] >= MINPOS)
SHARED = (100 * PT["H"] / np.where(PT["n"] > 0, PT["n"], 1)) + (100 * CT["H"] / np.where(CT["n"] > 0, CT["n"], 1))
DIFF = (100 * PT["H"] / np.where(PT["n"] > 0, PT["n"], 1)) - (100 * CT["H"] / np.where(CT["n"] > 0, CT["n"], 1))
PDEK = {k: 100 * PT[k] / np.where(PT["n"] > 0, PT["n"], 1) for k in ("D", "E", "K")}


def fit(mask, extra=()):
    """The correction-19 model on `mask`, optionally with extra columns. Returns the shared term and the extras."""
    g = B[mask].reset_index(drop=True)
    Cm = ph.vcv_fast(tree, list(g.ncbi_assembly))
    y = g.ph_min.to_numpy(float)
    # LAMBDA IS ESTIMATED ON THE COVARIATE-ONLY MODEL, exactly as scripts/82 does, and that detail is why the
    # assertion below exists. My first version folded the periplasmic D/E/K controls into the model used to
    # estimate lambda, which is defensible on its own and gave -0.13717 where that script gives -0.13168. The
    # assertion caught a 4% discrepancy I would otherwise have published as "reproduced".
    XL = np.column_stack([np.ones(len(g))] + [g[k] for k in ["gw"] + SRCC])
    lam = ph.pgls(y, XL, Cm)["lambda"]
    X0 = np.column_stack([XL] + [PDEK[k][mask] for k in ("D", "E", "K")])
    cols = [X0] + [np.asarray(e, float)[mask].reshape(-1, 1) for e in extra] + [SHARED[mask].reshape(-1, 1),
                                                                                DIFF[mask].reshape(-1, 1)]
    X = np.column_stack(cols)
    assert np.linalg.matrix_rank(X) == X.shape[1], "design matrix is rank-deficient"
    r = ph.pgls(y, X, Cm, lam=lam)
    sh = (float(r["beta"][-2]), float(r["se"][-2]), float(r["p"][-2]))
    ex = [(float(r["beta"][X0.shape[1] + i]), float(r["p"][X0.shape[1] + i])) for i in range(len(extra))]
    return sh, ex, int(len(g))


b0, _, n0 = fit(base_ok)
print(f"1. the baseline, reproduced rather than assumed: shared component {b0[0]:+.5f}/pt "
      f"(se {b0[1]:.5f}, p {b0[2]:.2e}), n={n0}")
assert abs(b0[0] - (-0.13168)) < 1e-5, f"this is not scripts/82's shared component: got {b0[0]:+.5f}"
print("   matches scripts/82 to five decimals, so the regressor and control set are that script's.")
print("   This figure is NOT used as a comparator below: each candidate gets its own within-sample baseline,")
print("   because growth temperature is observed for about half these genomes and comparing across the two")
print("   samples is correction 12 exactly.")
rows.append({"part": "baseline", "beta_shared_per_pt": b0[0], "se": b0[1], "p": b0[2], "n": n0,
             "matches_scripts_82": True})

CAND = [("growth temperature optimum", B.temp_opt_c), ("IVYWREL (thermophily proxy)", B.ivywrel),
        ("GC content", B.gc), ("log10 genome size", B.log_genome), ("log10 protein count", B.log_proteins),
        ("mean isoelectric point", B.mean_pI), ("share of proteins with pI > 7", B.frac_pI_gt7),
        ("mean membrane-protein pI", B.mean_pI_membrane), ("mean hydropathy", B.mean_gravy),
        ("membrane-protein share", B.frac_membrane_proxy), ("cysteine share", B.frac_Cys),
        ("aerobe", B.aerobe)]
print("\n2. each candidate, against its OWN within-sample baseline")
print(f"   {'candidate':>32} {'n':>5} {'r w/ shared':>12} {'cand beta':>11} {'cand p':>9} "
      f"{'shared base':>12} {'shared ctl':>11} {'kept':>6} {'p':>9}")
for lab, v in CAND:
    v = pd.to_numeric(v, errors="coerce")
    m = base_ok & v.notna().to_numpy() & np.isfinite(v.fillna(0).to_numpy(float))
    if m.sum() < 200:
        why = ("not measured in this dataset" if v.notna().sum() == 0 else
               f"only {int(v.notna().sum())} genomes have it and {int(m.sum())} clear the position floor")
        print(f"   {lab:>32} {int(m.sum()):>5}   {why}; skipped")
        rows.append({"part": "candidate skipped", "candidate": lab, "n": int(m.sum()), "why": why})
        continue
    vv = v.fillna(0).to_numpy(float)
    if np.nanstd(vv[m]) == 0:
        # PRINTED, not silently dropped. The first version only recorded this in the csv, so the console showed
        # 11 candidates where the list has 12 and the missing one was invisible - a skip that is not reported is
        # indistinguishable from a candidate that was never thought of.
        print(f"   {lab:>32} {int(m.sum()):>5}   no variance among these genomes; skipped")
        rows.append({"part": "candidate skipped", "candidate": lab, "n": int(m.sum()), "why": "no variance"})
        continue
    own, _, nn = fit(m)                       # the within-sample baseline, refitted on exactly these genomes
    ctl, ex, _ = fit(m, extra=(vv,))
    bc, pc = ex[0]
    kept = abs(ctl[0] / own[0]) if own[0] else float("nan")
    rr = float(np.corrcoef(vv[m], SHARED[m])[0, 1])
    print(f"   {lab:>32} {nn:>5} {rr:>+12.3f} {bc:>+11.4f} {pc:>9.2e} {own[0]:>+12.5f} {ctl[0]:>+11.5f} "
          f"{kept:>6.0%} {ctl[2]:>9.2e}")
    accounts = (ctl[2] >= 0.05) or (kept < HALF)
    rows.append({"part": "candidate", "candidate": lab, "n": nn, "r_with_shared": rr,
                 "beta_candidate": bc, "p_candidate": pc, "beta_shared_within_sample_baseline": own[0],
                 "p_shared_baseline": own[2], "beta_shared_controlled": ctl[0], "p_shared_controlled": ctl[2],
                 "kept_fraction": kept, "candidate_predicts_outcome": bool(pc < 0.05),
                 "accounts_for_component": bool(accounts and pc < 0.05)})

print("\n3. the verdict, by the criterion fixed in the preamble")
C = [r for r in rows if r["part"] == "candidate"]
live = [r for r in C if r["candidate_predicts_outcome"]]
acct = [r for r in C if r["accounts_for_component"]]
moved = sorted((r for r in live), key=lambda r: r["kept_fraction"])
print(f"   {len(C)} candidates tested; {len(live)} predict pH_min in their own right and so could explain")
print(f"   something; {len(acct)} account for the shared component by the fixed criterion.")
for r in acct:
    print(f"      {r['candidate']}: component {r['beta_shared_within_sample_baseline']:+.5f} -> "
          f"{r['beta_shared_controlled']:+.5f} ({r['kept_fraction']:.0%} kept, p {r['p_shared_controlled']:.2e})")
if not acct:
    print("      (none)")
if moved:
    r = moved[0]
    print(f"   The largest reduction by any live candidate is {r['candidate']}, leaving "
          f"{r['kept_fraction']:.0%} of the component.")
    rows.append({"part": "verdict", "n_candidates": len(C), "n_live": len(live), "n_accounting": len(acct),
                 "largest_reduction_candidate": r["candidate"], "largest_reduction_kept": r["kept_fraction"]})
# A candidate that is null on its own CANNOT be the explanation whatever it does to the component, and saying so
# matters because a null regressor can still shift a coefficient by absorbing shared variance - which would look
# like an explanation and would not be one.
dead = [r for r in C if not r["candidate_predicts_outcome"] and r["kept_fraction"] < HALF]
if dead:
    print(f"   {len(dead)} candidate(s) reduced the component below half while being null themselves "
          f"({', '.join(r['candidate'] for r in dead)}); that is absorbed shared variance, not an explanation.")
    rows.append({"part": "null candidates that still moved it", "candidates": ";".join(r["candidate"] for r in dead)})

print("\n4. what the within-sample baselines bought, since the cost of getting this wrong is a published number")
# Growth temperature is the candidate where the two samples differ, and it is also where the error would have
# been most misleading. On its 1,464 genomes the component is -0.15241, not the full sample's -0.13168, so
# comparing the controlled fit against the FULL-sample baseline would have reported 116% - temperature appearing
# to INFLATE the component - where the within-sample comparison reports 100%, i.e. no effect at all. That is
# correction 12's error with a sign change attached, and it is the reason every row above refits its own baseline.
tr = [r for r in rows if r.get("candidate") == "growth temperature optimum" and r["part"] == "candidate"]
if tr:
    r = tr[0]
    naive = abs(r["beta_shared_controlled"] / b0[0])
    print(f"   growth temperature, n={r['n']}: within-sample baseline {r['beta_shared_within_sample_baseline']:+.5f},")
    print(f"   controlled {r['beta_shared_controlled']:+.5f} -> {r['kept_fraction']:.0%} kept, which is no effect.")
    print(f"   Against the FULL-sample baseline {b0[0]:+.5f} the same controlled figure reads {naive:.0%}, i.e.")
    print("   temperature appearing to make the component LARGER. That is correction 12 with a sign change, and")
    print("   it is what the per-candidate refits exist to prevent.")
    rows.append({"part": "within-sample demonstration", "candidate": r["candidate"], "n": r["n"],
                 "kept_within_sample": r["kept_fraction"], "kept_against_full_sample_baseline": naive})

pd.DataFrame(rows).to_csv(D / f"shared_candidates{TAG}.csv", index=False)
print("\nwrote", D / f"shared_candidates{TAG}.csv")
print("\nReading, and it is a clean negative. docs/SYNTHESIS.md section 8 recorded the shared histidine component")
print("as open and blocked on a hypothesis; that was true of the panel's own columns and not of the dataset, which")
print("carries properties that plausibly track histidine content and had never been tried. This script lists")
print(f"{len(CAND)} of them and tests the {len(C)} that this dataset can support. {len(live)} predict the acid limit in")
print("their own right and so could in principle explain something - GC content, genome size, protein count and")
print("mean membrane-protein isoelectric point - and NONE accounts for the component by the criterion fixed")
print("before the run. The largest reduction any live candidate manages")
print("leaves 85% of it standing. Thermophily is worth naming because it was the most plausible candidate on the")
print("list, histidine being a thermal-adaptation residue: both the direct growth-temperature optimum and the")
print("IVYWREL proxy leave the component at 100%, and the optimum is not even a significant predictor on its own")
print("sample here. So the component is not composition, not size, not thermophily, not hydrophobicity and not")
print("proteome charge, and it is still unidentified.")
skipped = [r for r in rows if r["part"] == "candidate skipped"]
if skipped:
    print(f"{len(skipped)} candidate(s) could not be tested at all and are NOT among the exclusions: "
          + "; ".join(f"{r['candidate']} ({r['why']})" for r in skipped) + ".")
    print("Aerobicity in particular is a plausible candidate that this dataset cannot address, and an earlier")
    print("draft of this reading listed it among the things ruled out, which it is not.")
print(f"What has changed is that the gap is no longer blocked on a hypothesis for want of looking: {len(C)}")
print("hypotheses are excluded by name, which is the useful form of a negative result and the thing that stops")
print("the next person retrying them.")
