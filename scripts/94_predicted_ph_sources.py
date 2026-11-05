"""Are there datasets of PREDICTED ph_min, and can this project use them? -> DATA_DIR/predicted_ph_sources.csv
THE QUESTION, AS ASKED. TEMPURA and the other trait databases were checked in docs/DATABASES.md for MEASURED growth pH. The separate question is whether anyone publishes
PREDICTED ph_min for genomes, which would be attractive here because the measured sample is the binding constraint: 137 validated rows, 3 of them high confidence. The
answer is that predicted ph_min exists at scale and this project cannot use a single value of it as an outcome, for a reason that is not a quality objection but a
circularity: every published genome-based ph_min predictor predicts FROM the predictor this project tests.
WHY THIS SCRIPT MEASURES THE CIRCULARITY INSTEAD OF ASSERTING IT. GenomeSPOT (Barnum et al. 2024, doi 10.1101/2024.03.22.586313) predicts pH from amino acid frequencies
computed separately per predicted protein localization, with extracellular proteins identified by an HMM, and reports that its strongest residue-level pH correlate is "the
difference between glutamic acid frequency for extracellular and intracellular proteins" at rho 0.56. That is a description of this project's own measure: an
extracellular-minus-intracellular acidic composition contrast. A claim that strong should not rest on reading two feature lists and judging them similar, so part 2
constructs GenomeSPOT's stated feature IN THIS PROJECT'S DATA and measures how much of the headline it already contains. If the overlap is high, a regression of a
GenomeSPOT-derived ph_min on this project's predictor is guaranteed to come back positive whether or not the biology is real, and the exercise is worthless as evidence.
WHAT IS NOT CIRCULAR, AND IS WORTH MORE THAN THE PREDICTIONS WOULD HAVE BEEN. metaTraits (Podlesny et al., NAR 2026, doi 10.1093/nar/gkaf1241) integrates culture-based
trait records and ships MEASURED-ONLY bulk exports (*_summary_no_predictions.tsv.gz). Its culture sources for pH are BacDive, BV-BRC and JGI IMG/GOLD - and the latter two
are the two rows docs/DATABASES.md lists as "Not checked". So the search for predicted values turns up net-new MEASURED values instead. Part 3 counts them.
THE TRAP IN PART 3, which is this project's oldest one. metaTraits exports are aggregated at family, genus and species level; there is no per-strain bulk table. This
project's unit is a strain with an assembly, and a species-level pH attached to a strain genome is exactly the error that caused six Alicyclobacillus rows to be dropped
from BacDive. So part 3 counts what is AVAILABLE TO CURATE, never what is usable as-is, and part 4 measures the quality objections rather than listing them: how many rows
rest on a single observation, and how many carry values that a primary source would have to confirm before anyone believed them.
BOUND TYPE IS STILL ABSENT. docs/DATABASES.md's standing point is that no database records whether a minimum is a demonstrated limit or merely the lowest value tested.
metaTraits does not change that, so everything part 3 counts would enter at confidence `low` under this project's existing rule, and nothing here proposes relaxing it.
SCALE. Part 2's coefficients are per percentage point of pooled residue share, from ph.pgls directly, never ph.pgls_scan - correction 15.
Usage: python3 scripts/94_predicted_ph_sources.py   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2)"""
import os
import re
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
MT = ROOT / "data/raw/metatraits/ncbi_species_ph_minimum_no_predictions.tsv"
MINPOS = 200
rows = []
# data/raw/* is gitignored project-wide, so part 3's input is NOT in the repository and this script cannot be run from
# a fresh clone without re-fetching it. The provenance is recorded here rather than in a note somewhere else, because a
# derived count whose input cannot be regenerated is not a result:
#   curl -O https://www.bork.embl.de/~robbani/metatraits/ncbi_species_summary_no_predictions.tsv.gz
#   gzcat ncbi_species_summary_no_predictions.tsv.gz | awk -F'\t' 'NR==1 || $3=="pH minimum"' > <MT>
# The export is metaTraits v1 with culture-based records only; its database stamps read 2025-07-29. The "no_predictions"
# spelling is the whole point of the file and must not be dropped: the sibling *_summary_all.tsv.gz merges GenomeSPOT's
# predicted minima into the same columns, and part 2 is the demonstration that using those here would be circular.
if not MT.exists():
    sys.exit(f"missing {MT}\nre-fetch it with the two commands in this script's header comment (data/raw/* is gitignored)")

print("1. what each source actually is, and whether its ph_min can serve as an outcome here")
SOURCES = [
    # (name, has predicted ph_min, what the prediction is computed from, usable as an outcome, why)
    ("TEMPURA", False, "n/a - temperature only, 8,639 strains", False,
     "no pH of any kind; used here only to fill temp_opt_c"),
    ("GenomeSPOT", True, "per-localization amino acid frequencies, pI, R/(R+K); extracellular set by HMM", False,
     "CIRCULAR: predicts from this project's predictor (part 2 measures how much)"),
    ("metaTraits (prediction-based pH minimum)", True, "GenomeSPOT, the only model listed for this trait (756 training genomes)", False,
     "CIRCULAR: inherits GenomeSPOT's features wholesale"),
    ("MICROPHERRET / Traitar (via metaTraits)", False, "KO and PFAM matrices from eggNOG-mapper", False,
     "no pH trait; and gene content is the OTHER arm this project tests (scripts/92)"),
    ("Bac_Dive-AI (via metaTraits)", False, "undisclosed; contributes a discrete 'acidophilic' call, 0 training genomes reported", False,
     "no continuous ph_min; discrete call is not an outcome"),
    ("Madin et al. 2020 / Barberan et al. 2017", False, "n/a - measured, optimum only", False,
     "optimum only, no limit; already recorded in docs/DATABASES.md"),
    ("metaTraits (culture-based pH minimum)", False, "n/a - MEASURED, from BacDive + BV-BRC + JGI IMG/GOLD", True,
     "NOT a prediction and NOT circular; the real find, counted in part 3"),
]
print(f"   {'source':42s} {'pred?':6s} {'usable?':8s} note")
for name, pred, _, usable, why in SOURCES:
    print(f"   {name:42s} {('yes' if pred else 'no'):6s} {('yes' if usable else 'NO'):8s} {why}")
npred = sum(1 for _, p, _, _, _ in SOURCES if p)
nus = sum(1 for _, _, _, u, _ in SOURCES if u)
print(f"   {npred} of {len(SOURCES)} sources publish a predicted ph_min; {nus} source is usable as an outcome here, and it is not one of them")
assert nus == 1 and not any(p for _, p, _, u, _ in SOURCES if u), \
    "the inventory claims a predicted source is usable as an outcome; part 2 exists because that would need separate argument"

print("\n2. the circularity, measured rather than asserted: how much of the headline does GenomeSPOT's stated pH feature already contain?")
pan = pd.read_csv(D / "structure_panel_wide_det.csv")
_r = pd.read_csv(ROOT / "data/curated/structure_panel_refs.csv")
MEM = [str(x) for x in _r[_r.accepted & _r.membrane].label]
df = (pd.read_csv(D / "phenotypes_all.csv")
      .merge(pd.read_csv(D / "genome_features.csv")[["ncbi_assembly", "acidic_excess"]], on="ncbi_assembly")
      .merge(pd.read_csv(D / "compartment_charge.csv")[["ncbi_assembly", "secreted_DE", "cytoplasmic_DE"]], on="ncbi_assembly")
      .merge(pan, on="ncbi_assembly"))
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
B = df[(df.domain == "Bacteria") & df.ncbi_assembly.isin(tips) & df.ph_min.notna()].reset_index(drop=True)
SRC = ("thermobase", "acidb", "bacdive", "literature_auto")
for o in SRC:
    B[f"src_{o}"] = (B.origin == o).astype(float)
B["gw"] = B.acidic_excess * 100
SRCC = [f"src_{o}" for o in SRC if B[f"src_{o}"].nunique() > 1]


def pooled(cls, keys=("D", "E", "K", "R")):
    """scripts/82's pooling: sum counts and positions across membrane references, then divide."""
    tot = {k: np.zeros(len(B)) for k in ("n",) + tuple(keys)}
    used = 0
    for ref in MEM:
        if f"{ref}|{cls}|n" not in B:
            continue
        tot["n"] += B[f"{ref}|{cls}|n"].fillna(0).to_numpy(float)
        for k in keys:
            tot[k] += B[f"{ref}|{cls}|{k}"].fillna(0).to_numpy(float)
        used += 1
    return tot, used


tp, up = pooled("peri_exp")
tc, uc = pooled("cyto_exp")
ok = (tp["n"] >= MINPOS) & (tc["n"] >= MINPOS)
B = B[ok].reset_index(drop=True)
tp = {k: v[ok] for k, v in tp.items()}
tc = {k: v[ok] for k, v in tc.items()}
print(f"   the sample: {len(B)} genomes with >= {MINPOS} pooled positions on both sides ({up} periplasmic references, {uc} cytoplasmic)")

# THE PROJECT'S MEASURE, and the two GenomeSPOT-style measures. The GenomeSPOT feature is PROTEOME-WIDE and per
# predicted localization; this project's is PANEL-BASED and per exposed structural position. They are built from
# different things, which is the whole question: if they nonetheless carry the same signal, a predicted ph_min
# built on one of them cannot test the other.
X = {}
X["this project's measure, uncontrolled: panel periplasmic exposed D+E"] = 100.0 * (tp["D"] + tp["E"]) / tp["n"]
X["this project's contrast, uncontrolled: periplasmic minus cytoplasmic D+E"] = X["this project's measure, uncontrolled: panel periplasmic exposed D+E"] - 100.0 * (tc["D"] + tc["E"]) / tc["n"]
X["GenomeSPOT feature: proteome secreted minus cytoplasmic D+E"] = 100.0 * (B.secreted_DE - B.cytoplasmic_DE).to_numpy(float)
X["GenomeSPOT input: proteome secreted D+E"] = 100.0 * B.secreted_DE.to_numpy(float)
KEYS = list(X)
GS = "GenomeSPOT feature: proteome secreted minus cytoplasmic D+E"

print("\n   2a. correlation between this project's measure and GenomeSPOT's stated pH feature, across genomes")
for k in KEYS:
    if k == GS:
        continue
    r = float(np.corrcoef(X[k], X[GS])[0, 1])
    print(f"      r({k[:62]:62s}, GenomeSPOT feature) = {r:+.3f}")
r_head = float(np.corrcoef(X["this project's measure, uncontrolled: panel periplasmic exposed D+E"], X[GS])[0, 1])
r_contr = float(np.corrcoef(X["this project's contrast, uncontrolled: periplasmic minus cytoplasmic D+E"], X[GS])[0, 1])

print("\n   2b. does GenomeSPOT's feature carry the ph_min signal on its own? If it does, its predicted ph_min is not independent evidence.")
# NO K AND R CONTROL HERE, DELIBERATELY, and the figures below are therefore NOT the headline. The headline is
# periplasmic exposed D+E with lysine and arginine held fixed, +0.10540/pt (scripts/84, asserted in scripts/86);
# entering a residue share without its correlate is correction 14's error wherever the coefficient is read as that
# residue's effect. It is the right specification for THIS question only, because the quantity being interrogated is
# GenomeSPOT's feature as GenomeSPOT actually computes it - an uncontrolled compositional contrast - and a controlled
# comparison would be measuring something its predictor does not contain. Read these four numbers only against each
# other; none of them is quotable as an effect of acidic residues.
print("      (no K/R control in these four fits, deliberately - see the comment above; none of these is the headline,")
print("       which is +0.10540/pt with lysine and arginine held fixed)")
Cm = ph.vcv_fast(tree, list(B.ncbi_assembly))
y = B.ph_min.to_numpy(float)
# X0 is intercept + growth-optimum proxy + source dummies. Lambda is fitted ONCE on this covariate-only model and
# reused, which is scripts/85's correction: fitting it on a model that already contains the tested predictor makes
# the comparison between predictors depend on which predictor it was fitted with.
X0 = np.column_stack([np.ones(len(B)), B.gw.to_numpy(float)] + [B[c].to_numpy(float) for c in SRCC])
lam = ph.pgls(y, X0, Cm)["lambda"]
print(f"      Pagel's lambda {lam:.3f}, fitted on the covariate-only model and held fixed across the four fits")
fits = {}
for k in KEYS:
    f = ph.pgls(y, np.column_stack([X0, X[k]]), Cm, lam=lam)
    b, se, p = float(f["beta"][-1]), float(f["se"][-1]), float(f["p"][-1])
    fits[k] = (b, se, p)
    sd = float(np.std(X[k], ddof=1))
    print(f"      {k:54s} {b:+.5f}/pt (se {se:.5f}, p {p:.2e})  per SD {b * sd:+.4f}")
    rows.append({"part": "2b", "measure": k, "n": len(B), "beta_per_pt": b, "se": se,
                 "p": p, "sd": sd, "beta_per_sd": b * sd, "r_with_genomespot_feature": float(np.corrcoef(X[k], X[GS])[0, 1])})

# THE VERDICT, GATED. Circularity is established by BOTH halves: the measures must be correlated AND GenomeSPOT's
# feature must independently predict ph_min. A high correlation with a null feature would mean the shared part is
# noise; a predictive feature uncorrelated with the headline would be independent evidence and worth having.
gs_b, gs_se, gs_p = fits[GS]
shared = max(abs(r_head), abs(r_contr))
if gs_p < 0.05 and shared > 0.5:
    verdict = (f"CIRCULAR. GenomeSPOT's stated pH feature predicts ph_min here on its own ({gs_b:+.5f}/pt, p {gs_p:.2e}) and correlates with this "
               f"project's measure at r {shared:+.3f}. A predicted ph_min built from it cannot test the measure it is built from.")
elif gs_p >= 0.05:
    verdict = (f"NOT CIRCULAR BY THIS TEST: GenomeSPOT's feature does not predict ph_min in this sample (p {gs_p:.2e}), so its predictions "
               "cannot be inheriting the headline's signal. The usability question would then turn on label quality instead.")
else:
    verdict = (f"PARTIAL: GenomeSPOT's feature predicts ph_min ({gs_p:.2e}) but shares only r {shared:+.3f} with this project's measure, so its "
               "predictions are closer to independent evidence than to a restatement. Treat as a weak comparator, not an outcome.")
print(f"\n   verdict: {verdict}")
rows.append({"part": "2-verdict", "measure": "circularity", "n": len(B), "r_headline": r_head,
             "r_contrast": r_contr, "genomespot_p": gs_p, "verdict": verdict})

print("\n3. what the search turned up instead: MEASURED ph_min from two sources docs/DATABASES.md lists as not checked")
mt = pd.read_csv(MT, sep="\t")
mt = mt[mt.trait_name == "pH minimum"].copy()
mt["val"] = pd.to_numeric(mt.minimum, errors="coerce")
mt["obs"] = pd.to_numeric(mt.total_observations, errors="coerce")
mt["db"] = mt.databases.str.replace("=2025-07-29", "", regex=False)
print(f"   {len(mt)} species-level rows carry a measured pH minimum; {int(mt.val.notna().sum())} have a numeric value")
print("   by contributing database:")
for k, v in mt.db.value_counts().items():
    print(f"      {k:46s} {v:5d}")
no_bd = mt[~mt.db.str.contains("BacDive")]
print(f"   {len(no_bd)} rows come from JGI/GOLD or JGI/IMG with NO BacDive record, i.e. they are not reachable through the source this project already uses")


def norm(s):
    s = re.sub(r"\s*\(ex .*?\)", "", str(s)).strip()
    s = re.sub(r"^Candidatus ", "", s)
    return " ".join(s.split()[:2]).lower()


pheno = pd.read_csv(ROOT / "data/curated/phenotypes.csv")
have = {norm(x) for c in ("species", "strain_id") if c in pheno for x in pheno[c].dropna()}
mt["key"] = mt.taxon_name.map(norm)
mt["is_new"] = ~mt.key.isin(have)
print(f"\n   against this project's {len(pheno)} curated phenotype rows ({len(have)} distinct species keys):")
print(f"      {int((~mt.is_new).sum())} already present, {int(mt.is_new.sum())} NEW")
print("   new species by acidity, which is the band that matters:")
for thr in (3.0, 3.5, 4.0, 4.5, 5.0):
    n = mt[mt.is_new & (mt.val <= thr)]
    nj = n[~n.db.str.contains("BacDive")]
    print(f"      ph_min <= {thr}: {len(n):4d} new species ({len(nj):3d} of them absent from BacDive as well)")
    rows.append({"part": "3", "measure": f"new species with measured ph_min <= {thr}", "n": len(n), "n_non_bacdive": len(nj)})

print("\n4. the quality objections, measured. Nothing in part 3 is usable until these are answered per row.")
nm = mt[mt.is_new & mt.val.notna()]
one = int((nm.obs <= 1).sum())
print(f"   4a. {one} of {len(nm)} new rows ({100.0 * one / len(nm):.0f}%) rest on a SINGLE observation, so there is no within-species agreement to check")
print(f"       and no database records bound type, so every one of them enters at confidence `low` under this project's existing rule")
# 4b. THE VALUES THAT WOULD HAVE TO BE WRONG. A growth minimum at or below pH 2 is a strong claim: it is near the
# limit of the known bacterial range, and in this table it sits beside genera that are not acidophiles at all. This
# is the "lowest value tested" confusion docs/DATABASES.md warns about, surfacing as data.
# THE FLAG BELOW IS A NAME HINT AND IT UNDER-DETECTS, which is said here because the printed list would otherwise
# read as nine doubtful rows and it is not. At least two of the organisms it fails to flag are real acidophiles -
# Metallibacterium scheffleri came from acid mine drainage and Candidatus Aramenus sulfurataquae from an acidic hot
# spring - so for them a pH 2.0 minimum is plausible and the flag is simply wrong. A lineage name is not a phenotype,
# and the only thing this part establishes is WHICH rows need a primary source, never which ones are false.
ACIDOPHILE_HINT = ("acid", "sulfolobus", "picrophilus", "thermoplasma", "ferroplasma", "leptospirillum", "metallosphaera",
                   "saccharolobus", "acidilobus", "sulfurisphaera", "alicyclobacillus", "granulicella", "edaphobacter")
ext = nm[nm.val <= 2.0].copy()
ext["looks_acidophilic"] = ext.key.str.contains("|".join(ACIDOPHILE_HINT))
print(f"   4b. {len(ext)} new rows claim a growth minimum at or below pH 2.0; {int(ext.looks_acidophilic.sum())} sit in a lineage whose name suggests acidophily and "
      f"{int((~ext.looks_acidophilic).sum())} do not:")
for _, r in ext[~ext.looks_acidophilic].sort_values("val").iterrows():
    print(f"       {r.val:5.2f}  {r.taxon_name:44s} obs={int(r.obs):2d}  {r.db}")
print("       those are the rows a primary source has to confirm or kill, in BOTH directions: a pH 2.0 minimum for a marine")
print("       organism is far more likely to be the lowest pH anyone tested than a demonstrated limit, while at least two")
print("       of the nine are genuine acid-mine or hot-spring isolates that the name hint simply failed to recognise.")
rows.append({"part": "4", "measure": "new rows resting on a single observation", "n": one, "n_total": len(nm)})
rows.append({"part": "4", "measure": "new rows claiming ph_min <= 2.0 outside an acidophile-suggesting lineage",
             "n": int((~ext.looks_acidophilic).sum()), "n_total": len(ext)})

print("\n5. what this does and does not change")
print("   it CLOSES the 'possible comparator' line in docs/DATABASES.md: GenomeSPOT and everything downstream of it is not a")
print("   comparator for this project at all, and recording it as one invites exactly the circular test part 2 rules out.")
print("   it OPENS a curation route that did not exist in that table: two culture-based sources, reachable in one bulk file,")
print(f"   carrying {int(mt.is_new.sum())} species this project has never seen and {len(mt[mt.is_new & (mt.val <= 4.0)])} of them at ph_min <= 4.0.")
print("   it does NOT deliver a single usable row. The export is species-level and this project's unit is a strain with an")
print("   assembly, so each row needs a strain and a primary source before it is anything. The honest summary is that the")
print("   measured sample's bottleneck was never discovery - it was bound type, and no database has ever recorded it.")
print("   AND THERE IS ONE GENUINE GAIN THAT IS NOT DATA. GenomeSPOT is an independent group, on a different label set")
print("   (BacDive, 1,020 genomes with pH), arriving from the opposite direction - prediction rather than inference - at the")
print("   same quantity: an extracellular-minus-intracellular acidic composition contrast, their strongest residue-level pH")
print("   correlate at rho 0.56. That is external corroboration of the mechanism even though it is useless as a label, and it")
print("   is the first candidate from outside this dataset that docs/SYNTHESIS.md section 8 asks for.")
print("   THE CAVEAT ON THAT CORROBORATION, stated because it is easy to overclaim: their rho 0.56 is against pH OPTIMUM, not")
print("   ph_min, and pH is their weakest prediction (R2 0.48, RMSE 1.1 pH) and is 'systematically overpredicted below pH 5' -")
print("   which is the whole range this project is about. It supports the mechanism; it does not support a number.")

pd.DataFrame(rows).to_csv(D / "predicted_ph_sources.csv", index=False)
print(f"\nwrote {D / 'predicted_ph_sources.csv'}")
