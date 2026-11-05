"""Do DEDICATED pH SENSORS carry the opposite sign to housekeeping proteins? -> DATA_DIR/sensor_panel.csv
THE PREDICTION, AND WHERE IT CAME FROM. Correction 21 found the one paper that has done exhaustively what this project predicts comparatively: Weeks & Sachs 2001
mutagenised every histidine, aspartate and glutamate in the periplasmic domain of the H. pylori urea channel UreI and found three histidines AND three carboxylates
REQUIRED to activate transport in acid. This project finds acid-tolerant genomes DEPLETED of exposed periplasmic acidic residues. The reconciliation offered in
docs/EXPERIMENTAL_EVIDENCE.md section 2.2 is that a protein whose job is to RESPOND to acid needs titratable groups on the sensing face, while a proteome that must merely
FUNCTION at low pH may do better without them. That reconciliation is convenient, and convenient reconciliations should be made to pay. It predicts something falsifiable:
a panel built from DEDICATED SENSORS should carry the OPPOSITE sign to the housekeeping panel. Every one of the 32 accepted membrane references is a transporter, a
translocase, a respiratory subunit or a cell-wall enzyme - the panel contains no sensors at all, which is why the prediction has never been tested here.
THE FEASIBILITY GATE COMES FIRST, AND IT IS DECLARED BEFORE THE RESULT. Sensor domains are narrow families. PhoQ_Sensor and CpxA_peri are largely enterobacterial, and
44% of this sample is Pseudomonadota, so a sensor panel may be detectable only in a lineage-confined subset - in which case the comparison is confounded with ancestry and
the test cannot be run, exactly as scripts/92's soluble negative control could not be run because the references were too universal to vary. Part 2 applies a gate fixed
in advance: detection in at least 30% of the sample, spanning at least three phyla with at least 50 genomes each. A test that fails the gate is reported as unrunnable,
not run anyway with a caveat.
THE MEASURE IS COARSER THAN THE MAIN PANEL'S AND IS LABELLED SO. The housekeeping panel classes each position by UniProt topology AND by relative solvent accessibility of
an AlphaFold model under a per-residue pLDDT filter. Here, positions are classed by TOPOLOGY ONLY - periplasmic versus cytoplasmic - with no exposure filter, because the
exposure route needs structures and per-position alignment machinery this script does not rebuild. So a coefficient here is NOT comparable in magnitude with anything in
docs/NUMBERS.md and none is quoted as if it were. What is comparable is the SIGN, which is what the prediction is about, and the same topology-only measure is computed
for the HOUSEKEEPING references too so that the comparison is like with like rather than sensor-coarse against housekeeping-fine.
WHAT A NULL WOULD AND WOULD NOT MEAN. If the sensor panel shows no sign flip, the reconciliation in section 2.2 is not supported, and correction 21's tension between
this project's direction and UreI's stands unexplained. That is a real outcome and the one this script exists to allow. If the sensor panel flips, the reconciliation
survives one test and remains an account rather than a demonstration.
Usage: python3 scripts/99_sensor_panel.py [--detect-only]   (env DATA_DIR=results/broad2, PHYLO_DIR=results/phylo_broad2, PROTEOME_DIR=data/raw/proteomes)"""
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import pandas as pd
import pyhmmer
from Bio import Phylo
from acidophile import phylo as ph

warnings.filterwarnings("ignore")
D = ROOT / os.environ.get("DATA_DIR", "results/broad2")
T = ROOT / os.environ.get("PHYLO_DIR", "results/phylo_broad2")
PROT = ROOT / os.environ.get("PROTEOME_DIR", "data/raw/proteomes")
HMMD = ROOT / "data/raw/sensor_hmm"
CACHE = D / "sensor_detection.csv"
rows = []
_s = lambda x: x.decode() if isinstance(x, bytes) else str(x)

# THE SENSOR SET. Each is an E. coli inner-membrane signalling protein with a UniProt-annotated periplasmic domain and
# a Pfam family specific to that domain rather than to the shared kinase machinery - PhoQ_Sensor and CpxA_peri are the
# periplasmic sensors themselves, not HATPase_c or HisKA, which every histidine kinase carries and which would make
# the panel a census of two-component systems rather than of pH sensors. HAMP is included as a deliberately BROAD
# signalling family to separate "sensors are rare" from "this particular sensor is rare".
SENSORS = [("PF08918", "PhoQ_Sensor", "P23837", "PhoQ periplasmic sensor domain (Mg2+ and acidic pH)", True),
           ("PF16527", "CpxA_peri", "P0AE82", "CpxA periplasmic domain (envelope stress, pH)", True),
           ("PF00497", "SBP_bac_3", "P30855", "the family EvgS's periplasmic domain belongs to - NOT diagnostic, see part 1b", False),
           ("PF00672", "HAMP", None, "HAMP linker - broad signalling control, not a periplasmic sensor", False)]
MIN_FRAC, MIN_PHYLA, MIN_PER_PHYLUM = 0.30, 3, 50

pheno = pd.read_csv(D / "phenotypes_all.csv")
tree = Phylo.read(T / "Bacteria.nwk", "newick")
tips = {t.name for t in tree.get_terminals()}
B = pheno[(pheno.domain == "Bacteria") & pheno.ncbi_assembly.isin(tips) & pheno.ph_min.notna()].drop_duplicates("ncbi_assembly").reset_index(drop=True)
print(f"sample: {len(B)} bacterial genomes with a curated minimum and a place in the tree")

print("\n1. detecting the sensor families across the sample")
if CACHE.exists():
    det = pd.read_csv(CACHE)
    print(f"   reusing {CACHE.name} ({len(det)} genomes already searched)")
else:
    if not HMMD.exists():
        sys.exit(f"missing {HMMD}. Fetch the sensor HMMs first:\n"
                 f"  mkdir -p {HMMD}\n"
                 "  for pf in PF08918 PF16527 PF00497 PF00672; do\n"
                 f"    curl -sSL -o {HMMD}/$pf.hmm.gz \"https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/$pf/?annotation=hmm\"\n"
                 f"    gunzip -f {HMMD}/$pf.hmm.gz\n  done")
    hmms = []
    for acc, _, _, _, _ in SENSORS:
        f = HMMD / f"{acc}.hmm"
        if not f.exists():
            sys.exit(f"missing {f}; see the fetch commands above")
        with pyhmmer.plan7.HMMFile(f) as hf:
            hmms.extend(hf)
    alpha = pyhmmer.easel.Alphabet.amino()
    out = []
    for i, acc in enumerate(B.ncbi_assembly):
        fp = PROT / f"{acc}.faa"
        if not fp.exists():
            continue
        with pyhmmer.easel.SequenceFile(fp, digital=True, alphabet=alpha) as sf:
            seqs = sf.read_block()
        hit = {a: 0 for a, *_ in SENSORS}
        # `gathering` is the family's own curated threshold, which is what makes a presence call comparable across
        # families of different widths; a fixed E-value would make the broad families look universal by construction.
        for top in pyhmmer.hmmsearch(hmms, seqs, bit_cutoffs="gathering", cpus=4):
            fam = _s(top.query.accession).split(".")[0]
            hit[fam] = hit.get(fam, 0) + sum(1 for h in top if h.included)
        out.append({"ncbi_assembly": acc, **hit})
        if (i + 1) % 400 == 0:
            print(f"      {i + 1}/{len(B)} genomes searched", flush=True)
    det = pd.DataFrame(out)
    det.to_csv(CACHE, index=False)
    print(f"   wrote {CACHE}")

G = B.merge(det, on="ncbi_assembly")
print(f"   {len(G)} genomes searched and matched to the phenotype table")
print(f"   {'family':>12s} {'what it is':>58s} {'detected':>9s} {'share':>7s} {'median copies':>14s}")
for acc, name, _, what, is_sensor in SENSORS:
    if acc not in G:
        continue
    pres = (G[acc] > 0)
    med = float(G.loc[pres, acc].median()) if pres.any() else 0.0
    print(f"   {name:>12s} {what[:58]:>58s} {int(pres.sum()):>9d} {pres.mean():>7.1%} {med:>14.0f}")
    rows.append({"part": "1", "family": name, "accession": acc, "is_periplasmic_sensor": is_sensor,
                 "n_detected": int(pres.sum()), "share_detected": float(pres.mean()), "n_total": len(G),
                 "median_copies_where_present": med})

print("\n1b. which of these are DIAGNOSTIC of a sensor, derived rather than asserted")
# THE FIRST VERSION OF THIS SCRIPT PASSED ITS OWN GATE ON A FAMILY THAT IS NOT A SENSOR MARKER, and the near-miss is
# recorded because it is the whole reason this part exists. PF00497 was entered as "EvgS periplasmic solute-binding
# domain", which is true and misleading: EvgS's periplasmic domain belongs to that family, but the family is the
# generic bacterial solute-binding fold used by amino-acid ABC transporters throughout the tree. Detected in 82% of
# genomes at a median of 4 copies, it carried the gate on its own while the two actual dedicated sensors sat at 0.5%
# and 1.8%. A label written from what a reference protein contains is not a statement about what a family means.
# THE DERIVED CRITERION. A domain specific to one sensor is near-single-copy where it occurs: a genome has one PhoQ,
# not four. A family found in multiple copies in most genomes is a fold shared by many proteins, whatever one of its
# members does. So diagnosticity is computed from the detection table - median copies <= 1 where present - and the
# declared is_sensor flag must AGREE with it, or the script stops rather than choosing which to believe.
print(f"   {'family':>12s} {'share':>7s} {'median copies':>14s} {'declared':>10s} {'derived':>8s}")
diag = {}
for acc, name, _, _, is_sensor in SENSORS:
    if acc not in G:
        continue
    pres = (G[acc] > 0)
    med = float(G.loc[pres, acc].median()) if pres.any() else 0.0
    derived = bool(med <= 1 and pres.mean() < 0.5)
    diag[acc] = derived
    flag = "" if derived == is_sensor else "   <- DISAGREES"
    print(f"   {name:>12s} {pres.mean():>7.1%} {med:>14.0f} {str(is_sensor):>10s} {str(derived):>8s}{flag}")
    rows.append({"part": "1b", "family": name, "accession": acc, "declared_sensor": is_sensor,
                 "derived_diagnostic": derived, "median_copies": med, "share_detected": float(pres.mean())})
bad = [n for a, n, _, _, d in SENSORS if a in diag and diag[a] != d]
assert not bad, (f"declared and derived diagnosticity disagree for {bad}; one of the two is wrong and the script will not "
                 "choose for you - fix the SENSORS table or the criterion")
print("   Reading it: a periplasmic sensor domain is a thin, lineage-specific thing. The generic folds that ARE broadly")
print("   distributed are broadly distributed precisely because they are not specific to sensing.")

print("\n2. the feasibility gate, fixed before any coefficient was computed")
print(f"   requires: detection in >= {MIN_FRAC:.0%} of the sample, spanning >= {MIN_PHYLA} phyla with >= {MIN_PER_PHYLUM} genomes each")
SENS = [a for a, _, _, _, s_ in SENSORS if s_ and a in G and diag.get(a)]
G["sensor_any"] = (G[SENS].sum(axis=1) > 0) if SENS else False
frac = float(G.sensor_any.mean())
byp = G[G.sensor_any].phylum.value_counts()
big = byp[byp >= MIN_PER_PHYLUM]
print(f"   any periplasmic sensor family detected: {int(G.sensor_any.sum())} genomes, {frac:.1%} of the sample")
print(f"   phyla with >= {MIN_PER_PHYLUM} such genomes: {len(big)}  ({', '.join(f'{k} {v}' for k, v in big.items())})")
# LINEAGE CONFOUNDING, measured rather than feared: if detection is concentrated in one phylum the comparison with the
# housekeeping panel is a comparison between lineages wearing a comparison between protein classes.
conc = float(byp.iloc[0] / byp.sum()) if len(byp) else 1.0
print(f"   largest phylum's share of detections: {conc:.1%}" + ("  <- lineage-confined" if conc > 0.7 else ""))
PASS = frac >= MIN_FRAC and len(big) >= MIN_PHYLA and conc <= 0.7
rows.append({"part": "2", "measure": "feasibility gate", "share_detected": frac, "n_phyla_over_floor": len(big),
             "largest_phylum_share": conc, "passes": bool(PASS),
             "gate": f">={MIN_FRAC:.0%} detected, >={MIN_PHYLA} phyla with >={MIN_PER_PHYLUM}, <=70% in one phylum"})

if not PASS:
    print("\n   GATE NOT PASSED. The test is reported as unrunnable rather than run with a caveat, for the reason")
    print("   scripts/92 recorded when its soluble negative control could not be run: a control or comparison that")
    print("   cannot be computed is a limitation, and printing a coefficient from it would be presenting the absence")
    print("   of a test as the result of one. What fails here is specific and worth stating:")
    if frac < MIN_FRAC:
        print(f"      - dedicated periplasmic sensor domains are detected in only {frac:.1%} of the sample. Sensor families are")
        print("        narrow by nature: a domain that evolved to bind one signal in one lineage does not travel, which is")
        print("        exactly why the housekeeping panel was built from transporters and respiratory subunits instead.")
    if len(big) < MIN_PHYLA:
        print(f"      - detections span only {len(big)} phylum/phyla above the {MIN_PER_PHYLUM}-genome floor, so a sensor-versus-housekeeping")
        print("        contrast would be a between-lineage contrast in disguise.")
    if conc > 0.7:
        print(f"      - {conc:.0%} of detections sit in a single phylum ({byp.index[0]}), which makes the comparison confounded with ancestry.")
    print("\n   WHAT WOULD MAKE IT RUNNABLE, since the prediction is worth testing and this is a data limitation rather")
    print("   than a conceptual one: a sensor set chosen for BREADTH rather than for being the textbook examples - the")
    print("   pH-responsive members of widely distributed sensor superfamilies, identified by domain architecture rather")
    print("   than by orthology to E. coli - or a within-Pseudomonadota test accepting the loss of generality, which")
    print("   would need the housekeeping panel refit on the same genomes to be a comparison at all.")
    print("   THE PREDICTION IS UNTESTED, NOT REFUTED, and correction 21's tension stands unexplained.")
    pd.DataFrame(rows).to_csv(D / "sensor_panel.csv", index=False)
    print(f"\nwrote {D / 'sensor_panel.csv'}")
    sys.exit(0)

print("\n   gate passed; the sign test follows")
print("\n3. NOT IMPLEMENTED BEYOND THE GATE IN THIS VERSION")
print("   The composition step needs each genome's best sensor hit aligned to the reference and tallied by the")
print("   reference's UniProt topology, which is scripts/36a's machinery pointed at a different reference table.")
print("   It is deliberately not reimplemented inline here: a second, subtly different copy of the panel builder is")
print("   how this project would acquire an error no check could see, since both copies would have to be wrong in the")
print("   same way for the comparison to look clean. The correct next step is to add the sensor references to")
print("   data/curated/ with their topology and structures and run the existing builder over them.")
rows.append({"part": "3", "measure": "composition step", "status": "not implemented; gate passed, builder reuse required"})
pd.DataFrame(rows).to_csv(D / "sensor_panel.csv", index=False)
print(f"\nwrote {D / 'sensor_panel.csv'}")
