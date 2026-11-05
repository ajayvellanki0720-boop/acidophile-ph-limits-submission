"""Hand-curated pH values, read from the Europe PMC abstracts in data/curated/epmc_ph_review.csv.

Every number below was read from the quoted abstract statement, not taken from the regex. Convention:
  - A species-description sentence of the form "grew at pH a-b (optimum c)" is recorded as bound type
    `measured`, confidence `medium`: the abstract does not state that growth failed beyond a/b, so it cannot be `high`.
  - If the wording shows tested points rather than limits ("also occurred at", "media maintained between",
    "good growth at"), the bound is `lowest_tested`/`highest_tested` and the row stays out of the primary set.
  - Optimum ranges are stored as their midpoint.
Entries: pmid -> (ph_min, ph_opt, ph_max, bound_override, note). bound_override: None | "censored".
Usage: python3 scripts/00f_curate_literature.py   (writes data/curated/literature_curated.csv)"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
C = "censored"
CURATED = {
    "31995463": (3.5, 6.0, 8.0, None, ""),
    "17066338": (1.0, 1.35, 5.0, None, ""),
    "17625168": (0.35, 1.1, 3.0, None, "optimum 0.8-1.4; range 0.35-3.0"),
    "27902281": (4.0, 5.0, 5.5, None, ""),
    "25273513": (2.0, 3.5, 5.5, None, "multi-strain group; type strain MS8T"),
    "26475169": (3.5, 4.25, 8.0, None, ""),
    "19628615": (0.0, 1.5, 4.0, None, ""),
    "38805027": (2.0, 2.5, 4.5, None, ""),
    "42061222": (0.4, 2.55, 6.5, None, ""),
    "27692038": (4.5, 6.5, 7.0, None, ""),
    "19622665": (2.0, 3.5, 6.0, None, ""),
    "24711593": (3.5, 4.8, 6.5, None, ""),
    "35312474": (2.5, 3.5, 4.5, None, "A. curvatus protologue, type strain ALEF1T (full text PMC9558579); pH tested 1.0-5.0 in 0.5 steps"),
    "25505343": (3.0, 4.0, 6.0, None, ""),
    "26476812": (2.5, 3.5, 5.5, None, ""),
    "29557767": (1.5, 3.0, 4.5, None, ""),
    "18175712": (4.5, 5.5, 7.5, None, ""),
    "25037879": (2.0, 3.2, 6.0, None, ""),
    "15388732": (2.0, 4.0, 6.0, None, ""),
    "28857027": (4.5, 5.75, 7.5, C, "growth 'also occurred at pH 4.5 and 7.5': tested points"),
    "16957103": (4.5, 5.8, 8.0, None, ""),
    "19819996": (3.0, None, 8.0, C, "'good growth at pH 3.0-8.0'"),
    "20709915": (4.0, None, 8.0, C, "'good growth at pH 4.0-8.0'"),
    "19666801": (3.0, None, 7.0, C, "'good growth between pH 3.0 and 7.0'"),
    "25323596": (3.5, 5.5, 8.3, None, ""),
    "20358236": (3.6, 5.2, 5.5, None, ""),
    "34255623": (4.0, 5.5, 7.0, None, ""),
    "39869511": (4.0, 6.5, 7.5, None, ""),
    "40875059": (3.5, 5.25, 7.5, None, ""),
    "15388730": (5.2, 6.0, None, None, "'minimal pH for growth of 5.2'"),
    "28945535": (3.0, 4.75, 7.0, None, "optimum 4.0-5.5 for type strain 4G-K17T"),
    "26503862": (3.5, 4.0, 5.5, None, ""),
    "39196616": (3.6, 5.05, 7.3, None, ""),
    "34235583": (4.0, 5.0, 6.0, None, ""),
    "30775959": (3.5, 4.75, 7.0, None, ""),
    "28875923": (4.5, 5.5, 6.5, None, ""),
    "23771625": (4.5, 5.35, 6.2, None, ""),
    "23435243": (4.5, None, 6.8, None, ""),
    "22753525": (4.5, 5.35, 7.2, None, ""),
    "11594618": (5.0, None, 7.0, None, ""),
    "25404480": (1.5, None, 6.5, None, "species description covers several strains"),
    "37755432": (1.0, 2.5, 6.0, None, ""),
    "13130000": (4.5, 5.5, 7.0, None, ""),
    "14742473": (4.2, 5.75, 7.5, None, ""),
    "10843033": (4.5, 5.25, 7.0, None, "range from text 'growth between pH 4.5 and pH 7.0'"),
    "29858827": (4.0, 6.5, 9.0, None, ""),
    "41989600": (3.9, 3.9, 8.3, None, "abstract states optimum pH 3.9"),
    "16738120": (4.8, 5.35, 7.0, None, ""),
    "23125319": (3.8, 5.25, 9.5, None, ""),
    "21335501": (3.2, 4.9, 7.1, None, ""),
    "41525133": (5.0, None, 7.0, None, ""),
    "30896386": (4.0, 6.5, 9.0, None, ""),
    "26297487": (4.0, 5.0, 9.0, None, ""),
    "35951006": (1.0, 3.0, 5.5, None, ""),
    "34236956": (1.5, 3.0, 5.0, None, ""),
    "31958046": (0.4, 3.5, 5.5, None, ""),
    "32375939": (4.3, 6.6, 7.4, None, ""),
    "26908287": (4.6, 6.05, 8.1, None, ""),
    "29671720": (2.5, 3.75, 6.0, None, ""),
    "21186292": (4.5, None, 7.5, None, ""),
    "24096353": (3.5, 5.0, 7.0, None, ""),
    "23990651": (6.0, 6.0, 7.0, None, ""),
    "16341463": (4.0, None, 7.5, None, ""),
    "25323594": (4.0, 5.0, 6.0, None, ""),
    "16432746": (1.75, 2.75, 3.0, C, "'media maintained at pH between 1.75 and 3.0': tested interval"),
    "12892143": (2.3, 3.75, 5.4, None, ""),
    "20639230": (2.0, 3.0, 4.5, None, ""),
    "19965997": (4.0, 4.45, 6.0, None, ""),
    "10843038": (1.3, 1.7, 2.2, None, ""),
    "31355740": (2.5, 4.0, 5.0, None, ""),
    "21551329": (3.2, 4.95, 6.6, None, ""),
    "21460138": (3.0, 4.75, 7.5, None, ""),
    "17978217": (1.8, 3.0, 4.0, None, ""),
    "39190454": (0.5, 1.0, 4.0, None, ""),
    "38792683": (3.5, 4.0, 6.0, None, ""),
    "10425774": (2.3, 3.95, 6.4, None, ""),
    "9734041": (2.6, 4.0, 5.9, None, ""),
    "37342988": (3.5, 4.75, 5.5, None, ""),
    "37166571": (4.0, 5.0, 7.5, None, ""),
    "32658637": (3.5, 5.5, 7.5, None, ""),
    "27902294": (3.5, 5.5, 6.5, None, ""),
    "26651491": (4.5, 5.75, 7.5, None, ""),
    "25724745": (4.2, 6.25, 7.1, None, ""),
    "24867171": (4.0, 5.75, 7.2, None, ""),
    "24048862": (4.1, 6.5, 7.8, None, ""),
    "19622650": (3.1, 4.75, 6.5, None, ""),
    "11837310": (4.2, None, 7.2, None, ""),
    "19196778": (3.8, 5.75, 7.2, None, ""),
    "17978240": (4.2, 5.6, 7.5, None, ""),
    "19837732": (4.5, 5.75, 7.5, None, ""),
    "18450711": (4.2, 5.6, 7.5, None, ""),
    "29627890": (4.5, None, 10.5, None, ""),
    "22228140": (5.0, 7.0, 10.0, None, "marine neutrophile; included because ph_min <= 5"),
    "10555364": (2.0, 4.0, 6.0, None, ""),
    "7801633": (1.0, 2.0, 5.0, None, ""),
    "15653874": (4.2, 5.2, 7.0, C, "'(initial) pH'"),
}
# Second batch: values read from full texts / targeted searches. (taxon, pmid, min, opt, max, bound_override, note, temp_opt)
EXTRA = [
    ("Metallosphaera cuprina", "21057050", 2.5, 3.5, 5.5, None, "type strain Ar-4T", 65.0),
    ("Acididesulfobacillus acetoxydans", "35391737", 3.8, 5.0, 6.5, None, "unable to grow at pH 7.0 (full text)", None),
    ("Methylacidimicrobium thermophilum", "33643272", 1.5, 4.0, 5.5, None, "full text: 'not able to grow at pH below 1.5'; growth still seen at 5.5", 50.0),
    ("Granulicella tundricola", "25197431", 3.5, 5.0, 6.5, None, "genome-report statement for type strain MP5ACTX9T", None),
    ("Kyrpidia tusciae", "22180816", 4.2, None, 7.5, None, "genome-report statement for T2T (as Bacillus tusciae)", 55.0),
    ("Streptacidiphilus hamsterleyensis", "23989983", 4.5, 5.5, 6.0, None, "", None),
    ("Actinospica durhamensis", "26026797", 4.0, 5.5, 6.0, None, "", None),
    ("Herminiimonas arsenitoxidans", "29101383", 3.0, 5.5, 9.0, None, "", None),
    ("Frigoriglobus tundricola", "32847778", 4.2, 5.25, 6.8, None, "", None),
    ("Limnoglobus roseus", "31800383", 5.0, 6.5, 7.5, None, "", None),
    ("Marinobacter subterrani", "26236300", 5.0, 6.0, 9.0, None, "", None),
    ("Methanosuratincola petrocarbonis", "40637718", 5.0, 6.25, 8.0, None, "", 55.0),
    ("Endosaccharibacter trunci", "39005926", 3.0, 4.75, 8.0, None, "two species in one paper, assigned by title order (ASSUMPTION)", None),
    ("Rhizosaccharibacter radicis", "39005926", 3.5, 4.25, 8.5, None, "two species in one paper, assigned by title order (ASSUMPTION)", None),
    ("Alicyclobacillus sendaiensis", "12892130", 2.5, 5.5, 6.5, None, "original-description abstract: pH range for growth 2.5-6.5, optimum 5.5 (matches ThermoBase)", 55.0),
    ("Sulfodiicoccus acidiphilus", "28629504", 1.4, 3.25, 5.5, None, "", None),
    ("Bryobacter aggregatus", "19651730", 4.5, 6.0, 7.2, None, "", None),
    ("Acidicapsa dinghuensis", "29873630", 4.0, 4.75, 6.5, None, "", None),
    ("Thiomonas islandica", "20173005", 4.0, 4.0, 7.0, None, "", 45.0),
    ("Methylocystis silviterrae", "34913862", 4.5, 6.25, 7.5, None, "", None),
    ("Acidibacter ferrireducens", "25116055", 2.5, 3.75, 4.5, None, "", None),
    ("Methylomonas paludis", "23159751", 3.8, 6.1, 7.3, None, "", None),
    ("Methylocystis bryophila", "22707532", 4.2, 6.25, 7.6, None, "", None),
    ("Acidisoma cladoniae", "39320523", 4.0, 5.25, 7.5, None, "", None),
    ("Mesoaciditoga lauensis", "23959829", 4.1, 5.6, 6.0, None, "", 58.0),
    ("Methylorosula polaris", "21551336", 4.0, 5.75, 7.8, None, "", None),
    ("Thiobacter aerophilum", "39597641", 4.8, 5.35, 7.0, None, "", None),
    ("Alicyclobacillus acidiphilus", "12361274", 2.5, 3.0, 5.5, None, "", None),
    ("Alicyclobacillus pomorum", "13130044", 3.0, 4.25, 6.0, None, "", None),
    ("Clostridium uliginosum", "11411680", 4.0, 6.5, 9.0, None, "", None),
    ("Alicyclobacillus herbarius", "11837292", 3.5, 4.75, 6.0, None, "", None),
    ("Ferroplasma cupricumulans", "16721487", 0.4, 1.1, None, "censored", "'lowest pH at which growth was observed was 0.4'", None),
    ("Alicyclobacillus ferrooxydans", "19060079", 2.0, 3.0, 6.0, None, "", 28.0),
    ("Desulfurella amilsii", "26704766", 3.0, 6.25, 7.0, None, "", None),
    ("Burkholderia dabaoshanensis", "23226514", 4.0, None, 9.0, None, "", None),
    ("Methanobacterium paludis", "24449792", 4.7, None, 7.1, None, "", None),
    ("Caenibacillus caldisaponilyticus", "27117268", 4.8, 6.2, 8.1, None, "", None),
]
# Lowest pH actually tested (from full-text methods). If tested_lo < ph_min the minimum is a verified limit, not a tested edge.
TESTED_LO = {"38805027": 0.5, "39196616": 2.5, "39869511": 3.0, "35312474": 1.0, "39190454": 0.0}
EXPLICIT_MIN = {"33643272": "full text: not able to grow at pH below 1.5"}

# Reviewed and deliberately left out (reason): optimum-only statements, isolation-site pH mistaken for a range,
# co-cultures, non-acidophiles with neutral optimum and min >5, multi-strain ambiguity, eukaryote, virus.
EXCLUDED = {
    "10758893": "no species-level range (single growth point)",
    "19643887": "optimum only",
    "17720836": "optimum only",
    "19196772": "optimum only, hedged",
    "42377362": "Candidatus, no taxon match",
    "39005906": "min >5, not acidophilic",
    "28068215": "min >5, neutral optimum",
    "28553696": "min >5, neutral optimum",
    "24591424": "alkaline",
    "29939126": "min >5, neutral optimum",
    "23064351": "min >5, neutral optimum",
    "41506132": "multi-strain abstract with different optima",
    "28646634": "min >5, neutral optimum",
    "22544786": "min >5, neutral optimum",
    "29465339": "optimum only, hedged",
    "35993221": "obligate ectosymbiont, co-culture",
    "37791995": "optimum only",
    "32926197": "min >5",
    "25701846": "min >5, neutral optimum",
    "22984136": "min >5, neutral optimum",
    "19605731": "min >5, neutral optimum",
    "12148613": "optimum only",
    "39805719": "min >5",
    "30631314": "co-culture",
    "16957104": "isolation-site pH, not growth range",
    "17625173": "isolation-site pH, not growth range",
    "21462554": "eukaryote",
    "16357204": "virus",
    "16348110": "optimum only",
    "39348174": "range not stated unambiguously in abstract",
    "PPR285527": "preprint duplicate of 34235583",
    "24860109": "optimum only",
    "28829024": "optimum only",
    "24425740": "min >5",
    "22529844": "isolation-site pH",
    "21097638": "isolation-site pH",
    "16957124": "optimum only",
    "17329771": "isolation-site pH",
    "19620354": "two species share one pH statement",
    "34683374": "two species share one pH statement",
    "22058325": "four species share one pH statement",
    "20118293": "several species/strains share one pH statement",
}

rev = pd.read_csv(ROOT / "data/curated/epmc_ph_review.csv", dtype={"pmid": str}).drop_duplicates("pmid")
rows = []
for pmid, (lo, opt, hi, ov, note) in CURATED.items():
    r = rev[rev.pmid == pmid]
    if r.empty:
        print("MISSING pmid", pmid); continue
    r = r.iloc[0]
    taxon = r.taxon
    if pmid == "35312474":
        rows.append({"taxon": "Alicyclobacillus mengziensis", "pmid": pmid, "ph_min": 2.0, "ph_opt": 4.5, "ph_max": 5.5,
                     "bound": "measured", "assumption": False, "note": "A. mengziensis protologue, type strain S30H14T (full text PMC9558579); pH tested 1.0-5.0 in 0.5 steps",
                     "statement": r.statement, "temp_opt_c": 30.0})
    rows.append({"taxon": taxon, "pmid": pmid, "ph_min": lo, "ph_opt": opt, "ph_max": hi,
                 "bound": "censored" if ov else "measured", "assumption": "ASSUMPTION" in note, "note": note, "statement": r.statement, "temp_opt_c": r.temp_opt})
for tx, pmid, lo, opt, hi, ov, note, topt in EXTRA:
    rows.append({"taxon": tx, "pmid": pmid, "ph_min": lo, "ph_opt": opt, "ph_max": hi, "bound": "censored" if ov else "measured",
                 "assumption": "ASSUMPTION" in note, "note": note, "statement": "", "temp_opt_c": topt})
TAXON_FIX = {"39190454": "Oxyplasma meridianum", "37166571": "Lichenifustis flavocetrariae"}   # names read from the paper titles
out = pd.DataFrame(rows)
out["taxon"] = [TAXON_FIX.get(p, tx) if pd.isna(tx) else tx for p, tx in zip(out.pmid, out.taxon)]
assert out.taxon.notna().all(), out[out.taxon.isna()]
out["tested_lo"] = out.pmid.map(TESTED_LO)
out["min_verified"] = ((out.tested_lo < out.ph_min) | out.pmid.isin(EXPLICIT_MIN)) & (out.bound == "measured")
out.to_csv(ROOT / "data/curated/literature_curated.csv", index=False)
print(len(out), "curated rows;", (out.bound == "measured").sum(), "measured bounds;", out.ph_opt.notna().sum(), "with optimum")
unreviewed = set(rev.pmid) - set(CURATED) - set(EXCLUDED)
print("review rows neither curated nor excluded:", len(unreviewed))
