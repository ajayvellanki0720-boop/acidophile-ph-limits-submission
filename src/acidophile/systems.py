"""Stage 3: catalog of candidate acid-limiting systems and annotation-based copy counting.

Detection here is annotation-keyword matching and is PROVISIONAL. As with the thermophile trait calls,
no system may be used in a claim until HMM / reference-sequence validation (with a coverage check) is run.
"""
from __future__ import annotations

import re

import pandas as pd

# category -> what failing here would mean for the lower pH limit
CATEGORIES = {
    "donnan_k": "reversed (inside-positive) membrane potential via K+ uptake limits proton influx",
    "proton_handling": "proton extrusion / respiratory and ATP-synthase proton handling",
    "acid_consumption": "decarboxylase / urease / deiminase systems that consume protons (acid-tolerant neutralophile style)",
    "membrane_barrier": "low proton permeability membranes (tetraether, hopanoid, cyclopropane)",
    "envelope": "S-layer, exopolysaccharide, cell wall",
    "proteostasis_repair": "chaperones, proteases, DNA repair",
    "oxidative_metal": "oxidative stress and metal resistance (low pH mobilises metals)",
    "acid_energy": "iron/sulfur oxidation energetics available only in acid",
    "organic_acid": "weak-acid (uncoupler) catabolism",
}

# name: (category, gene regexes, product regexes, product exclusion regex or None)
SYSTEMS = {
    "kdp": ("donnan_k", [r"^kdp[ABCDF]$"], [r"potassium-transporting ATPase"], None),
    "trk_ktr": ("donnan_k", [r"^(trk|ktr)[ABCDEH]$"], [r"potassium uptake protein|trk system potassium"], None),
    "kef": ("donnan_k", [r"^kef[ABCFG]$"], [r"glutathione-regulated potassium-efflux"], None),
    "f_atp_synthase": ("proton_handling", [r"^atp[ABCDEFGHI]$"], [r"F0F1|F1F0|ATP synthase (subunit|F1|F0|epsilon|gamma|alpha|beta)"], r"V-type|A-type|archaeal"),
    "v_a_atp_synthase": ("proton_handling", [r"^(ntp|atp)[ABCDEFIKH]$"], [r"V-type ATP synthase|A-type ATP synthase|archaeal ATP synthase"], None),
    "nuo_complex_i": ("proton_handling", [r"^nuo[A-N]$"], [r"NADH-quinone oxidoreductase|NADH dehydrogenase \(quinone\)"], None),
    "cytochrome_oxidase": ("proton_handling", [r"^(cox|cyo|cyd|cta|qox)[A-F]$"], [r"cytochrome (c|bd|bo3|aa3|ba3)? ?oxidase|quinol oxidase"], None),
    "p_type_atpase": ("proton_handling", [r"^(cop[AB]|zntA|cadA|ctpA|pma)$"], [r"P-type (ATPase|Cu|Zn|Cd|HMA)|heavy metal translocating P-type"], None),
    "gad": ("acid_consumption", [r"^gad[ABC]$"], [r"glutamate decarboxylase|glutamate:gamma-aminobutyrate antiporter"], None),
    "adi": ("acid_consumption", [r"^adi[AC]$"], [r"arginine decarboxylase|arginine:agmatine antiporter"], None),
    "cad": ("acid_consumption", [r"^cad[AB]$"], [r"lysine decarboxylase|lysine:cadaverine antiporter"], r"cadmium"),
    "urease": ("acid_consumption", [r"^ure[ABCDEFG]$"], [r"urease"], None),
    "adi_arc": ("acid_consumption", [r"^arcA$"], [r"arginine deiminase"], None),
    "carbonic_anhydrase": ("acid_consumption", [r"^(can|cynT|cah)$"], [r"carbonic anhydrase"], None),
    "gdgt_ring_synthase": ("membrane_barrier", [r"^(grsAB|tes|grs[AB])$"], [r"GDGT|tetraether|glycerol dibiphytanyl glycerol tetraether"], None),
    "hopanoid_shc": ("membrane_barrier", [r"^shc$"], [r"squalene-hopene cyclase|squalene--hopene cyclase"], None),
    "cfa": ("membrane_barrier", [r"^cfa$"], [r"cyclopropane[- ]fatty[- ]acyl[- ]phospholipid synthase"], None),
    "s_layer": ("envelope", [r"^slp[A-Z]?$"], [r"S-layer"], None),
    "eps_biosynth": ("envelope", [r"^(eps|wz[abcxy])[A-Z]?$"], [r"exopolysaccharide|polysaccharide biosynthesis|capsular"], None),
    "chaperones": ("proteostasis_repair", [r"^(groE[LS]|dnaK|dnaJ|grpE|clpB|htpG|hsp20|ibpA|thsA)$"], [r"chaperon|heat shock protein"], None),
    "proteases": ("proteostasis_repair", [r"^(lon|clp[PXA]|ftsH|hslV|hslU)$"], [r"ATP-dependent (Clp|Lon|zinc metallo)protease|proteasome"], None),
    "dna_repair": ("proteostasis_repair", [r"^(recA|uvr[ABC]|mutS|mutL|phr|rad[A-D]|nth|xthA)$"], [r"excinuclease|DNA repair|photolyase"], None),
    "oxidative": ("oxidative_metal", [r"^(sod[ABC]|kat[AEG]|ahpC|trx[AB]|msrA|prx)$"], [r"superoxide dismutase|catalase|peroxiredoxin|thioredoxin"], None),
    "copper_resistance": ("oxidative_metal", [r"^(cop[ABCDZ]|cus[ABCF]|cut[ACE])$"], [r"copper (resistance|efflux|chaperone)"], None),
    "arsenic_resistance": ("oxidative_metal", [r"^ars[ABCDR]$"], [r"arsenical|arsenate reductase|arsenite"], None),
    "cdf_metal_efflux": ("oxidative_metal", [r"^(czc[ABCD]|cdf|zitB|fieF)$"], [r"cation diffusion facilitator|cobalt-zinc-cadmium"], None),
    "rusticyanin_iron_ox": ("acid_energy", [r"^(rus|cyc1|cyc2|cyt579)$"], [r"rusticyanin|iron oxidase|cytochrome 579"], None),
    "sulfur_oxidation": ("acid_energy", [r"^(sox[A-Z]|sqr|tetH|sor|doxD)$"], [r"sulfur oxygenase|sulfide:quinone|tetrathionate hydrolase|thiosulfate oxidation"], None),
    "organic_acid_catabolism": ("organic_acid", [r"^(acs|fdh[A-Z]?|fdo[A-Z]|pta|ackA)$"], [r"acetyl-CoA synthetase|formate dehydrogenase|phosphate acetyltransferase|acetate kinase"], None),
}

_COMPILED = {
    name: (cat, [re.compile(g) for g in genes], [re.compile(p, re.I) for p in prods],
           re.compile(ex, re.I) if ex else None)
    for name, (cat, genes, prods, ex) in SYSTEMS.items()
}


def match_system(gene, product):
    """Return the list of system names a single annotation row matches."""
    gene = "" if pd.isna(gene) else str(gene)
    product = "" if pd.isna(product) else str(product)
    hits = []
    for name, (_, gres, pres, ex) in _COMPILED.items():
        if ex and ex.search(product):
            continue
        if any(g.match(gene) for g in gres) or any(p.search(product) for p in pres):
            hits.append(name)
    return hits


def count_systems(annotations: pd.DataFrame) -> pd.DataFrame:
    """annotations: genome_id, protein_id, gene, product. One row per genome; copy numbers + per-1000-gene rates."""
    need = {"genome_id", "protein_id", "gene", "product"}
    if not need <= set(annotations.columns):
        raise ValueError(f"annotations need columns {sorted(need)}")
    rows = []
    for gid, g in annotations.groupby("genome_id"):
        counts = dict.fromkeys(SYSTEMS, 0)
        for r in g.itertuples():
            for s in match_system(r.gene, r.product):
                counts[s] += 1
        row = {"genome_id": gid, "n_genes": g["protein_id"].nunique()}
        for s, c in counts.items():
            row[f"cn_{s}"] = c
            row[f"rate_{s}"] = 1000.0 * c / row["n_genes"]
        rows.append(row)
    return pd.DataFrame(rows)


def pathway_scores(counts: pd.DataFrame) -> pd.DataFrame:
    """Per category: fraction of member systems present (>0 copies) -> pw_<category>."""
    out = counts[["genome_id"]].copy()
    for cat in CATEGORIES:
        members = [s for s, v in SYSTEMS.items() if v[0] == cat]
        out[f"pw_{cat}"] = (counts[[f"cn_{s}" for s in members]] > 0).mean(axis=1)
    return out
