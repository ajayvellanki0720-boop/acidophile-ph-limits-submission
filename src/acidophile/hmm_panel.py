"""HMM-based system detection (Pfam gathering cutoffs) used to validate and replace the keyword calls in systems.py.

Rules per system: `any` = protein hits any listed family; `all_protein` = one protein must hit all listed families;
`require_genome` = families that must also occur somewhere in the genome (e.g. KdpC for KdpA).
Families here are what Pfam can resolve. Gad vs Adi/Cad decarboxylases and several others are only separable at
family level, which is recorded in `CAVEATS`. Systems with no clean Pfam model are listed in `UNVALIDATABLE`."""
from __future__ import annotations

import gzip
import urllib.request
from pathlib import Path

PANEL = {
    "kdp": dict(any=["PF03814"], require_genome=["PF02669"]),
    "trk_ktr": dict(any=["PF02386"]),
    "kef": dict(all_protein=[["PF00999", "PF02254"]]),
    "f_atp_synthase": dict(any=["PF00119"]),
    "v_a_atp_synthase": dict(any=["PF01991"]),
    "nuo_complex_i": dict(any=["PF00346"]),
    "cytochrome_oxidase": dict(any=["PF00115", "PF01654"]),
    "p_type_atpase": dict(all_protein=[["PF00122", "PF00403"]]),
    "gad": dict(any=["PF00282"]),
    "adi_cad": dict(any=["PF01276"]),
    "urease": dict(any=["PF00449"]),
    "adi_arc": dict(any=["PF02274"]),
    "carbonic_anhydrase": dict(any=["PF00484", "PF00194"]),
    "hopanoid_shc": dict(any=["PF13243"]),
    "cfa": dict(any=["PF02353"]),
    "s_layer": dict(any=["PF00395"]),
    "chaperones": dict(any=["PF00118", "PF00012", "PF00011", "PF02861"]),
    "proteases": dict(any=["PF05362", "PF01434", "PF00574"]),
    "dna_repair": dict(any=["PF00154", "PF00488"]),
    "oxidative": dict(any=["PF00081", "PF00080", "PF00199", "PF00578"]),
    "cdf_metal_efflux": dict(any=["PF01545"]),
    "arsenic_resistance": dict(any=["PF03960"]),
    "rusticyanin_iron_ox": dict(any=["PF00127"]),
    "copper_resistance": dict(any=["PF00394"]),
}
UNVALIDATABLE = ["gdgt_ring_synthase", "sulfur_oxidation", "eps_biosynth", "organic_acid_catabolism", "cad", "adi"]
# keyword system -> HMM system it is compared against
KW_TO_HMM = {"kdp": "kdp", "trk_ktr": "trk_ktr", "kef": "kef", "f_atp_synthase": "f_atp_synthase",
             "v_a_atp_synthase": "v_a_atp_synthase", "nuo_complex_i": "nuo_complex_i",
             "cytochrome_oxidase": "cytochrome_oxidase", "p_type_atpase": "p_type_atpase", "gad": "gad",
             "urease": "urease", "adi_arc": "adi_arc", "carbonic_anhydrase": "carbonic_anhydrase",
             "hopanoid_shc": "hopanoid_shc", "cfa": "cfa", "s_layer": "s_layer", "chaperones": "chaperones",
             "proteases": "proteases", "dna_repair": "dna_repair", "oxidative": "oxidative",
             "cdf_metal_efflux": "cdf_metal_efflux", "arsenic_resistance": "arsenic_resistance",
             "rusticyanin_iron_ox": "rusticyanin_iron_ox", "copper_resistance": "copper_resistance"}
CAVEATS = {
    "gad": "PF00282 is the whole group-II PLP decarboxylase family (includes histidine/DOPA decarboxylases): family-level only",
    "adi_cad": "PF01276 covers Orn/Lys/Arg decarboxylases: Adi and Cad are not separable",
    "nuo_complex_i": "PF00346 also occurs in some NiFe hydrogenases",
    "rusticyanin_iron_ox": "PF00127 is the whole cupredoxin family, broader than rusticyanin",
    "copper_resistance": "PF00394 is multicopper oxidase generally, not only CopA",
    "chaperones": "Cpn60/HSP70/HSP20/Clp_N pooled; near-universal",
    "cfa": "PF02353 also covers mycolic-acid methyltransferases",
}


def all_families():
    s = set()
    for r in PANEL.values():
        s.update(r.get("any", []))
        s.update(r.get("require_genome", []))
        for g in r.get("all_protein", []):
            s.update(g)
    return sorted(s)


def fetch_hmms(directory: Path):
    directory.mkdir(parents=True, exist_ok=True)
    for pf in all_families():
        f = directory / f"{pf}.hmm"
        if not f.exists():
            f.write_bytes(gzip.decompress(urllib.request.urlopen(
                f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pf}?annotation=hmm", timeout=60).read()))
    return [directory / f"{pf}.hmm" for pf in all_families()]


def apply_rules(protein_families: dict) -> dict:
    """protein_families: {protein_id: set(PFxxxxx)} -> {system: copy number}."""
    genome_fams = set().union(*protein_families.values()) if protein_families else set()
    out = {}
    for sys_, r in PANEL.items():
        if any(g not in genome_fams for g in r.get("require_genome", [])):
            out[sys_] = 0; continue
        n = 0
        for fams in protein_families.values():
            if any(f in fams for f in r.get("any", [])) or any(all(f in fams for f in grp) for grp in r.get("all_protein", [])):
                n += 1
        out[sys_] = n
    return out
