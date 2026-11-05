"""Stage 2: broad genomic signatures.

The halophile plan's acidic excess (D+E)-(K+R) is kept with its sign unassumed: acid-adapted
proteomes are not known to mirror halophiles, so the direction is something to be measured.
IVYWREL is carried as a thermal control because many acidophiles are also thermophiles.
"""
from __future__ import annotations

import math
from collections import Counter
from pathlib import Path

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
MIN_PROTEIN_LENGTH = 20
IVYWREL = set("IVYWREL")
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
      "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
      "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
# EMBOSS pKa set
PKA_POS = {"K": 10.8, "R": 12.5, "H": 6.5}
PKA_NEG = {"D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}
NTERM, CTERM = 8.6, 3.6


def parse_fasta(path):
    seqs, name, buf = [], None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.upper())
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def isoelectric_point(seq: str) -> float:
    counts = Counter(seq)

    def net(ph):
        pos = 1 / (1 + 10 ** (ph - NTERM)) + sum(counts[a] / (1 + 10 ** (ph - p)) for a, p in PKA_POS.items())
        neg = 1 / (1 + 10 ** (CTERM - ph)) + sum(counts[a] / (1 + 10 ** (p - ph)) for a, p in PKA_NEG.items())
        return pos - neg

    lo, hi = 0.0, 14.0
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if net(mid) > 0 else (lo, mid)
    return (lo + hi) / 2


def gravy(seq: str) -> float:
    return float(np.mean([KD[a] for a in seq if a in KD])) if seq else float("nan")


def has_tm_helix(seq: str, window=19, cutoff=1.6) -> bool:
    """Kyte-Doolittle proxy for a membrane protein. Replace with DeepTMHMM before any claim."""
    v = np.array([KD.get(a, 0.0) for a in seq])
    if len(v) < window:
        return False
    return bool((np.convolve(v, np.ones(window) / window, mode="valid") >= cutoff).any())


def proteome_features(proteins: list) -> dict:
    """proteins: list of (id, sequence). Residue-weighted composition plus per-protein pI summaries."""
    seqs = [s.rstrip("*") for _, s in proteins]
    seqs = [s for s in seqs if len(s) >= MIN_PROTEIN_LENGTH and sum(a in AA for a in s) / len(s) >= 0.95]
    if not seqs:
        raise ValueError("no usable proteins")
    total = Counter()
    for s in seqs:
        total.update(a for a in s if a in AA)
    n = sum(total.values())
    f = {a: total[a] / n for a in AA}
    pis = np.array([isoelectric_point(s) for s in seqs])
    gr = np.array([gravy(s) for s in seqs])
    memb = np.array([has_tm_helix(s) for s in seqs])
    out = {f"frac_{a}": f[a] for a in AA}
    out.update({
        "n_proteins": len(seqs),
        "acidic_excess": (f["D"] + f["E"]) - (f["K"] + f["R"]),
        "frac_acidic": f["D"] + f["E"], "frac_basic": f["K"] + f["R"],
        "ivywrel": sum(f[a] for a in IVYWREL),
        "mean_pI": float(pis.mean()), "median_pI": float(np.median(pis)),
        "frac_pI_gt7": float((pis > 7).mean()), "frac_pI_lt5": float((pis < 5).mean()),
        "mean_gravy": float(gr.mean()),
        "frac_membrane_proxy": float(memb.mean()),
        "mean_pI_membrane": float(pis[memb].mean()) if memb.any() else float("nan"),
        "mean_pI_soluble": float(pis[~memb].mean()) if (~memb).any() else float("nan"),
        "frac_Cys": f["C"], "frac_His": f["H"],
    })
    return out


def genome_features(fasta_path) -> dict:
    """Genome size, GC and dinucleotide odds ratios rho_XY = f_XY / (f_X f_Y), strand-symmetrised."""
    seq = "".join(s for _, s in parse_fasta(fasta_path))
    seq = "".join(c for c in seq if c in "ACGT")
    n = len(seq)
    if n < 2:
        raise ValueError("empty genome")
    c = Counter(seq)
    fx = {b: c[b] / n for b in "ACGT"}
    di = Counter(seq[i:i + 2] for i in range(n - 1))
    out = {"genome_size_bp": n, "gc": fx["G"] + fx["C"]}
    for x in "ACGT":
        for y in "ACGT":
            out[f"rho_{x}{y}"] = (di[x + y] / (n - 1)) / (fx[x] * fx[y])
    return out


def gc3_from_cds(cds_fasta) -> float:
    gc = tot = 0
    for _, s in parse_fasta(cds_fasta):
        for i in range(2, len(s) - 2, 3):
            if s[i] in "ACGT":
                tot += 1
                gc += s[i] in "GC"
    return gc / tot if tot else float("nan")


def his_breakdown(proteins: list) -> dict:
    """Where the histidines are: clustered (another His within +-3 residues: HH, HxH, HxxH, HxxxH, the usual metal-binding spacings) vs isolated,
    and in membrane-proxy vs soluble proteins. Fractions are of all residues (clustered/isolated) or of the residues in that protein class."""
    import re
    seqs = [s.rstrip("*") for _, s in proteins]
    seqs = [s for s in seqs if len(s) >= MIN_PROTEIN_LENGTH and sum(a in AA for a in s) / len(s) >= 0.95]
    tot = his_c = his_i = 0
    res_m = his_m = res_s = his_s = n_rich = 0
    for s in seqs:
        L = len(s); h = [i for i, a in enumerate(s) if a == "H"]
        hs = set(h); clustered = sum(1 for i in h if any((i + d) in hs for d in (-3, -2, -1, 1, 2, 3)))
        tot += L; his_c += clustered; his_i += len(h) - clustered
        if has_tm_helix(s):
            res_m += L; his_m += len(h)
        else:
            res_s += L; his_s += len(h)
        n_rich += (len(h) / L) >= 0.05
    return {"his_clustered_frac": his_c / tot, "his_isolated_frac": his_i / tot, "his_frac_membrane": his_m / res_m if res_m else float("nan"),
            "his_frac_soluble": his_s / res_s if res_s else float("nan"), "frac_proteins_his_rich": n_rich / len(seqs), "his_total_check": (his_c + his_i) / tot}


def signal_peptide_like(seq: str) -> bool:
    """Crude Sec-type signal-peptide heuristic (no dedicated predictor is installed): starts with M, a basic residue in the first 6 positions,
    a 9-residue window with mean Kyte-Doolittle >= 1.8 starting between positions 5 and 25, and no further TM-like segment after position 35.
    Approximate: use for class-level comparisons only."""
    if len(seq) < 40 or seq[0] != "M" or not any(a in "KR" for a in seq[:6]):
        return False
    v = np.array([KD.get(a, 0.0) for a in seq[:45]])
    w = np.convolve(v, np.ones(9) / 9, mode="valid")
    if not (w[5:25] >= 1.8).any():
        return False
    return not has_tm_helix(seq[35:])


def compartment_composition(proteins: list, marker_names: set) -> dict:
    """Residue-weighted charge composition by protein class: membrane (TM helix), secreted-like (signal peptide, no TM), cytoplasmic (the rest),
    the conserved core (ribosomal-protein markers, identical protein families in every genome) and the non-core soluble remainder."""
    acc = {k: np.zeros(5) for k in ("membrane", "secreted", "cytoplasmic", "core", "noncore_soluble", "all")}   # residues, D, E, K, R
    counts = dict.fromkeys(acc, 0)
    for name, s in proteins:
        s = s.rstrip("*")
        if len(s) < MIN_PROTEIN_LENGTH or sum(a in AA for a in s) / len(s) < 0.95:
            continue
        vec = np.array([len(s), s.count("D"), s.count("E"), s.count("K"), s.count("R")], float)
        sp = signal_peptide_like(s)                       # tested first: a cleavable signal peptide's hydrophobic stretch must not make the protein "membrane"
        tm = (not sp) and has_tm_helix(s)
        cls = "secreted" if sp else ("membrane" if tm else "cytoplasmic")
        for k in (cls, "all"):
            acc[k] += vec; counts[k] += 1
        if name in marker_names:
            acc["core"] += vec; counts["core"] += 1
        elif not tm and not sp:
            acc["noncore_soluble"] += vec; counts["noncore_soluble"] += 1
    out = {}
    for k, v in acc.items():
        n = max(v[0], 1.0)
        out[f"{k}_acidic_excess"] = (v[1] + v[2] - v[3] - v[4]) / n
        out[f"{k}_DE"] = (v[1] + v[2]) / n
        out[f"{k}_KR"] = (v[3] + v[4]) / n
        out[f"{k}_n_proteins"] = counts[k]
        out[f"{k}_share_residues"] = v[0] / max(acc["all"][0], 1.0)
    return out
