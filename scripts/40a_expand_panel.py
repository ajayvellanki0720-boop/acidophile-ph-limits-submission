"""Build an EXPANDED structural reference panel, data-driven rather than hand-picked -> data/curated/structure_panel_refs.csv (+ models in data/raw/structures).
The panel in 36a was seven membrane proteins chosen by me, so the result could in principle depend on that choice. This selects references systematically:
  1. every reviewed Escherichia coli K-12 protein in UniProt with at least 4 transmembrane helices and at least 3 periplasmic and 3 cytoplasmic topological domains, carrying
     exactly one Pfam family (single-family proteins only, so the alignment to the reference is unambiguous);
  2. one reference per Pfam, the one with the most annotated topological domains;
  3. the Pfam HMM fetched from InterPro and required to hit at least MINPREV of a random sample of bacterial proteomes, so that the same protein really can be compared
     across species - this is what excludes the many narrowly distributed transporters among the candidates;
  4. the AlphaFold model fetched and required to have median pLDDT >= 70 and at least MINPOS exposed periplasmic positions surviving the per-residue pLDDT >= 70 filter,
     since those are the positions the test actually uses.
Soluble references (no topology needed, classed surface against buried) are added from a fixed conserved set as controls.
Usage: python3 scripts/40a_expand_panel.py [sample_genomes] [processes]"""
import gzip
import json
import sys
import urllib.parse
import urllib.request
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import numpy as np
import pandas as pd
import pyhmmer
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1

NSAMP = int(sys.argv[1]) if len(sys.argv) > 1 else 150
PROCS = int(sys.argv[2]) if len(sys.argv) > 2 else 8
MINPREV, MINPOS, MINPLDDT = 0.80, 15, 70.0
S = ROOT / "data/raw/structures"; S.mkdir(parents=True, exist_ok=True)
HMMD = ROOT / "data/raw/hmm_panel"; HMMD.mkdir(parents=True, exist_ok=True)
SOLUBLE = {"EF-Tu": ("P0CE47", "PF00009"), "GroEL": ("P0A6F5", "PF00118"), "RecA": ("P0A7G6", "PF00154"),
           "ATPsynAlpha": ("P0ABB0", "PF00006"), "DnaK": ("P0A6Y8", "PF00012")}
MAXASA = dict(A=129, R=274, N=195, D=193, C=167, Q=225, E=223, G=104, H=224, I=197, L=201, K=236, M=224, F=240, P=159, S=155, T=172, W=285, Y=263, V=174)
OUTSIDE = {"Periplasmic", "Extracellular"}
_s = lambda x: x.decode() if isinstance(x, bytes) else str(x)


def uniprot_candidates():
    q = {"query": "organism_id:83333 AND reviewed:true AND ft_transmem:*", "format": "json", "size": 500,
         "fields": "accession,protein_name,length,ft_transmem,ft_topo_dom,xref_pfam"}
    u = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(q)
    res = []
    while u:
        r = urllib.request.urlopen(u, timeout=120); d = json.load(r); res += d["results"]
        link = r.headers.get("Link", ""); u = link.split("<")[1].split(">")[0] if 'rel="next"' in link else None
    out = []
    for x in res:
        fts = x.get("features", [])
        tm = [f for f in fts if f["type"] == "Transmembrane"]; td = [f for f in fts if f["type"] == "Topological domain"]
        peri = [f for f in td if f.get("description") in OUTSIDE]; cyto = [f for f in td if f.get("description") == "Cytoplasmic"]
        pf = [c["id"] for c in x.get("uniProtKBCrossReferences", []) if c["database"] == "Pfam"]
        if len(tm) >= 4 and len(peri) >= 3 and len(cyto) >= 3 and len(pf) == 1:
            out.append({"acc": x["primaryAccession"], "pfam": pf[0], "n_topo": len(peri) + len(cyto), "n_tm": len(tm),
                        "name": x.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "?")[:60]})
    best = {}
    for c in out:
        if c["pfam"] not in best or c["n_topo"] > best[c["pfam"]]["n_topo"]:
            best[c["pfam"]] = c
    return out, list(best.values())


def get_hmm(pf):
    f = HMMD / f"{pf}.hmm"
    if not f.exists():
        for d in ("data/raw/hmm_markers", "data/raw/hmm"):
            if (ROOT / d / f"{pf}.hmm").exists():
                f.write_bytes((ROOT / d / f"{pf}.hmm").read_bytes()); return f
        try:
            f.write_bytes(gzip.decompress(urllib.request.urlopen(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pf}?annotation=hmm", timeout=90).read()))
        except Exception as e:
            print("   HMM fetch failed", pf, str(e)[:50], flush=True); return None
    return f


HM = None
PFL = None


def init(pfs):
    global HM, PFL
    HM, PFL = [], []
    for pf in pfs:
        p = HMMD / f"{pf}.hmm"
        if p.exists():
            with pyhmmer.plan7.HMMFile(p) as hf:
                HM.append(next(iter(hf))); PFL.append(pf)


def prevalence(a):
    alpha = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(ROOT / "data/raw/proteomes" / f"{a}.faa", digital=True, alphabet=alpha) as sf:
        seqs = sf.read_block()
    hit = set()
    for top in pyhmmer.hmmsearch(HM, seqs, bit_cutoffs="gathering", cpus=1):
        if any(h.included for h in top):
            hit.add(_s(top.query.accession).split(".")[0] if top.query.accession else _s(top.query.name))
    return hit


def model_ok(acc, membrane):
    pdb = S / f"AF-{acc}.pdb"; js = S / f"{acc}.json"
    try:
        if not pdb.exists():
            m = json.load(urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=90))[0]
            pdb.write_bytes(urllib.request.urlopen(m["pdbUrl"], timeout=180).read())
        if not js.exists():
            js.write_bytes(urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=60).read())
    except Exception as e:
        return None, f"model fetch failed: {str(e)[:40]}"
    st = PDBParser(QUIET=True).get_structure(acc, pdb)
    ShrakeRupley().compute(st[0], level="R")
    aa, rsa, pl = [], [], []
    for r in st[0]["A"]:
        a1 = seq1(r.get_resname())
        if a1 not in MAXASA:
            a1 = "G"
        aa.append(a1); rsa.append(r.sasa / MAXASA[a1]); pl.append(next(r.get_atoms()).get_bfactor())
    pl = np.array(pl); rsa = np.array(rsa)
    if np.median(pl) < MINPLDDT:
        return None, f"median pLDDT {np.median(pl):.0f}"
    if not membrane:
        return int(((pl >= MINPLDDT) & (rsa >= 0.25)).sum()), "ok"
    topo = ["?"] * len(aa)
    for f in json.loads(js.read_text())["features"]:
        if f["type"] not in ("Transmembrane", "Topological domain"):
            continue
        s0, s1 = f["location"]["start"]["value"], f["location"]["end"]["value"]
        lab = "TM" if f["type"] == "Transmembrane" else ("peri" if f.get("description") in OUTSIDE else "cyto" if f.get("description") == "Cytoplasmic" else None)
        if lab:
            for i in range(s0 - 1, min(s1, len(aa))): topo[i] = lab
    n = sum(1 for i in range(len(aa)) if topo[i] == "peri" and pl[i] >= MINPLDDT and rsa[i] >= 0.25)
    return (n, "ok") if n >= MINPOS else (n, f"only {n} exposed periplasmic positions")


if __name__ == "__main__":
    allc, uniq = uniprot_candidates()
    print(f"{len(allc)} candidate proteins -> {len(uniq)} distinct single-Pfam families", flush=True)
    got = [c for c in uniq if get_hmm(c["pfam"])]
    print(f"{len(got)} with an HMM available", flush=True)
    PFL_SHARED = [c["pfam"] for c in got]
    accs = [a for a in pd.read_csv(ROOT / "results/broad2/phenotypes_all.csv").ncbi_assembly.dropna().unique()
            if (ROOT / "data/raw/proteomes" / f"{a}.faa").exists()]
    rng = np.random.default_rng(40); samp = list(rng.choice(accs, size=min(NSAMP, len(accs)), replace=False))
    init(PFL_SHARED)
    print(f"testing prevalence of {len(PFL)} families on {len(samp)} genomes", flush=True)
    with Pool(PROCS, initializer=init, initargs=(PFL_SHARED,)) as pl_:
        hits = list(pl_.imap_unordered(prevalence, samp, chunksize=4))
    prev = {pf: np.mean([pf in h for h in hits]) for pf in PFL}
    rows = []
    for c in got:
        pv = prev.get(c["pfam"], 0.0)
        if pv < MINPREV:
            continue
        n, why = model_ok(c["acc"], True)
        rows.append({"label": c["acc"], "uniprot": c["acc"], "pfam": c["pfam"], "membrane": True, "prevalence": pv,
                     "exposed_peri_positions": n, "status": why, "n_tm": c["n_tm"], "name": c["name"]})
        print(f"  {c['acc']} {c['pfam']} prev {pv:.2f} -> {why} | {c['name'][:48]}", flush=True)
    for lab, (acc, pf) in SOLUBLE.items():
        if not get_hmm(pf):
            continue
        n, why = model_ok(acc, False)
        rows.append({"label": lab, "uniprot": acc, "pfam": pf, "membrane": False, "prevalence": np.nan,
                     "exposed_peri_positions": n, "status": why, "n_tm": 0, "name": lab})
    R = pd.DataFrame(rows)
    R["accepted"] = (R.status == "ok")
    R.to_csv(ROOT / "data/curated/structure_panel_refs.csv", index=False)
    print(f"\n{int(R.accepted.sum())} accepted references ({int((R.accepted & R.membrane).sum())} membrane, {int((R.accepted & ~R.membrane).sum())} soluble)")
    print(R[R.accepted][["label", "pfam", "membrane", "prevalence", "exposed_peri_positions", "n_tm", "name"]].to_string(index=False))
    print("\nrejected:", R[~R.accepted][["label", "pfam", "status"]].to_string(index=False))
