"""Search each proteome with the Pfam panel (gathering cutoffs) and compare with the keyword calls.
Usage: python3 scripts/05_hmm_systems.py phenotypes.csv PROTEOME_DIR keyword_counts.csv|none OUTDIR
Writes OUTDIR/hmm_counts.csv (hm_cn_*, hm_rate_*, hm_pw_*; rows already present are reused) and, if keyword counts are given, OUTDIR/detection_validation.csv."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import pandas as pd
import pyhmmer
from acidophile import hmm_panel as hp

_s = lambda x: x.decode() if isinstance(x, bytes) else str(x)
ph, prot_dir, kw_csv, out = sys.argv[1:5]
out = Path(out); out.mkdir(parents=True, exist_ok=True)
root = Path(__file__).resolve().parents[1]
hmms = []
for f in hp.fetch_hmms(root / "data/raw/hmm"):
    with pyhmmer.plan7.HMMFile(f) as hf:
        hmms.extend(hf)
alpha = pyhmmer.easel.Alphabet.amino()
prev = pd.read_csv(out / "hmm_counts.csv") if (out / "hmm_counts.csv").exists() else pd.DataFrame()
have = set(prev.genome_id) if len(prev) else set()
rows = []
for acc in pd.read_csv(ph)["ncbi_assembly"]:
    if acc in have or not (Path(prot_dir) / f"{acc}.faa").exists():
        continue
    with pyhmmer.easel.SequenceFile(Path(prot_dir) / f"{acc}.faa", digital=True, alphabet=alpha) as sf:
        seqs = sf.read_block()
    pf = {}
    for top in pyhmmer.hmmsearch(hmms, seqs, bit_cutoffs="gathering", cpus=4):
        fam = _s(top.query.accession).split(".")[0]
        for h in top:
            if h.included:
                pf.setdefault(_s(h.name), set()).add(fam)
    cn = hp.apply_rules(pf)
    rows.append({"genome_id": acc, "n_genes": len(seqs), **{f"hm_cn_{k}": v for k, v in cn.items()},
                 **{f"hm_rate_{k}": 1000 * v / len(seqs) for k, v in cn.items()}})
    print(acc, "done", flush=True)
hm = pd.concat([prev, pd.DataFrame(rows)], ignore_index=True) if len(prev) else pd.DataFrame(rows)
cats = {}
from acidophile.systems import SYSTEMS
for k, v in hp.KW_TO_HMM.items():
    cats.setdefault(SYSTEMS[k][0], []).append(v)
for cat, members in cats.items():
    hm[f"hm_pw_{cat}"] = (hm[[f"hm_cn_{m}" for m in members]] > 0).mean(axis=1)
hm.to_csv(out / "hmm_counts.csv", index=False)

if kw_csv.lower() == "none":
    sys.exit(0)
kw = pd.read_csv(kw_csv).set_index("genome_id"); hm = hm.set_index("genome_id").loc[kw.index]
val = []
for k, h in hp.KW_TO_HMM.items():
    a, b = kw[f"cn_{k}"] > 0, hm[f"hm_cn_{h}"] > 0
    tp, fp, fn, tn = int((a & b).sum()), int((a & ~b).sum()), int((~a & b).sum()), int((~a & ~b).sum())
    val.append({"system": k, "hmm_prev": b.mean(), "kw_prev": a.mean(), "kw_true_pos": tp, "kw_false_pos": fp, "kw_false_neg": fn,
                "kw_precision": tp / (tp + fp) if tp + fp else None, "kw_recall": tp / (tp + fn) if tp + fn else None,
                "copy_corr": kw[f"cn_{k}"].corr(hm[f"hm_cn_{h}"]), "caveat": hp.CAVEATS.get(h, "")})
pd.DataFrame(val).to_csv(out / "detection_validation.csv", index=False)
print(pd.DataFrame(val).round(2).drop(columns="caveat").to_string(index=False))
