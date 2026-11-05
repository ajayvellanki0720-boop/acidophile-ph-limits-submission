import numpy as np
import pandas as pd
import pytest
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from acidophile.lineage import contrast_features, matched_pairs
from acidophile.models import compare_models, make_splitter
from acidophile.phenotype import (REQUIRED_COLUMNS, censored_rows, classify_acid_class, primary_set,
                                  undissociated_fraction, validate)
from acidophile.proteome import isoelectric_point, proteome_features
from acidophile.synthetic import make_synthetic
from acidophile.systems import count_systems, match_system, pathway_scores


def _row(**kw):
    base = dict.fromkeys(REQUIRED_COLUMNS, np.nan)
    base.update(strain_id="A", ncbi_assembly="GCA_000000001.1", ph_min=1.0, ph_min_bound="measured",
                ph_opt=2.0, ph_max=4.0, ph_max_bound="measured", acid_identity="H2SO4",
                ph_measured_at_c=25, confidence="high", domain="Archaea", family="F")
    base.update(kw)
    return base


def test_pi_tracks_biopython_on_realistic_proteins():
    from scipy.stats import spearmanr
    rng = np.random.default_rng(0)
    aa = list("ACDEFGHIKLMNPQRSTVWY")
    w = np.array([8.5, 1.1, 5.5, 6.2, 3.9, 7.4, 2.2, 5.9, 5.0, 10.0, 2.4, 4.0, 4.4, 4.4, 5.5, 5.8, 5.4, 7.1, 1.3, 2.9])
    seqs = []
    for i in range(60):
        pw = w * rng.uniform(0.4, 1.6, 20)  # vary composition so pI spans acidic to basic
        seqs.append("".join(rng.choice(aa, 250, p=pw / pw.sum())))
    ours = [isoelectric_point(s) for s in seqs]
    ref = [ProteinAnalysis(s).isoelectric_point() for s in seqs]
    assert spearmanr(ours, ref)[0] > 0.95 and np.mean(np.abs(np.array(ours) - ref)) < 0.5


def test_pi_direction():
    assert isoelectric_point("D" * 30) < 4 < 10 < isoelectric_point("K" * 30)


def test_proteome_features_acidic_excess():
    f = proteome_features([("p1", "DEDEDEDEDEDEDEDEDEDEKA"), ("p2", "ACDEFGHIKLMNPQRSTVWYAC")])
    assert f["acidic_excess"] > 0 and 0 < f["mean_pI"] < 14


def test_validate_flags_and_never_repairs():
    df = pd.DataFrame([_row(), _row(strain_id="B", ph_min=5, ph_opt=2),
                       _row(strain_id="C", ph_min_bound="grew_lower"), _row(strain_id="D", ncbi_assembly=np.nan)])
    out, issues = validate(df)
    assert out.loc[1, "ph_min"] == 5  # not silently fixed
    probs = set(zip(issues.strain_id, issues.field))
    assert ("B", "ph_min/ph_opt/ph_max") in probs and ("C", "ph_min_bound") in probs and ("D", "ncbi_assembly") in probs


def test_censored_excluded_from_primary():
    df = pd.DataFrame([_row(), _row(strain_id="B", ph_min_bound="lowest_tested")])
    out, _ = validate(df)
    assert list(primary_set(out).strain_id) == ["A"] and list(censored_rows(out).strain_id) == ["B"]


def test_weak_acid_chemistry():
    assert undissociated_fraction(4.76, 4.76) == pytest.approx(0.5)
    assert undissociated_fraction(2.0, 4.76) > 0.99 > undissociated_fraction(7.0, 4.76)


def test_acid_class():
    assert classify_acid_class(1.5, 4, 0) == "extreme_acidophile"
    assert classify_acid_class(7.0, 9, 4.0) == "acid_tolerant"
    assert classify_acid_class(7.0, 9, 6.0) == "non_acidophile"


def test_system_matching_and_exclusions():
    assert match_system("kdpA", "") == ["kdp"]
    assert "cad" not in match_system("cadA", "cadmium-translocating P-type ATPase")
    assert "f_atp_synthase" not in match_system("", "V-type ATP synthase subunit A")
    assert "v_a_atp_synthase" in match_system("", "V-type ATP synthase subunit A")


def test_count_and_pathway():
    ann = pd.DataFrame({"genome_id": ["g1"] * 4 + ["g2"], "protein_id": list("abcde"),
                        "gene": ["kdpA", "kdpB", "gadB", "", ""], "product": ["", "", "", "hypothetical", "hypothetical"]})
    c = count_systems(ann)
    assert c.set_index("genome_id").loc["g1", "cn_kdp"] == 2
    pw = pathway_scores(c).set_index("genome_id")
    assert pw.loc["g1", "pw_donnan_k"] > 0 and pw.loc["g2", "pw_donnan_k"] == 0


def test_blocked_splits_never_share_family():
    df = make_synthetic()
    for scheme in ("family_blocked", "leave_one_lineage_out"):
        g = df["family"] if scheme == "family_blocked" else df["class"]
        for tr, te in make_splitter(g, scheme).split(df, groups=g):
            assert not set(g.iloc[tr]) & set(g.iloc[te])


def test_planted_signal_recovered_and_null_not_claimed():
    res, _ = compare_models(make_synthetic(signal=True))
    assert res.set_index("model").loc["M4_gene_content", "beats_baseline"]
    null, _ = compare_models(make_synthetic(signal=False, seed=1))
    assert not null.set_index("model").loc["M2_acidic_proteome", "beats_baseline"]


def test_too_few_groups_refused():
    with pytest.raises(ValueError):
        compare_models(make_synthetic(n_families=2))


def test_matched_pairs_and_contrast():
    df = make_synthetic(per_family=8)
    p = matched_pairs(df, rank="genus", opt_tol=0.8, min_min_diff=0.8)
    assert (p.d_ph_min < 0).all()
    assert contrast_features(df, p, ["rate_kdp"], min_pairs=1).shape[0] <= 1


def test_hmm_rules_conjunction_and_genome_requirement():
    from acidophile.hmm_panel import apply_rules
    pf = {"p1": {"PF03814"}, "p2": {"PF00122", "PF00403"}, "p3": {"PF00122"}}
    r = apply_rules(pf)
    assert r["kdp"] == 0                     # KdpA without KdpC anywhere in the genome is not a Kdp system
    assert r["p_type_atpase"] == 1           # only the protein with both E1-E2 and HMA counts
    pf["p4"] = {"PF02669"}
    assert apply_rules(pf)["kdp"] == 1


def _lineage(n=16, effect=0.0, seed=3):
    rng = np.random.default_rng(seed)
    feat = rng.normal(size=n)
    opt = rng.uniform(1, 5, n)
    return pd.DataFrame({"ph_opt": opt, "ph_min": opt - 1.0 - effect * feat + rng.normal(0, 0.15, n),
                         "mean_pI": feat, "genome_size_bp": rng.normal(size=n), "ivywrel": rng.normal(size=n)})


def test_within_lineage_detects_planted_signal_and_not_null():
    from acidophile.within_lineage import analyse
    tab, s = analyse(_lineage(effect=0.8), n_perm=1000)
    assert tab.iloc[0].feature == "mean_pI" and tab.iloc[0].p_family_wise < 0.05
    _, s0 = analyse(_lineage(effect=0.0, seed=5), n_perm=1000)
    assert s0["p_family_wise_max"] > 0.05


def test_within_lineage_reports_detectable_effect_shrinking_with_n():
    from acidophile.within_lineage import analyse
    _, small = analyse(_lineage(n=8), n_perm=600)
    _, large = analyse(_lineage(n=40), n_perm=600)
    assert small["detectable_abs_rho_family_wise"] > large["detectable_abs_rho_family_wise"]


def _balanced_tree(depth=4):
    from io import StringIO
    from Bio import Phylo

    def build(d, prefix=""):
        return f"{prefix}:1.0" if d == 0 else "(" + build(d - 1, prefix + "a") + "," + build(d - 1, prefix + "b") + "):1.0"
    t = Phylo.read(StringIO(build(depth) + ";"), "newick")
    for i, c in enumerate(t.get_terminals()):
        c.name = f"t{i}"
    return t


def test_vcv_is_symmetric_psd_with_sister_covariance():
    from acidophile.phylo import vcv_from_tree
    t = _balanced_tree(3)
    tips = [c.name for c in t.get_terminals()]
    C = vcv_from_tree(t, tips)
    assert np.allclose(C, C.T) and np.linalg.eigvalsh(C).min() > -1e-9
    assert C[0, 0] > C[0, 1] > C[0, 2] > C[0, 4] == 0      # self > sister > cousin > across the root


def test_pgls_recovers_planted_effect_and_is_calibrated_under_null():
    from acidophile.phylo import pgls, vcv_from_tree
    t = _balanced_tree(4)
    tips = [c.name for c in t.get_terminals()]
    C = vcv_from_tree(t, tips); L = np.linalg.cholesky(C + 1e-9 * np.eye(len(tips)))
    rng = np.random.default_rng(0)
    hits, betas = 0, []
    for _ in range(300):
        x = L @ rng.normal(size=len(tips))
        y_null = L @ rng.normal(size=len(tips))
        hits += pgls(y_null, np.column_stack([np.ones(len(tips)), x]), C)["p"][1] < 0.05
        y = 0.8 * x + L @ rng.normal(size=len(tips)) * 0.3
        betas.append(pgls(y, np.column_stack([np.ones(len(tips)), x]), C)["beta"][1])
    assert hits / 300 < 0.12                     # nominal 0.05; ML lambda makes it slightly liberal at n=16
    assert abs(np.mean(betas) - 0.8) < 0.1


def test_nj_groups_sister_sequences():
    from acidophile.phylo import distance_matrix, nj_tree
    M = np.array([list("AAAAAAAAAACCCCCCCCCC" * 3), list("AAAAAAAAAGCCCCCCCCCC" * 3),
                  list("TTTTTTTTTTGGGGGGGGGG" * 3), list("TTTTTTTTTAGGGGGGGGGG" * 3)])
    t = nj_tree(distance_matrix(M), ["A", "B", "C", "D"])
    clades = {frozenset(c.name for c in cl.get_terminals()) for cl in t.get_nonterminals()}
    assert frozenset({"A", "B"}) in clades or frozenset({"C", "D"}) in clades


def test_vcv_fast_matches_reference():
    from acidophile.phylo import vcv_fast, vcv_from_tree
    t = _balanced_tree(4)
    tips = [c.name for c in t.get_terminals()]
    assert np.allclose(vcv_fast(t, tips), vcv_from_tree(t, tips))


def test_fast_distance_and_nj_match_reference():
    from acidophile.phylo import distance_matrix, distance_matrix_fast, nj_tree, nj_tree_fast
    rng = np.random.default_rng(4)
    base = rng.integers(0, 20, 400)
    rows = []
    for k in range(14):
        r = base.copy(); idx = rng.choice(400, 20 + 6 * k, replace=False); r[idx] = rng.integers(0, 20, len(idx))
        rows.append(["ACDEFGHIKLMNPQRSTVWY"[v] for v in r])
    M = np.array(rows); M[3, :40] = "-"; M[9, 100:160] = "-"
    names = [f"t{i}" for i in range(14)]
    assert np.allclose(distance_matrix(M), distance_matrix_fast(M), atol=1e-9)

    def splits(t):
        tips = frozenset(x.name for x in t.get_terminals()); out = set()
        for c in t.get_nonterminals():
            s = frozenset(x.name for x in c.get_terminals())
            if 1 < len(s) < len(tips) - 1:
                out.add(min(s, tips - s, key=lambda x: (len(x), sorted(x))))
        return out
    assert splits(nj_tree(distance_matrix(M), names)) == splits(nj_tree_fast(distance_matrix_fast(M), names))


def test_pgls_scan_matches_pgls_at_fixed_lambda():
    from acidophile.phylo import pgls, pgls_scan, vcv_from_tree
    t = _balanced_tree(5)
    tips = [c.name for c in t.get_terminals()]
    C = vcv_from_tree(t, tips); rng = np.random.default_rng(8); n = len(tips)
    x, f = rng.normal(size=n), rng.normal(size=n)
    y = 1 + 0.5 * x + 0.3 * f + rng.normal(size=n) * 0.5
    Xb = np.column_stack([np.ones(n), x])
    b, se, p, lam = pgls_scan(C, y, Xb, f[:, None], lam=0.4)
    ref = pgls(y, np.column_stack([Xb, (f - f.mean()) / f.std()]), C, lam=0.4)
    assert abs(b[0] - ref["beta"][-1]) < 1e-8 and abs(se[0] - ref["se"][-1]) < 1e-8 and abs(p[0] - ref["p"][-1]) < 1e-8


def test_his_breakdown_partitions_histidines():
    from acidophile.proteome import his_breakdown
    seqs = [("a", "M" + "AHAAHAAAHHAAAAAAAAAAAAAAAAH" + "K" * 10), ("b", "A" * 40 + "H" + "A" * 40)]
    r = his_breakdown(seqs)
    assert abs(r["his_clustered_frac"] + r["his_isolated_frac"] - r["his_total_check"]) < 1e-12
    assert r["his_isolated_frac"] > 0 and r["his_clustered_frac"] > r["his_isolated_frac"]


def test_pic_contrasts_known_values_and_cherries():
    from io import StringIO
    from Bio import Phylo
    from acidophile.phylo import pic_contrasts, regress_through_origin
    t = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
    x = np.array([[1.0], [3.0], [4.0], [10.0]])
    c, ch = pic_contrasts(t, ["A", "B", "C", "D"], x)
    assert len(c) == 3 and ch.sum() == 2
    assert sorted(np.round(np.abs(c[ch][:, 0]), 6)) == sorted(np.round([2 / np.sqrt(2), 6 / np.sqrt(2)], 6))   # (x1-x2)/sqrt(l1+l2)
    # contrasts of a trait regressed on itself give slope 1
    b, se, p = regress_through_origin(c[:, 0], c[:, [0]])
    assert abs(b[0] - 1) < 1e-9


def test_pgls_scan_returns_per_sd_and_not_per_unit():
    """The identity five separate scale errors in this project violated, pinned here rather than only inside scripts/57.

    pgls_scan standardises its tested column internally, so its beta is already per standard deviation of that column.
    The per-unit slope is therefore beta / sd, and multiplying beta by sd converts it back to per unit. The test uses a
    predictor whose SD is far from 1 so that the two scales cannot be confused, and asserts both that the identity holds
    and that the returned value is NOT the per-unit slope - without the second assertion the test would pass under the
    very mistake it exists to catch.
    """
    from acidophile.phylo import pgls, pgls_scan, vcv_from_tree
    t = _balanced_tree(6)
    tips = [c.name for c in t.get_terminals()]
    C = vcv_from_tree(t, tips)
    rng = np.random.default_rng(17)
    n = len(tips)
    sd_f = 7.3                                  # deliberately far from 1
    f = rng.normal(scale=sd_f, size=n)
    y = 2 + 0.04 * f + rng.normal(size=n) * 0.3
    Xb = np.column_stack([np.ones(n)])
    beta, se, p, lam = pgls_scan(C, y, Xb, f[:, None], lam=0.4)
    per_unit = pgls(y, np.column_stack([Xb, f]), C, lam=0.4)["beta"][-1]
    assert abs(beta[0] - per_unit * f.std()) < 1e-8          # the identity
    assert abs(beta[0] / f.std() - per_unit) < 1e-8          # dividing, not multiplying, recovers per unit
    # and the two scales differ by the predictor's realised SD (7.7 here, not the nominal 7.3), so swapping them could not
    # round into agreement: a per-unit return would fail the identity above by a factor of nearly 8
    assert abs(beta[0] / per_unit - f.std()) < 1e-6 and beta[0] / per_unit > 5

    # CORRECTION 15's direction, which the assertions above do not cover. Corrections 2, 4 and 6 all multiplied a
    # per-SD beta by an SD and reported per unit. The opposite mistake is to leave the per-SD beta alone and CALL it
    # per unit, which is what seven scripts did by storing it in a column named `beta_per_pt`: the value is never
    # touched, so no arithmetic identity is violated and nothing fails. Only the size of the discrepancy reveals it.
    # Pinned here as the factor a "per point" claim must be divided by, with the predictor on the scale those scripts
    # actually use - a percentage share, whose SD is a few points, so the mislabelling is a several-fold error and
    # never a rounding one.
    share = rng.uniform(1.0, 12.0, size=n)                   # a percentage share, as scripts/36c and 40c pass it
    ys = 1 + 0.009 * share + rng.normal(size=n) * 0.2
    b_scan, *_ = pgls_scan(C, ys, Xb, share[:, None], lam=0.4)
    b_point = pgls(ys, np.column_stack([Xb, share]), C, lam=0.4)["beta"][-1]
    assert abs(b_scan[0] / share.std() - b_point) < 1e-8, "per point is the scan beta divided by the share's SD"
    assert abs(b_scan[0]) > 2.0 * abs(b_point), "mislabelling a per-SD beta as per-point is a several-fold error"


def test_optimum_as_range_detector_reproduces_the_hand_read_verdicts():
    """The detector behind scripts/61 against the fifteen statements read by hand, with their recorded verdicts.

    It must score every one: eight optimum-read-as-range statements flagged, seven correct statements left alone. The
    first draft of the detector lacked the width condition and flagged three of the correct ones, which is what this
    guards against.
    """
    from acidophile.parser_flags import flag
    cases = [
        ("Strain KSL-104(T) grew optimally at pH 7.0-8.0 and 30 degrees C. It was characterized chemotaxonomically as ha", True),
        ("Strain ISL-14(T) grew optimally at pH 7.0-8.0, 30 degrees C and in the presence of approximately 2 % (w/v) NaC", True),
        ("Strain YIM 48875(T) grew optimally at 28 °C, at pH 6.0-8.0 and in the presence of 2 % (w/v) NaCl", True),
        ("Strain WL0062T grew optimally at 28 °C, pH 7.0-8.0 and with 1.0-3.0% (w/v) NaCl", True),
        ("Strain W-BA2(T) grew optimally at 25 °C, at pH 7.0-8.0 and in the presence of approximately 2.0-3.0% (w/v) NaC", True),
        ("Strain BS-12M(T) grew optimally at 35 °C, at pH 7.0-8.0 and in the presence of 2.0% (w/v) NaCl", True),
        ("Strain 2T4P-2-4T grew optimally without NaCl at 28-30 °C, pH 7.0-8.0. The 16S rRNA gene sequence of strain", True),
        ("Growth optimum was observed to be at 28-37 °C, pH 5.5-6.5 and 2-5% (w/v) NaCl", True),
        ("e was aerobic and exhibited optimal growth from 25-40 °C at pH 6.0-10.0 with 0-1% (w/v) NaCl", False),
        ("een 15 and 35 degrees C (optimum, 30 degrees C) and between pH 6.0 and 9.0 (optimum, pH 7.5-8.5)", False),
        ("ew anaerobically at 20-45 °C (optimally at 37-40 °C) and at pH 6.0-10.0 (optimally at pH 6.0-8.0)", False),
        ("Growth was observed at 40-65 °C (optimum at 60 °C), pH 6-9 (optimum at 6.5-8) and 0.5-5% (optimum at 1-2%) (w/v) NaCl", False),
        ("cterial growth was observed at 10-37 °C (optimum 25 °C) and pH 6.0-8.0 (optimum pH 6.0)", False),
        ("Growth occurred at pH 5.5-8.5 (optimum pH 7.0) and at 10-45 °C (optimum 30 °C)", False),
        ("Strain DHOB07T grew at 10-37 °C, pH 4-7 and 0-0.5 % (w/v) NaCl, with an optimum at 28 °C, pH 5-5.5 and 0% (w/v) NaCl on R2A medium", False),
    ]
    wrong = [(s, want, flag(s)) for s, want in cases if flag(s) != want]
    assert not wrong, f"detector disagrees with the hand-read verdicts on {len(wrong)}: {[w[0][:60] for w in wrong]}"
    assert flag(None) is False and flag("Grown in rich medium at 30 degrees C") is False
