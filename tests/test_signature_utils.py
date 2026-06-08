"""
Synthetic smoke test for signature_utils (precision-modeling requirement).

Proves, without any private data, that:
  - SBS96 channels + SNV classification are correct (incl. pyrimidine folding),
  - the reference loader reorders/normalises correctly,
  - the refit recovers known exposures on simulated mixtures,
  - bootstrap separates present vs absent signatures,
  - the audit emits ok / insufficient_evidence / poor_reconstruction / confounded
    in the situations each label is meant to describe.

Run:  python3 tests/test_signature_utils.py    (needs only numpy)
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin", "utils"))
import signature_utils as su


def _make_block_reference(k=4, names=None, support=96):
    """k distinct signatures, each concentrated on its own block of channels."""
    names = names or [f"SBS_test{i+1}" for i in range(k)]
    W = np.full((96, k), 0.001)
    block = support // k
    for j in range(k):
        lo, hi = j * block, (j + 1) * block
        W[lo:hi, j] += 1.0
    W = W / W.sum(axis=0)
    return names, W


def test_channels_and_classify():
    ch = su.sbs96_channels()
    assert len(ch) == 96 and len(set(ch)) == 96, "expected 96 unique channels"
    assert ch[0] == "A[C>A]A" and ch[-1] == "T[T>G]T"
    assert su.classify_snv("C", "A", "ACA") == "A[C>A]A"
    # purine ref folds to pyrimidine: A>G in TAC  ->  G[T>C]A
    assert su.classify_snv("A", "G", "TAC") == "G[T>C]A"
    # guards
    assert su.classify_snv("C", "C", "ACA") is None        # not a substitution
    assert su.classify_snv("CC", "A", "ACA") is None        # not an SNV
    assert su.classify_snv("C", "A", "AGA") is None        # centre != ref
    assert su.classify_snv("C", "A", "ANA") is None        # ambiguous base
    print("ok  channels + classify_snv")


def test_reference_loader_reorders_and_normalises():
    names, W = _make_block_reference(k=2)
    channels = su.sbs96_channels()
    order = list(range(96))
    rng = np.random.default_rng(0)
    rng.shuffle(order)  # write rows out of canonical order on purpose
    fd, path = tempfile.mkstemp(suffix=".tsv")
    with os.fdopen(fd, "w") as fh:
        fh.write("Type\t" + "\t".join(names) + "\n")
        for i in order:
            fh.write(channels[i] + "\t" + "\t".join(f"{W[i,j]:.6f}" for j in range(2)) + "\n")
    got_names, got_W = su.load_signature_reference(path)
    os.remove(path)
    assert got_names == names
    assert got_W.shape == (96, 2)
    assert np.allclose(got_W.sum(axis=0), 1.0), "columns must normalise to 1"
    assert np.allclose(got_W, W, atol=1e-5), "rows must be back in canonical order"
    print("ok  reference loader")


def test_refit_recovers_known_exposures():
    names, W = _make_block_reference(k=4)
    true_frac = np.array([0.5, 0.3, 0.2, 0.0])
    total = 4000
    expected_spectrum = W @ (true_frac * total)
    counts = np.random.default_rng(1).poisson(expected_spectrum).astype(float)

    h = su.refit_exposures(counts, W, seed=0)
    frac = h / h.sum()
    assert su.cosine_similarity(frac, true_frac) > 0.95, frac
    assert abs(h.sum() - counts.sum()) < 1e-6, "exposures must conserve total mass"
    recon = su.reconstruct(h, W)
    assert su.cosine_similarity(counts, recon) > 0.97
    print(f"ok  refit recovers exposures (frac={np.round(frac,3)})")


def test_bootstrap_separates_present_from_absent():
    names, W = _make_block_reference(k=4)
    true_frac = np.array([0.6, 0.4, 0.0, 0.0])
    counts = np.random.default_rng(2).poisson(W @ (true_frac * 5000)).astype(float)
    boot = su.bootstrap_exposures(counts, W, n_boot=80, seed=0)
    assert boot["lo"][0] > 0.05 and boot["lo"][1] > 0.02, "present sigs must be stable"
    assert boot["lo"][2] < 0.02 and boot["lo"][3] < 0.02, "absent sigs must be unstable"
    print("ok  bootstrap separates present/absent")


def test_audit_labels():
    names, W = _make_block_reference(k=4)

    # 1) clean, high-burden, well-reconstructed, non-artifact -> ok
    counts = np.random.default_rng(3).poisson(W @ np.array([0.5,0.3,0.2,0.0]) * 3000).astype(float)
    h = su.refit_exposures(counts, W, seed=0)
    boot = su.bootstrap_exposures(counts, W, n_boot=60, seed=0)
    a = su.audit_fit(counts, W, names, h, boot, tumor_only=True)
    assert a["label"] == "ok", a["label"]
    assert a["n_stable_signatures"] >= 2

    # 2) too few mutations -> insufficient_evidence
    low = np.random.default_rng(4).poisson(W @ np.array([0.5,0.5,0,0]) * 20).astype(float)
    hl = su.refit_exposures(low, W, seed=0)
    bl = su.bootstrap_exposures(low, W, n_boot=40, seed=0)
    al = su.audit_fit(low, W, names, hl, bl, min_snv=50, tumor_only=True)
    assert al["label"] == "insufficient_evidence", al["label"]

    # 3) reference cannot represent the spectrum -> poor_reconstruction
    names2, W2 = _make_block_reference(k=2, support=48)   # signatures cover ch 0..47 only
    spec = np.zeros(96); spec[48:] = 100.0                # all signal in uncovered channels
    h2 = su.refit_exposures(spec, W2, seed=0)
    b2 = su.bootstrap_exposures(spec, W2, n_boot=30, seed=0)
    a2 = su.audit_fit(spec, W2, names2, h2, b2, min_snv=50, tumor_only=True)
    assert a2["label"] == "poor_reconstruction", (a2["label"], a2["reconstruction_cosine"])

    # 4) artifact-dominated, well reconstructed, high burden -> confounded
    names3 = ["SBS27", "SBS_test2", "SBS_test3", "SBS_test4"]   # SBS27 is an artefact sig
    _, W3 = _make_block_reference(k=4, names=names3)
    counts3 = np.random.default_rng(5).poisson(W3 @ np.array([0.8,0.2,0,0]) * 3000).astype(float)
    h3 = su.refit_exposures(counts3, W3, seed=0)
    b3 = su.bootstrap_exposures(counts3, W3, n_boot=40, seed=0)
    a3 = su.audit_fit(counts3, W3, names3, h3, b3, tumor_only=True)
    assert a3["label"] == "confounded", (a3["label"], a3["artifact_exposure_fraction"])

    print("ok  audit labels (ok / insufficient / poor_recon / confounded)")


if __name__ == "__main__":
    test_channels_and_classify()
    test_reference_loader_reorders_and_normalises()
    test_refit_recovers_known_exposures()
    test_bootstrap_separates_present_from_absent()
    test_audit_labels()
    print("\nALL SIGNATURE-UTIL TESTS PASSED")
