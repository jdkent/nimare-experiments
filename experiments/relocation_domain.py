"""Candidate fixes for the relocation null's spatial-density mismatch.

Established: relocating foci uniformly over the analysis mask puts fewer studies on each
covered voxel than the observed foci do (2.27 vs 1.94 at span 20; 3.84 vs 1.94 at span 12).
With "none" the standard error is a pure function of that multiplicity, so the null's |z| runs
low -- it assigns p = 0.038 to the truth's own 0.95 point -- and the test is anticonservative.
This is the uniform-relocation assumption every coordinate-based null makes, made visible by a
statistic that lives far out in the tail.

Three candidates, measured against the same global-null data:

    uniform     what the estimator does now: any in-mask voxel, equally likely.
    hull        relocate only within the bounding box of the observed foci, dilated by the
                coverage radius. Signal-independent in the sense that matters -- it follows
                where the literature looked, not where it found things -- and it is the
                analysis-mask fix (use a gray-matter mask) done automatically.
    density     draw from a heavily smoothed density of all observed foci. Matches the observed
                geometry most closely, at the risk of importing the signal into the null.

Reported: mean studies per covered voxel (the quantity that has to match), the realised
rejection rate under a global null, and -- for hull and density -- the power cost at a true
focal effect, since a null that is harder to beat is only a fix if it still detects something.
"""
import json, sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from pathlib import Path
from scipy import ndimage
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.meta.cbma.effectsize import _p_from_histogram, _null_bin_edges, _NULL_MAX_Z

OUT = Path("/tmp/claude-0/reloc"); OUT.mkdir(parents=True, exist_ok=True)
SHAPE = (12, 12, 12)
AFF = np.diag([4.0, 4.0, 4.0, 1.0]); AFF[:3, 3] = -22.0
MASK = str(OUT / "mask.nii.gz")
nib.save(nib.Nifti1Image(np.ones(SHAPE, np.int32), AFF), MASK)


def build(seed, span=20.0, n_coord=15, effect=None):
    rng = np.random.default_rng(seed)
    studies = []
    for k in range(n_coord):
        n = int(rng.integers(20, 40))
        pts = []
        if effect is not None:
            pts.append({"space": "MNI",
                        "coordinates": [float(v) for v in np.array(effect) + rng.normal(0, 4, 3)],
                        "values": [{"kind": "Z", "value": float(rng.uniform(4.5, 6.5))}]})
        pts += [{"space": "MNI", "coordinates": [float(v) for v in rng.uniform(-span, span, 3)],
                 "values": [{"kind": "Z", "value": float(rng.uniform(3.3, 5.0))}]}
                for _ in range(4)]
        studies.append({"id": f"c{k}", "name": f"c{k}", "metadata": {"sample_sizes": [n]},
            "analyses": [{"id": f"c{k}-1", "name": "1", "metadata": {"sample_sizes": [n]},
                "points": pts, "images": []}]})
    path = OUT / f"ss_{seed}_{span:.0f}_{effect is not None}.json"
    path.write_text(json.dumps({"id": "m", "name": "m", "studies": studies}))
    return Studyset(str(path))


def domain(est, kind):
    """Voxel indices a relocated focus may land on, plus their probabilities."""
    ijk_all = est._in_mask_ijk()
    foci = est._focus_table_[["i", "j", "k"]].values.astype(int)
    if kind == "uniform":
        return ijk_all, None
    if kind == "hull":
        pad = int(np.ceil((est.coverage_radius or 2 * est.fwhm)
                          / min(est.masker.mask_img.header.get_zooms()[:3])))
        lo, hi = foci.min(0) - pad, foci.max(0) + pad
        keep = np.all((ijk_all >= lo) & (ijk_all <= hi), axis=1)
        return ijk_all[keep], None
    if kind == "density":
        vol = np.zeros(SHAPE, dtype=float)
        np.add.at(vol, tuple(foci.T), 1.0)
        # smoothing far wider than the pooling kernel: follows where the literature looked,
        # not the peaks it found
        vol = ndimage.gaussian_filter(vol, sigma=3.0) + 1e-3
        p = vol[tuple(ijk_all.T)]
        return ijk_all, p / p.sum()
    raise ValueError(kind)


def null_histogram(est, kind, n_iters, seed):
    ijk, p = domain(est, kind)
    rng = np.random.default_rng(seed)
    hist = np.zeros(len(_null_bin_edges()) - 1, dtype=float)
    mult = []
    for _ in range(n_iters):
        t = est._focus_table_.copy()
        idx = rng.choice(len(ijk), size=len(t), p=p) if p is not None \
            else rng.integers(0, len(ijk), size=len(t))
        t[["i", "j", "k"]] = ijk[idx]
        fit, z = est._statistic(t, est._sample_sizes_, est._thresholds_, est._image_studies_)
        mult.append(float(fit["n_studies"][fit["covered"]].mean()))
        counts, _ = np.histogram(np.clip(np.abs(z), 0, _NULL_MAX_Z), bins=_null_bin_edges())
        hist += counts
    return hist, float(np.mean(mult))


S = int(sys.argv[1]) if len(sys.argv) > 1 else 6
ITERS = int(sys.argv[2]) if len(sys.argv) > 2 else 100
SPAN = 12.0
KINDS = ("uniform", "hull", "density")
print(f"{S} datasets per cell, {ITERS} relocations each, foci spanning +/-{SPAN:.0f} mm of a "
      f"+/-22 mm mask.")
print("n_obs / n_null are mean studies per covered voxel; they should agree.\n", flush=True)
print(f"{'selection':>14s} {'domain':>9s} {'n_obs':>6s} {'n_null':>7s} {'FPR':>7s} "
      f"{'power':>7s}", flush=True)

for sel in ("none", "zero-inflated"):
    acc = {k: {"fpr": [], "power": [], "nobs": [], "nnull": []} for k in KINDS}
    for seed in range(S):
        for effect in (None, (0.0, 0.0, 0.0)):
            est = CBES(fwhm=10.0, mask=MASK, null_method="none", selection_model=sel,
                       threshold="study-min", peak_bias="per-study", seed=seed)
            est.fit(build(seed, span=SPAN, effect=effect))
            fit, z = est._statistic(est._focus_table_, est._sample_sizes_,
                                    est._thresholds_, est._image_studies_)
            centre = np.argmin(np.abs(est._in_mask_ijk() - np.array([5, 5, 5])).sum(1))
            for kind in KINDS:
                hist, m = null_histogram(est, kind, ITERS, 5000 + seed)
                p = _p_from_histogram(np.abs(z), hist)
                if effect is None:
                    acc[kind]["nobs"].append(float(fit["n_studies"][fit["covered"]].mean()))
                    acc[kind]["nnull"].append(m)
                    acc[kind]["fpr"].append(float(np.mean(p < 0.05)))
                else:
                    acc[kind]["power"].append(bool(p[centre] < 0.05 / p.size))
    for kind in KINDS:
        a = acc[kind]
        print(f"{sel:>14s} {kind:>9s} {np.mean(a['nobs']):6.2f} {np.mean(a['nnull']):7.2f} "
              f"{np.mean(a['fpr']):7.4f} {np.mean(a['power']):7.2f}", flush=True)
    print(flush=True)
