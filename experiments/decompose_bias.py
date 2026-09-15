"""Where does the coordinates-only bias in `g` actually come from?

The coverage run puts it at +0.506 on a truth of 0.800 -- a factor of 1.63. Everything in the
PR attributes that to peak-height selection: a table reports the maximum of a cluster, which
exceeds the truth at the cluster's centre. But a back-of-envelope estimate of the winner's curse
at this signal (true z about 4.4, a handful of effectively independent points per cluster) gives
a reported z near 5.2 and so a g near 0.95, not 1.31. About a third of the bias is unaccounted
for, and guessing which stage adds it has been a reliable way to be wrong.

So the bias is decomposed by measuring each stage on the same collections:

  1. truth at the read-out voxel                      what we are trying to recover
  2. mean reported g among foci landing near it       what the tables actually say
  3. pooled estimate with selection_model='none'      weighting and the kernel, no selection
  4. pooled estimate with the selection model on      the shipped path

Stage 2 minus 1 is the winner's curse plus the localisation offset -- the part no estimator can
see. Stage 3 minus 2 is what the inverse-variance weighting and the 10 mm kernel do. Stage 4
minus 3 is what the zero-inflated censored likelihood does, and its sign is the interesting
part: a censoring correction should pull the estimate *down* toward the unobserved mass below
the threshold, so if it pushes up, something is the wrong way round.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import ndimage
from joblib import Parallel, delayed
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
import reporting

SHAPE, VOXEL_MM, SMOOTH_VOX, RADIUS_VOX = (30, 30, 30), 4.0, 0.8, 2.5
PEAK_G, N_SIMS, N_STUDIES = 0.8, 40, 12
NEAR_MM = 10.0                      # "lands near the read-out voxel", the estimator's kernel
AFF = np.eye(4); AFF[:3, :3] *= VOXEL_MM
AFF[:3, 3] = -VOXEL_MM * (np.array(SHAPE) - 1) / 2.0
MASK = nib.Nifti1Image(np.ones(SHAPE, dtype=np.int32), AFF)
MASK_BOOL = np.ones(SHAPE, dtype=bool); ZOOMS = np.full(3, VOXEL_MM)
SITE_IJK = [(10, 15, 15), (20, 15, 15), (15, 9, 18), (15, 21, 12), (15, 15, 22)]
SITE_WEIGHT = [1.0, 1.0, 0.75, 0.75, 0.5]
READ_AT = SITE_IJK[0]


def truth_field():
    out = np.zeros(SHAPE)
    grid = np.indices(SHAPE).astype(float)
    for (i, j, k), w in zip(SITE_IJK, SITE_WEIGHT):
        d2 = (grid[0] - i) ** 2 + (grid[1] - j) ** 2 + (grid[2] - k) ** 2
        out = np.maximum(out, PEAK_G * w * np.exp(-d2 / (2 * RADIUS_VOX**2)))
    return out


TRUTH = truth_field()
TRUE_G = float(TRUTH[READ_AT])
READ_MM = nib.affines.apply_affine(AFF, np.asarray(READ_AT))


def one(seed):
    rng = np.random.default_rng(seed)
    studies, reported_g, truth_at_focus = [], [], []
    for k in range(N_STUDIES):
        n = int(rng.integers(20, 41))
        noise = ndimage.gaussian_filter(rng.standard_normal(SHAPE), SMOOTH_VOX)
        noise *= 1.0 / (noise.std() + 1e-12)
        gmap = TRUTH + noise / np.sqrt(n)
        foci, _ = reporting.report_peaks((gmap * np.sqrt(n))[MASK_BOOL], MASK_BOOL, SHAPE,
                                         ZOOMS, "cluster", "max")
        meta = {"sample_sizes": [n]}
        points = []
        for ijk, zv in foci:
            mm = nib.affines.apply_affine(AFF, ijk)
            points.append({"space": "MNI", "coordinates": [float(v) for v in mm],
                           "values": [{"kind": "Z", "value": float(zv)}]})
            if np.linalg.norm(mm - READ_MM) <= NEAR_MM:
                reported_g.append(zv / np.sqrt(n))       # g implied by the reported statistic
                truth_at_focus.append(TRUTH[tuple(ijk)])  # truth where the focus actually is
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": points}]})
    if len(reported_g) < 2:
        return None

    pos = int(np.ravel_multi_index(READ_AT, SHAPE))
    ss = Studyset({"id": "d", "name": "d", "studies": studies}, target=None, mask=MASK)
    out = {}
    for name, selection in (("none", "none"), ("shipped", None)):
        kwargs = dict(fwhm=10.0, mask=MASK, null_method="none", use_images=False,
                      peak_bias=None)
        if selection == "none":
            kwargs["selection_model"] = "none"
        est = CBES(**kwargs)
        g = est.fit(ss).get_map("g", return_type="array").ravel()[pos]
        out[name] = float(g)
    return (float(np.mean(reported_g)), float(np.mean(truth_at_focus)),
            out["none"], out["shipped"], len(reported_g) / N_STUDIES)


if __name__ == "__main__":
    print(f"{N_STUDIES} coordinate-only studies, {N_SIMS} replications, "
          f"read-out voxel truth g = {TRUE_G:.3f}\n")
    rows = [r for r in Parallel(n_jobs=8)(delayed(one)(s) for s in range(N_SIMS))
            if r is not None]
    a = np.array(rows, dtype=float)
    reported, at_focus, none, shipped, share = (a[:, i] for i in range(5))
    print(f"{'stage':44s} {'value':>7s} {'cumulative':>11s} {'step':>7s}")
    print(f"{'1. truth at the read-out voxel':44s} {TRUE_G:7.3f} {0.0:+11.3f} {'':>7s}")
    print(f"{'   truth where the foci actually landed':44s} {at_focus.mean():7.3f} "
          f"{at_focus.mean()-TRUE_G:+11.3f} {at_focus.mean()-TRUE_G:+7.3f}   <- localisation")
    print(f"{'2. mean g implied by the reported peaks':44s} {reported.mean():7.3f} "
          f"{reported.mean()-TRUE_G:+11.3f} {reported.mean()-at_focus.mean():+7.3f}"
          f"   <- winner's curse")
    print(f"{'3. pooled, selection_model=none':44s} {none.mean():7.3f} "
          f"{none.mean()-TRUE_G:+11.3f} {none.mean()-reported.mean():+7.3f}"
          f"   <- weighting + kernel")
    print(f"{'4. pooled, shipped selection model':44s} {shipped.mean():7.3f} "
          f"{shipped.mean()-TRUE_G:+11.3f} {shipped.mean()-none.mean():+7.3f}"
          f"   <- selection model")
    print(f"\nfoci landing within {NEAR_MM:.0f} mm of the read-out voxel: "
          f"{share.mean():.2f} per study")
    print("\nA positive step at stage 4 would mean the censored likelihood pushes the estimate")
    print("up, when correcting for unobserved mass below a threshold should pull it down.")
