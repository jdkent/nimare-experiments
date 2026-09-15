"""Which part of the pooling step adds +0.324 to the coordinates-only bias?

The stage decomposition put as much of the bias after the reported values are in hand as in the
winner's curse itself: the mean g implied by the reported peaks near the read-out voxel is 0.993,
and the pooled estimate there is 1.318. Half the bias is in the weighting and the kernel, which
is code rather than reporting practice.

Before believing that, the first job is to rule out the dull explanation: that the two stages are
averaging different things. Stage 2 took the mean over foci within 10 mm, chosen by hand; the
estimator's kernel support is whatever ``fwhm`` implies and reaches whatever it reaches. If the
estimator is simply averaging a different set of foci, there is no bias to chase.

So everything here is computed from the estimator's *own* contributions at the read-out voxel --
the same foci, the same weights, the same variances -- and the alternatives differ only in which
part of the weighting is applied:

  plain mean            every contributing focus, equal weight
  kernel only           the spatial weight, no inverse-variance term
  precision only        the inverse-variance term, no spatial weight
  both (the estimator)  what CBES reports

If 'plain mean' already sits near 1.3 then the gap is the set of foci, not the weighting, and the
lead was illusory. If it sits near 1.0 and one of the weighted versions climbs, that names the
mechanism.
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


def one(seed, fwhm):
    rng = np.random.default_rng(seed)
    studies = []
    for k in range(N_STUDIES):
        n = int(rng.integers(20, 41))
        noise = ndimage.gaussian_filter(rng.standard_normal(SHAPE), SMOOTH_VOX)
        noise *= 1.0 / (noise.std() + 1e-12)
        z = (TRUTH + noise / np.sqrt(n)) * np.sqrt(n)
        foci, _ = reporting.report_peaks(z[MASK_BOOL], MASK_BOOL, SHAPE, ZOOMS, "cluster", "max")
        meta = {"sample_sizes": [n]}
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": [
                {"space": "MNI",
                 "coordinates": [float(v) for v in nib.affines.apply_affine(AFF, ijk)],
                 "values": [{"kind": "Z", "value": float(zv)}]} for ijk, zv in foci]}]})
    est = CBES(fwhm=fwhm, mask=MASK, null_method="none", use_images=False, peak_bias=None,
               selection_model="none")
    result = est.fit(Studyset({"id": "p", "name": "p", "studies": studies},
                              target=None, mask=MASK))
    pos = int(np.ravel_multi_index(READ_AT, SHAPE))
    shipped = float(result.get_map("g", return_type="array").ravel()[pos])

    # The estimator's own contributions at this voxel: same foci, same weights, same variances.
    contributions, _, _ = est._accumulate(est._focus_table_, None)
    g_here, w_here, var_here = [], [], []
    for _, cols, weights, g, var_g in contributions:
        hit = cols == pos
        if hit.any():
            g_here.append(g[hit]); w_here.append(weights[hit]); var_here.append(var_g[hit])
    if not g_here:
        return None
    g_here = np.concatenate(g_here)
    w_here = np.concatenate(w_here)
    var_here = np.concatenate(var_here)

    plain = float(np.mean(g_here))
    kernel_only = float(np.sum(w_here * g_here) / np.sum(w_here))
    precision_only = float(np.sum(g_here / var_here) / np.sum(1.0 / var_here))
    both = float(np.sum(w_here / var_here * g_here) / np.sum(w_here / var_here))
    return plain, kernel_only, precision_only, both, shipped, g_here.size


if __name__ == "__main__":
    print(f"{N_STUDIES} coordinate-only studies, {N_SIMS} replications, "
          f"truth at the read-out voxel {TRUE_G:.3f}\n")
    print(f"{'fwhm':>5s} {'foci':>6s} {'plain mean':>11s} {'kernel only':>12s} "
          f"{'precision only':>15s} {'both':>7s} {'shipped':>8s}")
    for fwhm in (4.0, 10.0, 16.0):
        rows = [r for r in Parallel(n_jobs=8)(delayed(one)(s, fwhm) for s in range(N_SIMS))
                if r is not None]
        a = np.array(rows, dtype=float)
        print(f"{fwhm:5.1f} {a[:, 5].mean():6.1f} {a[:, 0].mean():11.3f} {a[:, 1].mean():12.3f} "
              f"{a[:, 2].mean():15.3f} {a[:, 3].mean():7.3f} {a[:, 4].mean():8.3f}", flush=True)
    print("\n'plain mean' near the truth with a weighted column above it names the weighting.")
    print("'plain mean' already high means the gap is which foci reach the voxel, not weighting.")
