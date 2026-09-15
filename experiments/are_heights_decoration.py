"""How much does the estimator lose if the reported heights carry no within-study information?

Section 17 claims the signal in a coordinate table is in whether and where foci appear, not in
how large the reported statistics are; section 18 claims the prevalence/magnitude split is
identified by the spread of study *power* (sample size, threshold) rather than by the heights.
Together they imply the statistic column is close to decoration: the estimator needs it to put
the fit on a g scale at all, but the variation *within* a study should buy almost nothing.

That is directly testable by destroying exactly that variation and nothing else. Three inputs
over the same collections:

  as reported      the table a paper would print
  study-flattened  every focus in a study relabelled with that study's own mean height, so
                   between-study variation survives and within-study variation is gone
  all-flattened    every focus in the collection relabelled with the grand mean height, so
                   no height information of any kind survives

Locations, counts, study membership, sample sizes and thresholds are untouched in all three.
If 'study-flattened' matches 'as reported', the within-study heights are decoration. If
'all-flattened' also matches, the heights are decoration entirely and the estimator is really a
count-and-location method wearing a magnitude model.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import ndimage, stats
from joblib import Parallel, delayed
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
import reporting

SHAPE, VOXEL_MM, SMOOTH_VOX, RADIUS_VOX = (30, 30, 30), 4.0, 0.8, 2.5
PEAK_G, N_SIMS, N_STUDIES = 0.8, 30, 16
AFF = np.eye(4); AFF[:3, :3] *= VOXEL_MM
AFF[:3, 3] = -VOXEL_MM * (np.array(SHAPE) - 1) / 2.0
MASK = nib.Nifti1Image(np.ones(SHAPE, dtype=np.int32), AFF)
MASK_BOOL = np.ones(SHAPE, dtype=bool); ZOOMS = np.full(3, VOXEL_MM)
SITE_IJK = [(10, 15, 15), (20, 15, 15), (15, 9, 18), (15, 21, 12), (15, 15, 22)]
SITE_WEIGHT = [1.0, 0.875, 0.75, 0.625, 0.5]


def truth_field():
    out = np.zeros(SHAPE)
    grid = np.indices(SHAPE).astype(float)
    for (i, j, k), w in zip(SITE_IJK, SITE_WEIGHT):
        d2 = (grid[0] - i) ** 2 + (grid[1] - j) ** 2 + (grid[2] - k) ** 2
        out = np.maximum(out, PEAK_G * w * np.exp(-d2 / (2 * RADIUS_VOX**2)))
    return out


TRUTH = truth_field()
TRUTH_VEC = TRUTH.ravel()


def build(records, mode):
    """Assemble a studyset, optionally flattening the reported heights."""
    if mode == "all-flattened":
        grand = np.mean([z for _, pts in records for _, z in pts])
    studies = []
    for k, (n, pts) in enumerate(records):
        if mode == "study-flattened":
            value = np.mean([abs(z) for _, z in pts])
            used = [(ijk, np.sign(z) * value) for ijk, z in pts]
        elif mode == "all-flattened":
            used = [(ijk, np.sign(z) * abs(grand)) for ijk, z in pts]
        else:
            used = pts
        meta = {"sample_sizes": [n]}
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": [
                {"space": "MNI",
                 "coordinates": [float(v) for v in nib.affines.apply_affine(AFF, ijk)],
                 "values": [{"kind": "Z", "value": float(z)}]} for ijk, z in used]}]})
    return Studyset({"id": "h", "name": "h", "studies": studies}, target=None, mask=MASK)


def one(seed):
    rng = np.random.default_rng(seed)
    records = []
    for _ in range(N_STUDIES):
        n = int(rng.integers(20, 41))
        noise = ndimage.gaussian_filter(rng.standard_normal(SHAPE), SMOOTH_VOX)
        noise *= 1.0 / (noise.std() + 1e-12)
        z = (TRUTH + noise / np.sqrt(n)) * np.sqrt(n)
        foci, _ = reporting.report_peaks(z[MASK_BOOL], MASK_BOOL, SHAPE, ZOOMS, "cluster", "max")
        if foci:
            records.append((n, foci))
    if len(records) < 4:
        return None
    out = []
    for mode in ("as reported", "study-flattened", "all-flattened"):
        est = CBES(fwhm=10.0, mask=MASK, null_method="none", use_images=False, peak_bias=None)
        res = est.fit(build(records, mode))
        g = res.get_map("g", return_type="array").ravel()
        pi = res.get_map("prevalence", return_type="array").ravel()
        cov = np.isfinite(g) & (g != 0)
        out.append((stats.pearsonr(g[cov], TRUTH_VEC[cov])[0],
                    float(np.mean([g[int(np.ravel_multi_index(s, SHAPE))] for s in SITE_IJK])),
                    float(np.mean([pi[int(np.ravel_multi_index(s, SHAPE))] for s in SITE_IJK])),
                    stats.pearsonr(pi[cov], TRUTH_VEC[cov])[0]))
    return out


if __name__ == "__main__":
    print(f"{N_STUDIES} coordinate-only studies, {N_SIMS} replications")
    print(f"mean true g at the five sites: "
          f"{np.mean([TRUTH[s] for s in SITE_IJK]):.3f}\n")
    rows = [r for r in Parallel(n_jobs=8)(delayed(one)(s) for s in range(N_SIMS))
            if r is not None]
    a = np.array(rows, dtype=float)            # (sims, modes, metrics)
    names = ("as reported", "study-flattened", "all-flattened")
    print(f"{'height input':18s} {'r(g, truth)':>12s} {'mean g at sites':>16s} "
          f"{'mean pi at sites':>17s} {'r(pi, truth)':>13s}")
    for m, name in enumerate(names):
        print(f"{name:18s} {a[:, m, 0].mean():12.3f} {a[:, m, 1].mean():16.3f} "
              f"{a[:, m, 2].mean():17.3f} {a[:, m, 3].mean():13.3f}", flush=True)
    print("\nDifferences from 'as reported', paired across replications:")
    for m, name in enumerate(names[1:], start=1):
        d = a[:, m, 0] - a[:, 0, 0]
        t = stats.ttest_rel(a[:, m, 0], a[:, 0, 0])
        print(f"  {name:18s} r changes {d.mean():+.4f} "
              f"(sd {d.std(ddof=1):.4f}, paired p {t.pvalue:.3f})")
    print("\nIf flattening the heights within a study costs almost nothing, the statistic")
    print("column is decoration and the estimator is a count-and-location method.")
