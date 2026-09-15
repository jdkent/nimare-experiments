"""How much does the fit depend on the reporting threshold it was given or inferred?

Section 22 established that `prevalence` is the fraction of covering studies that reported,
inflated by however much of the silence the censoring term can explain -- and that the inflation
is governed entirely by where the fitted magnitude sits relative to the *assumed* cutoff. So the
cutoff is not a nuisance parameter here; it is the thing that sets the answer.

The protocol says to hand the real threshold in through metadata rather than letting it be
inferred from the smallest reported value. Every bed this session let it be inferred, because
`threshold="study-min"` is the default. This measures what that cost.

Under cluster-extent reporting the reported focus is a cluster *maximum*, well above the
cluster-forming cut, so the smallest reported value is the smallest cluster maximum and
`study-min` should infer a cutoff above the truth -- which would inflate the prevalence.

Scored against a known prevalence that ranges over the sweep, because in a bed where the true
prevalence is 1 an inflated estimate is accidentally closer to the truth and the bed cannot tell
a right answer from a lucky one.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from joblib import Parallel, delayed
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.utils import mm2vox
import reporting

SHAPE, VOXEL_MM, SMOOTH_VOX, RADIUS_VOX = (30, 30, 30), 4.0, 0.8, 2.5
N_SIMS, N_STUDIES = int(os.environ.get("NSIMS", 30)), 20
AFF = np.eye(4); AFF[:3, :3] *= VOXEL_MM
AFF[:3, 3] = -VOXEL_MM * (np.array(SHAPE) - 1) / 2.0
MASK = nib.Nifti1Image(np.ones(SHAPE, dtype=np.int32), AFF)
MASK_BOOL = np.ones(SHAPE, dtype=bool); ZOOMS = np.full(3, VOXEL_MM)

SITE_IJK = [(9, 15, 15), (21, 15, 15), (15, 9, 18), (15, 21, 12)]
SITE_PREVALENCE = [0.25, 0.50, 0.75, 1.00]
TRUE_MU = 0.70


def blob(centre):
    grid = np.indices(SHAPE).astype(float)
    d2 = sum((grid[i] - centre[i]) ** 2 for i in range(3))
    return TRUE_MU * np.exp(-d2 / (2 * RADIUS_VOX**2))


BLOBS = [blob(c) for c in SITE_IJK]


def one(seed, threshold):
    rng = np.random.default_rng(seed)
    studies = []
    for k in range(N_STUDIES):
        n = int(rng.integers(20, 41))
        has = rng.random(len(SITE_IJK)) < np.array(SITE_PREVALENCE)
        field = np.zeros(SHAPE)
        for b, h in zip(BLOBS, has):
            if h:
                field = np.maximum(field, b)
        t_map = reporting.study_t_field(field, n, SMOOTH_VOX, rng, shape=SHAPE)
        foci, _ = reporting.report_peaks(t_map[MASK_BOOL], MASK_BOOL, SHAPE, ZOOMS,
                                         "cluster", "max")
        meta = {"sample_sizes": [n]}
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": [
                {"space": "MNI",
                 "coordinates": [float(v) for v in nib.affines.apply_affine(AFF, ijk)],
                 "values": [{"kind": "T", "value": float(tv)}]} for ijk, tv in foci]}]})
    est = CBES(fwhm=10.0, mask=MASK, null_method="none", use_images=False, peak_bias=None,
               threshold=threshold)
    res = est.fit(Studyset({"id": "th", "name": "th", "studies": studies},
                           target=None, mask=MASK))
    pi = res.get_map("prevalence", return_type="array").ravel()
    g = res.get_map("g", return_type="array").ravel()
    idx = [int(np.ravel_multi_index(c, SHAPE)) for c in SITE_IJK]
    cut = getattr(est, "_cutoffs_z_", None)
    return (np.array([pi[i] for i in idx]), np.array([g[i] for i in idx]),
            float(np.nanmedian(cut.values)) if cut is not None else np.nan)


#: The cluster-forming cut the bed actually applied, on the z scale the estimator works on.
FORMING_Z = reporting.CLUSTER_FORMING_Z

if __name__ == "__main__":
    reporting.assert_statistic_convention(
        reporting.study_t_field(np.zeros(SHAPE), 30, SMOOTH_VOX,
                                np.random.default_rng(5), shape=SHAPE), 30, "T")
    print("statistic convention check passed: studies report a t on n - 1 degrees of freedom")
    print(f"{N_STUDIES} coordinate-only studies, {N_SIMS} replications, true mu {TRUE_MU}")
    print(f"the bed's cluster-forming cut is z = {FORMING_Z:.4f}\n")
    print(f"{'threshold setting':26s} {'median cut':>11s}   "
          + " ".join(f"pi@{p:.2f}" for p in SITE_PREVALENCE) + f"   {'mean g':>7s}")
    print(f"{'(the truth)':26s} {'':>11s}   "
          + " ".join(f"{p:7.2f}" for p in SITE_PREVALENCE) + f"   {TRUE_MU:7.2f}")
    for label, threshold in (("study-min (the default)", "study-min"),
                             ("pooled-min", "pooled-min"),
                             ("the true forming cut", FORMING_Z),
                             ("the library default 3.2905", 3.2905267314919255)):
        rows = [r for r in Parallel(n_jobs=8)(delayed(one)(s, threshold)
                                              for s in range(N_SIMS)) if r is not None]
        pi = np.array([r[0] for r in rows]); g = np.array([r[1] for r in rows])
        cut = np.array([r[2] for r in rows])
        print(f"{label:26s} {np.nanmean(cut):11.3f}   "
              + " ".join(f"{v:7.3f}" for v in pi.mean(0))
              + f"   {g.mean():7.3f}", flush=True)
    print("\nPredicted: study-min infers a cut above the forming one, because a reported focus")
    print("is a cluster maximum, and a high cut inflates the prevalence.")
