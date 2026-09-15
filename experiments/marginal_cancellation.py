"""Is the bias cancellation in `g_marginal` structural, or an accident of one regime?

The estimator's docstring says `g_marginal = g * prevalence` works because two biases cancel:
`g` runs high (peak-height selection) while `prevalence` runs low (a study that had the effect
but missed its threshold is pulled toward "did not have it"). The evidence was an HCP collection
whose true prevalence was 1 while `prevalence` read 0.68, and whose `g_marginal` matched an IBMA.

The coverage bed contradicts that. There, too, true prevalence is 1 -- but `prevalence` reads
0.994, so there is nothing for it to absorb and `g_marginal` inherits the whole +0.5 bias in `g`.
The difference between the two beds is *reporting density*: in the coverage bed the site is
strong and nearly every study reports a focus near it, while in the HCP collection most studies
were silent at any given voxel.

That suggests the cancellation is not structural at all. `prevalence` is pushed down by silence,
`g` is pushed up by selection, and only the first depends on how often studies stay silent. Where
reporting is dense, `prevalence` saturates at 1, the downward push vanishes, and `g_marginal` is
as wrong as `g`. Where reporting is sparse the two can offset. If that is right, `g_marginal` is
*least* trustworthy exactly at the strongest voxels -- the ones readers interpret.

The test varies both axes independently in one fit: four sites crossing two true prevalences
with two true magnitudes, so reporting density varies for reasons the truth distinguishes. The
target is the marginal effect, true prevalence times true magnitude, which is what an IBMA over
the same studies would estimate.
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
N_SIMS, N_STUDIES = 40, 24
AFF = np.eye(4); AFF[:3, :3] *= VOXEL_MM
AFF[:3, 3] = -VOXEL_MM * (np.array(SHAPE) - 1) / 2.0
MASK = nib.Nifti1Image(np.ones(SHAPE, dtype=np.int32), AFF)
MASK_BOOL = np.ones(SHAPE, dtype=bool); ZOOMS = np.full(3, VOXEL_MM)

#: (voxel, true magnitude among studies that have it, true prevalence)
SITES = [
    ((8, 15, 15),  0.80, 1.00),   # strong, universal   -> dense reporting
    ((22, 15, 15), 0.80, 0.40),   # strong, uncommon    -> sparse reporting, strong when present
    ((15, 8, 18),  0.45, 1.00),   # weak, universal     -> sparse reporting, weak when present
    ((15, 22, 12), 0.45, 0.40),   # weak, uncommon      -> very sparse
]


def blob(centre, amplitude):
    grid = np.indices(SHAPE).astype(float)
    d2 = sum((grid[i] - centre[i]) ** 2 for i in range(3))
    return amplitude * np.exp(-d2 / (2 * RADIUS_VOX**2))


BLOBS = [blob(c, mu) for c, mu, _ in SITES]
TRUE_MU = np.array([mu for _, mu, _ in SITES])
TRUE_PI = np.array([pi for _, _, pi in SITES])
TRUE_MARGINAL = TRUE_MU * TRUE_PI


def one(seed):
    rng = np.random.default_rng(seed)
    studies, reported_near = [], np.zeros(len(SITES))
    for k in range(N_STUDIES):
        n = int(rng.integers(20, 41))
        has = rng.random(len(SITES)) < TRUE_PI
        field = np.zeros(SHAPE)
        for b, h in zip(BLOBS, has):
            if h:
                field = np.maximum(field, b)
        noise = ndimage.gaussian_filter(rng.standard_normal(SHAPE), SMOOTH_VOX)
        noise *= 1.0 / (noise.std() + 1e-12)
        gmap = field + noise / np.sqrt(n)
        foci, _ = reporting.report_peaks((gmap * np.sqrt(n))[MASK_BOOL], MASK_BOOL, SHAPE,
                                         ZOOMS, "cluster", "max")
        meta = {"sample_sizes": [n]}
        points = []
        for ijk, zv in foci:
            mm = nib.affines.apply_affine(AFF, ijk)
            points.append({"space": "MNI", "coordinates": [float(v) for v in mm],
                           "values": [{"kind": "Z", "value": float(zv)}]})
            for s, (centre, _, _) in enumerate(SITES):
                if np.linalg.norm(np.asarray(ijk) - np.asarray(centre)) <= 3:
                    reported_near[s] += 1
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": points}]})

    est = CBES(fwhm=10.0, mask=MASK, null_method="none", use_images=False, peak_bias=None)
    res = est.fit(Studyset({"id": "m", "name": "m", "studies": studies},
                           target=None, mask=MASK))
    g = res.get_map("g", return_type="array").ravel()
    pi = res.get_map("prevalence", return_type="array").ravel()
    marg = (res.get_map("g_marginal", return_type="array").ravel()
            if "g_marginal" in res.maps else g * pi)
    idx = [int(np.ravel_multi_index(c, SHAPE)) for c, _, _ in SITES]
    return (np.array([g[i] for i in idx]), np.array([pi[i] for i in idx]),
            np.array([marg[i] for i in idx]), reported_near / N_STUDIES)


if __name__ == "__main__":
    print(f"{N_STUDIES} coordinate-only studies, {N_SIMS} replications\n")
    rows = [r for r in Parallel(n_jobs=6)(delayed(one)(s) for s in range(N_SIMS))
            if r is not None]
    g = np.array([r[0] for r in rows]); pi = np.array([r[1] for r in rows])
    marg = np.array([r[2] for r in rows]); dens = np.array([r[3] for r in rows])
    names = ["strong, universal", "strong, uncommon", "weak, universal", "weak, uncommon"]
    print(f"{'site':20s} {'rep dens':>8s} {'true mu':>8s} {'g':>7s} {'g bias':>7s} "
          f"{'true pi':>8s} {'pi':>6s} {'true marg':>9s} {'g_marg':>7s} {'marg bias':>9s}")
    for s, name in enumerate(names):
        print(f"{name:20s} {dens[:, s].mean():8.2f} {TRUE_MU[s]:8.2f} {g[:, s].mean():7.3f} "
              f"{g[:, s].mean()-TRUE_MU[s]:+7.3f} {TRUE_PI[s]:8.2f} {pi[:, s].mean():6.3f} "
              f"{TRUE_MARGINAL[s]:9.3f} {marg[:, s].mean():7.3f} "
              f"{marg[:, s].mean()-TRUE_MARGINAL[s]:+9.3f}", flush=True)
    print("\nIf the cancellation were structural, 'marg bias' would be small at every site.")
    print("If it depends on reporting density, it will be worst where 'rep dens' is highest --")
    print("that is, at the strongest and most consistent voxels, the ones readers interpret.")
