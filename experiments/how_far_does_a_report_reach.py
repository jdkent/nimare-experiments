"""How far should a reported focus assert |g| >= c? Not 0 mm, and not the silence radius.

The autopsy of a dead hypothesis. `alpha`-flattening the report limb -- the probe for "a report
is over-credited because a reported peak is a local maximum, not any exceedance" -- moved the
bias the wrong way and monotonically:

    alpha      1.0     1.5     2.0     3.0     5.0     0.7     0.5     0.3     0.1
    bias     -0.039  -0.024  -0.008  +0.011  +0.024  -0.045  -0.049  -0.055  -0.061
    rmse      0.0697  0.0762  0.0977  0.1338  0.1696  0.0733  0.0771  0.0826  0.0909

So the report limb is not over-credited, alpha = 1 is the rmse optimum, and buying the bias with
alpha costs 40% of rmse. The knob is the wrong lever.

What the residual bias looks like instead. A study's peak is *displaced* from the effect -- a
local maximum sits where the noise happened to help. So at the true focus, a study that plainly
detected the region usually names a neighbouring voxel, which the model then reads as **no
information at all** (sign 0) rather than as a detection. The studies still carrying an indicator
there are the ones whose whole 20 mm neighbourhood came up empty -- a sample enriched for
failures. That is a mechanism for a downward bias exactly where the effect is, which is where it
was measured.

Asserting the report across the full silence radius was already measured and is catastrophic
(rmse 0.457): 20 mm is the *reporting extent*, and a peak 18 mm away says nothing about here.
But there is an untested middle. A peak's position carries a *localisation* error of a few mm,
far smaller than the extent it stands in for, so asserting the indicator over that error is
nearly as safe as asserting it at the named voxel while recovering the studies the model is
currently discarding.

Sweep it. Two radii, not one: `coverage_radius` for the silence and `_report_radius` for the
report, the first fixed at 20 mm throughout.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
import tempfile
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn.maskers import NiftiMasker
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma.effectsize import CBES

SHAPE, ZOOMS, EXTENT, BLOB = (21, 21, 21), 4.0, 40.0, 10.0
AFFINE = np.array([[ZOOMS, 0, 0, -EXTENT], [0, ZOOMS, 0, -EXTENT],
                   [0, 0, ZOOMS, -EXTENT], [0, 0, 0, 1.0]])
MASK = nib.Nifti1Image(np.ones(SHAPE, np.int32), AFFINE)
MASKER = NiftiMasker(MASK).fit()
TRUE_G = 0.5
N_SEEDS = int(os.environ.get("NSEEDS", 8))
RADII = [float(r) for r in os.environ.get("RADII", "0,4,6,8,12,20").split(",")]

grid = np.stack(np.indices(SHAPE), -1) * ZOOMS - EXTENT
sd = BLOB / (2 * np.sqrt(2 * np.log(2)))
TRUTH = MASKER.transform(nib.Nifti1Image(
    (TRUE_G * np.exp(-(grid ** 2).sum(-1) / (2 * sd ** 2))).astype(np.float32), AFFINE)).ravel()
STRATA = [("quiet", TRUTH < 0.05), ("middle", (TRUTH >= 0.05) & (TRUTH < 0.25)),
          ("effect", TRUTH >= 0.25)]

if __name__ == "__main__":
    built = []
    for seed in range(N_SEEDS):
        built.append(create_effect_size_coordinate_studyset(
            [(0, 0, 0)], effect_sizes=TRUE_G, n_studies=20, sample_size=(20, 40), tau=0.1,
            seed=seed, simulate_field=True, n_image_studies=2,
            image_dir=tempfile.mkdtemp(), noise_extent=EXTENT, field_zooms=ZOOMS,
            blob_fwhm=BLOB))

    print(f"true g = {TRUE_G} at the focus, {N_SEEDS} collections, 20 studies, 2 images, "
          f"silence radius 20 mm\n")
    print(f"{'report radius':>14} " + " ".join(f"{'rmse ' + n:>12}" for n, _ in STRATA)
          + " " + " ".join(f"{'bias ' + n:>12}" for n, _ in STRATA))
    results = {}
    for radius in RADII:
        maps = []
        for ss in built:
            est = CBES(mask=MASK, null_method="none", threshold="reporting_threshold")
            est._report_radius = radius
            res = est.fit(ss)
            maps.append(np.abs(res.get_map("g", return_type="array").ravel()))
        results[radius] = maps
        rmse = [np.mean([float(np.sqrt(np.mean((m[sel] - TRUTH[sel]) ** 2))) for m in maps])
                for _, sel in STRATA]
        bias = [np.mean([float(np.mean(m[sel] - TRUTH[sel])) for m in maps])
                for _, sel in STRATA]
        label = "named voxel" if radius == 0 else f"{radius:g} mm"
        print(f"{label:>14} " + " ".join(f"{c:12.4f}" for c in rmse)
              + " " + " ".join(f"{b:+12.4f}" for b in bias))

    print("\nPaired against the named voxel, rmse per stratum:")
    base = results[RADII[0]]
    for radius in RADII[1:]:
        line = []
        for _, sel in STRATA:
            a = np.array([float(np.sqrt(np.mean((m[sel] - TRUTH[sel]) ** 2)))
                          for m in results[radius]])
            b = np.array([float(np.sqrt(np.mean((m[sel] - TRUTH[sel]) ** 2))) for m in base])
            pv = stats.ttest_rel(a, b).pvalue if np.any(a != b) else 1.0
            line.append(f"{np.mean(a - b):+8.4f} (p {pv:.4f})")
        print(f"  {radius:5g} mm  " + "  ".join(line))
