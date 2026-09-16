"""False-positive rate of the redesigned CBES under a global null. Never measured since.

Every error rate on record was taken under the old *arrangement* null, which shuffled reported
effect sizes among the foci of an analysis. That null no longer exists: the redesign reads no
magnitude from a table, so the shipped null (`permute-images`) scrambles each image study's values
among that study's own voxels and holds the silence pattern fixed. The mechanism is different, the
old figures do not transfer, and the familywise defect they exposed was a property of a null that
has been deleted.

So this re-establishes the rate from scratch. Nothing here has an effect anywhere: the field is
noise, the images are noise, and each study reports whatever noise peaks clear its own threshold.
Thresholds vary from study to study, as in a literature search.

Two quantities, easy to confuse, so both are reported:

  * **uncorrected** -- the mean *share of voxels* with p < 0.05. Should sit near 0.05.
  * **familywise** -- the share of *simulations* with any voxel surviving FWE correction. Should
    sit at or below 0.05.

No cap on peaks anywhere: each study reports whatever clears its own cut, so the count is a
consequence of the noise rather than a fixed parameter.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import logging; logging.getLogger("nimare").setLevel(logging.ERROR)
import tempfile
import numpy as np
import nibabel as nib
from nimare.correct import FWECorrector
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma.effectsize import CBES

N_SIMS = int(os.environ.get("NSIMS", 40))
N_ITERS = int(os.environ.get("NITERS", 200))
ZOOMS, EXTENT = 4.0, 36.0
HALF = int(np.ceil(EXTENT / ZOOMS))
SHAPE = (2 * HALF + 1,) * 3
AFFINE = np.array([[ZOOMS, 0, 0, -ZOOMS * HALF], [0, ZOOMS, 0, -ZOOMS * HALF],
                   [0, 0, ZOOMS, -ZOOMS * HALF], [0, 0, 0, 1.0]])
#: Voxels to drop from the boundary of the simulated volume before fitting. The field is
#: smoothed inside a finite cube, so the kernel runs off the edge there and leaves both the
#: magnitude and the peak density inflated in the outer shell -- `is_the_image_permutation_
#: exchangeable.py` measures rho(|g|, silent count) at -0.125 over the whole cube and -0.0005
#: beyond two voxels. That coupling makes the permutation null non-exchangeable for reasons
#: belonging to the bed rather than the estimator, so the rate has to be measured both ways.
ERODE = int(os.environ.get("ERODE", 0))
#: Whether the family-wise correction extrapolates its tail with a generalized Pareto fit.
TAIL = os.environ.get("TAIL", "1") == "1"
#: "permute-images" scatters each image's values; "spatial-images" rearranges them with the
#: image's own autocorrelation preserved, which is the thing the family-wise rate is sensitive to.
NULL_METHOD = os.environ.get("NULL", "permute-images")
_volume = np.ones(SHAPE, np.int32)
if ERODE > 0:
    _volume[:] = 0
    _volume[ERODE:-ERODE, ERODE:-ERODE, ERODE:-ERODE] = 1
MASK = nib.Nifti1Image(_volume, AFFINE)
# A literature's worth of different reporting conventions, cycled over the studies.
THRESHOLDS = [2.3263, 3.0902, 3.2905, 4.2649]

#: Arms are overridable, because the interesting comparisons turned out to be about the image
#: share of the roster rather than the study count on its own.
if os.environ.get("ARMS"):
    ARMS = [(f"{n} studies, {k} image(s)", int(n), int(k))
            for n, k in (pair.split(":") for pair in os.environ["ARMS"].split(","))]
else:
    ARMS = [("20 studies, 2 images", 20, 2), ("20 studies, 1 image", 20, 1),
            ("12 studies, 2 images", 12, 2)]

if __name__ == "__main__":
    print(f"global null, {N_SIMS} simulations, {N_ITERS} permutations, nominal 0.05, "
          f"null {NULL_METHOD}")
    print(f"mask: {int(_volume.sum())} of {_volume.size} voxels "
          f"({'eroded by %d voxels' % ERODE if ERODE else 'the whole cube'})")
    print(f"binomial standard error on the familywise rate: "
          f"{np.sqrt(0.05 * 0.95 / N_SIMS):.3f}\n")
    print(f"{'arm':>22} {'uncorrected':>12} {'familywise':>11} {'refused':>8} {'n':>4}")
    for label, n_studies, n_images in ARMS:
        shares, any_survived, refused, used = [], [], 0, 0
        distinct, cvs = [], []
        for sim in range(N_SIMS):
            ss = create_effect_size_coordinate_studyset(
                [(0, 0, 0)], effect_sizes=0.0, n_studies=n_studies, sample_size=(20, 40),
                seed=41000 + sim, simulate_field=True, n_image_studies=n_images,
                image_dir=tempfile.mkdtemp(), noise_extent=EXTENT, field_zooms=ZOOMS,
                blob_fwhm=10.0, threshold_z=[THRESHOLDS[i % len(THRESHOLDS)]
                                             for i in range(n_studies)])
            est = CBES(mask=MASK, n_iters=N_ITERS, seed=sim, null_method=NULL_METHOD,
                       threshold="reporting_threshold")
            try:
                res = est.fit(ss)
            except ValueError:
                continue
            p = res.get_map("p", return_type="array").ravel()
            if not np.isfinite(p).any() or np.allclose(p, 1.0):
                refused += 1      # the null was refused as degenerate; it cannot reject
                continue
            used += 1
            shares.append(float(np.mean(p < 0.05)))
            # The generalized-Pareto tail is the remaining suspect for the residual inflation
            # at small study counts, now that the permutation distribution is known to be well
            # spread there (200 distinct maxima of 200, cv 0.08). TAIL=0 refits without it.
            corrected = FWECorrector(
                method="montecarlo", n_iters=N_ITERS, tail_approximation=TAIL
            ).transform(res)
            name = "logp_level-voxel_corr-FWE_method-montecarlo"
            logp = corrected.get_map(name, return_type="array").ravel()
            any_survived.append(bool(np.nanmax(logp) > -np.log10(0.05)))
            # The guard's diagnostics are written by `correct_fwe_montecarlo`, which the
            # corrector runs on the result's own estimator rather than on `est`, so reading
            # them off `est` gives nan. Take them from whichever object actually carries them.
            for holder in (getattr(corrected, "estimator", None),
                           getattr(res, "estimator", None), est):
                nulls = getattr(holder, "null_distributions_", {}) or {}
                if "max_statistic_distinct_values" in nulls:
                    distinct.append(nulls["max_statistic_distinct_values"])
                    cvs.append(nulls.get("max_statistic_cv", np.nan))
                    break
            else:
                distinct.append(np.nan)
                cvs.append(np.nan)
        if not used:
            print(f"{label:>22}  every fit refused the null ({refused}/{N_SIMS})")
            continue
        print(f"  {label}: {np.nanmean(distinct):.0f} distinct maxima of {N_ITERS}, "
              f"cv {np.nanmean(cvs):.3f}", flush=True)
        print(f"{label:>22} {np.mean(shares):12.4f} {np.mean(any_survived):11.3f} "
              f"{refused:8d} {used:4d}")
