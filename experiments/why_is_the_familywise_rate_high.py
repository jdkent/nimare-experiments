"""Why is the familywise rate 0.375 when the uncorrected rate is 0.053?

`fpr_after_redesign.py` finds the voxelwise rate calibrated and the familywise rate seven times
nominal. A defect that spares the marginal and wrecks the maximum is a defect in the dependence
structure, and there is only one place it can come from.

The fit is voxelwise: `z` at a voxel depends only on that voxel's data. So permuting each image's
values among its own voxels does not change the *multiset* of image values, only which voxel each
one is paired with. If the indicator pattern were spatially constant, the permuted `z` map would
be a permutation of the observed one and the two maxima would be identical by construction. The
maxima can differ at all only because the indicator is held fixed while the image moves -- so the
observed map is scored on the *actual* pairing of image value to indicator, and every permutation
on a random pairing.

That makes the null's validity rest entirely on one assumption: **across voxels, a study's image
magnitude is independent of how many studies were silent there.** If large image values sit
preferentially where the indicator is unusual, the actual pairing reaches higher than a random one
and the observed maximum beats the permuted maxima more often than it should -- exactly the
failure observed, and in the right direction, which the alternative explanation does not manage.
Scrambling destroys spatial smoothness, which would give the permuted field *more* effective
independent tests and a *larger* maximum, making the test conservative rather than liberal.

Three things measured here, under the same global null:

  * the association itself, `corr(|g| at a voxel, number of silent studies there)`, per image
    study. The null needs this to be zero.
  * where the observed maximum sits: is its silent count unusual compared with the map?
  * the two maxima distributions, observed against permuted, so the gap is quantified rather
    than inferred.

If the association is zero and the maxima still disagree, this explanation is wrong and the
smoothness of the fixed indicator pattern is the next suspect.
"""
import os, sys, tempfile, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import logging; logging.getLogger("nimare").setLevel(logging.ERROR)
import numpy as np
import nibabel as nib
from scipy import stats
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma.effectsize import CBES

N_SIMS = int(os.environ.get("NSIMS", 6))
N_ITERS = int(os.environ.get("NITERS", 200))
N_STUDIES = int(os.environ.get("NSTUDIES", 20))
N_IMAGES = int(os.environ.get("NIMAGES", 2))
ZOOMS, EXTENT = 4.0, 36.0
HALF = int(np.ceil(EXTENT / ZOOMS))
SHAPE = (2 * HALF + 1,) * 3
AFFINE = np.array([[ZOOMS, 0, 0, -ZOOMS * HALF], [0, ZOOMS, 0, -ZOOMS * HALF],
                   [0, 0, ZOOMS, -ZOOMS * HALF], [0, 0, 0, 1.0]])
MASK = nib.Nifti1Image(np.ones(SHAPE, np.int32), AFFINE)
THRESHOLDS = [2.3263, 3.0902, 3.2905, 4.2649]


def silent_count(est, n_voxels, n_coordinate_studies):
    """Per voxel, how many coordinate studies were silent there."""
    _, (col, _, sign) = est._coverage_
    silent = col[sign > 0]
    return np.bincount(silent, minlength=n_voxels).astype(float)


if __name__ == "__main__":
    print(f"global null, {N_SIMS} realisations, {N_STUDIES} studies, {N_IMAGES} images, "
          f"{N_ITERS} permutations\n", flush=True)
    print(f"{'sim':>4} {'corr(|g|, silent)':>18} {'p':>8} | {'obs max':>8} "
          f"{'null 95th':>10} {'null max':>9} | {'silent at argmax':>17} {'map median':>11}")
    beats, rows = 0, []
    for sim in range(N_SIMS):
        # The same call the error-rate bed makes, so the diagnosis applies to those numbers.
        collection = create_effect_size_coordinate_studyset(
            [(0, 0, 0)], effect_sizes=0.0, n_studies=N_STUDIES, sample_size=(20, 40),
            seed=41000 + sim, simulate_field=True, n_image_studies=N_IMAGES,
            image_dir=tempfile.mkdtemp(), noise_extent=EXTENT, field_zooms=ZOOMS,
            blob_fwhm=10.0,
            threshold_z=[THRESHOLDS[i % len(THRESHOLDS)] for i in range(N_STUDIES)],
        )
        est = CBES(mask=MASK, null_method="permute-images", n_iters=N_ITERS, seed=sim,
                   threshold="reporting_threshold")
        res = est.fit(collection)
        z = np.abs(res.get_map("z", return_type="array").ravel())
        n_voxels = z.size
        counts = silent_count(est, n_voxels, N_STUDIES - N_IMAGES)

        # The association the null's validity rests on, measured on the image that carries the
        # most weight. Spearman rather than Pearson: the question is monotone association, and
        # g has heavy tails at small samples.
        images = getattr(est, "_image_studies_", None) or {}
        assoc, assoc_p = np.nan, np.nan
        if images:
            values, _, usable = list(images.values())[0]
            arr = np.abs(np.asarray(values, dtype=float).ravel())
            use = np.asarray(usable).ravel() & np.isfinite(arr)
            if use.sum() > 100 and np.ptp(counts[use]) > 0:
                assoc, assoc_p = stats.spearmanr(arr[use], counts[use])

        maxima = np.asarray(est.null_distributions_.get(
            "values_level-voxel_corr-fwe_method-montecarlo", []), dtype=float)
        obs_max = float(np.nanmax(z))
        null_95 = float(np.nanpercentile(maxima, 95)) if maxima.size else np.nan
        null_max = float(np.nanmax(maxima)) if maxima.size else np.nan
        beats += int(np.isfinite(null_95) and obs_max > null_95)
        at_argmax = float(counts[int(np.nanargmax(z))])
        print(f"{sim:4d} {assoc:+18.4f} {assoc_p:8.4f} | {obs_max:8.3f} {null_95:10.3f} "
              f"{null_max:9.3f} | {at_argmax:17.1f} {np.median(counts):11.1f}", flush=True)
        rows.append((assoc, obs_max, null_95, at_argmax, float(np.median(counts))))

    arr = np.array(rows, dtype=float)
    print(f"\nobserved max beat the permuted 95th in {beats}/{N_SIMS} realisations")
    print(f"mean association: {np.nanmean(arr[:, 0]):+.4f}")
    print(f"mean observed max {np.nanmean(arr[:, 1]):.3f} against permuted 95th "
          f"{np.nanmean(arr[:, 2]):.3f}")
    print(f"silent count at the observed maximum {np.nanmean(arr[:, 3]):.1f} against a map "
          f"median of {np.nanmean(arr[:, 4]):.1f}")
    print("\nA non-zero association means the permutation re-pairs image values with indicators "
          "that\nwere not independent of them, and the null is not exchangeable.")
