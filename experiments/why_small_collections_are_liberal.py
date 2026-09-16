"""Why is the family-wise rate 0.15 at twelve studies and 0.075 at twenty?

The boundary artefact is gone, the permutation distribution is not degenerate (200 distinct
maxima of 200, cv 0.08) and the generalized-Pareto tail makes no difference. The voxelwise rate
is nominal in every arm, so whatever is left acts only on the maximum.

Under a valid null the observed maximum is an exchangeable draw from the permuted maxima, so its
percentile within them is uniform on [0, 1] and exceeds 0.95 five times in a hundred. A null that
is too narrow shows up as that percentile piling up near 1. Reported here per arm, with the
pieces that would explain it:

  * **the percentile itself**, which says how badly and in which direction;
  * **the smoothness of the observed and permuted statistic maps**, because scrambling values
    among voxels destroys the image channel's spatial autocorrelation. That predicts the
    *opposite* sign -- a rougher permuted field has more effectively independent tests and a
    larger maximum -- so if the permuted maps are rougher and the test is still liberal, the
    smoothness is not the cause and is masking something bigger;
  * **the count of usable voxels in each**, since a maximum taken over fewer voxels is smaller,
    and a fit that returns fewer finite values under permutation than in the observed map would
    produce exactly this with no other mechanism needed.
"""
import os, sys, tempfile, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import logging; logging.getLogger("nimare").setLevel(logging.ERROR)
import numpy as np
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma.effectsize import CBES

N_SIMS = int(os.environ.get("NSIMS", 12))
N_ITERS = int(os.environ.get("NITERS", 200))
ERODE = int(os.environ.get("ERODE", 2))
ARMS = [tuple(int(v) for v in pair.split(":"))
        for pair in os.environ.get("ARMS", "12:2,20:2").split(",")]
ZOOMS, EXTENT = 4.0, 36.0
HALF = int(np.ceil(EXTENT / ZOOMS))
SHAPE = (2 * HALF + 1,) * 3
AFFINE = np.array([[ZOOMS, 0, 0, -ZOOMS * HALF], [0, ZOOMS, 0, -ZOOMS * HALF],
                   [0, 0, ZOOMS, -ZOOMS * HALF], [0, 0, 0, 1.0]])
THRESHOLDS = [2.3263, 3.0902, 3.2905, 4.2649]

_volume = np.zeros(SHAPE, np.int32)
_volume[ERODE:-ERODE, ERODE:-ERODE, ERODE:-ERODE] = 1
MASK = nib.Nifti1Image(_volume if ERODE else np.ones(SHAPE, np.int32), AFFINE)
MASK_BOOL = np.asarray(MASK.dataobj) > 0


def roughness(flat):
    """Mean squared first difference over the mask, as a scale-free measure of how rough."""
    volume = np.zeros(MASK_BOOL.shape)
    volume[MASK_BOOL] = flat
    total, count = 0.0, 0
    for axis in range(3):
        diff = np.diff(volume, axis=axis)
        keep = np.diff(MASK_BOOL.astype(int), axis=axis) == 0
        keep &= np.take(MASK_BOOL, range(volume.shape[axis] - 1), axis=axis)
        total += float((diff[keep] ** 2).sum())
        count += int(keep.sum())
    spread = float(np.var(flat))
    return (total / max(count, 1)) / max(spread, 1e-12)


if __name__ == "__main__":
    print(f"global null, {N_SIMS} simulations, {N_ITERS} permutations, "
          f"mask eroded by {ERODE} ({int(MASK_BOOL.sum())} voxels)\n", flush=True)
    print(f"{'arm':>18} {'pct of obs max':>15} {'>0.95':>7} {'obs max':>8} {'null 95th':>10} "
          f"{'rough obs':>10} {'rough null':>11} {'finite obs':>11} {'finite null':>12}")
    for n_studies, n_images in ARMS:
        pcts, obs_maxima, null_95 = [], [], []
        rough_obs, rough_null, finite_obs, finite_null = [], [], [], []
        for sim in range(N_SIMS):
            collection = create_effect_size_coordinate_studyset(
                [(0, 0, 0)], effect_sizes=0.0, n_studies=n_studies, sample_size=(20, 40),
                seed=41000 + sim, simulate_field=True, n_image_studies=n_images,
                image_dir=tempfile.mkdtemp(), noise_extent=EXTENT, field_zooms=ZOOMS,
                blob_fwhm=10.0,
                threshold_z=[THRESHOLDS[i % len(THRESHOLDS)] for i in range(n_studies)],
            )
            est = CBES(mask=MASK, null_method="permute-images", n_iters=N_ITERS, seed=sim,
                       threshold="reporting_threshold")
            result = est.fit(collection)
            observed = np.abs(result.get_map("z", return_type="array").ravel())
            maxima = np.asarray(est.null_distributions_.get(
                "values_level-voxel_corr-fwe_method-montecarlo", []), dtype=float)
            if not maxima.size:
                continue
            top = float(np.nanmax(observed))
            pcts.append(float(np.mean(maxima <= top)))
            obs_maxima.append(top)
            null_95.append(float(np.nanpercentile(maxima, 95)))
            finite_obs.append(float(np.isfinite(observed).mean()))
            rough_obs.append(roughness(np.nan_to_num(observed)))

            # One permuted map, refitted here so its roughness and finite share are visible.
            rng = np.random.default_rng(10_000 + sim)
            _, z_null = est._statistic(
                est._focus_table_, est._sample_sizes_, est._reported_thresholds_,
                est._permute_image_values(rng),
            )
            null_flat = np.abs(z_null)
            finite_null.append(float(np.isfinite(null_flat).mean()))
            rough_null.append(roughness(np.nan_to_num(null_flat)))

        if not pcts:
            continue
        label = f"{n_studies} studies, {n_images} img"
        print(f"{label:>18} {np.mean(pcts):15.3f} {np.mean(np.array(pcts) > 0.95):7.3f} "
              f"{np.mean(obs_maxima):8.3f} {np.mean(null_95):10.3f} "
              f"{np.mean(rough_obs):10.4f} {np.mean(rough_null):11.4f} "
              f"{np.mean(finite_obs):11.4f} {np.mean(finite_null):12.4f}", flush=True)

    print("\nA mean percentile above 0.5 means the observed maximum sits high in its own null. "
          "If the\npermuted maps are rougher, smoothness is not the cause: that predicts the "
          "opposite sign.")
