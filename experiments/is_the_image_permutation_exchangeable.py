"""Is `permute-images` exchangeable? The smallest test that can answer it.

`fpr_after_redesign.py` finds the voxelwise rate calibrated at 0.053 and the familywise rate at
0.375 against a nominal 0.05. A defect that spares the marginal and wrecks the maximum lives in
the dependence structure, and for this estimator there is exactly one candidate.

The fit is voxelwise -- `z` at a voxel depends only on that voxel's data -- so permuting an
image's values among its own voxels does not change the multiset of values, only which voxel each
is paired with. Were the indicator pattern spatially constant, the permuted `z` map would be a
permutation of the observed one and the maxima would agree by construction. They can differ only
because the indicator is held fixed while the image moves. So the observed map is scored on the
*actual* pairing of image value to indicator and every permutation on a random one, and the
null's validity reduces to a single question:

    across voxels, is a study's image magnitude independent of how many studies were silent
    there?

That question needs no CBES, no permutations and no correction -- only the simulator's own
studyset and a spatial query, which is why it belongs here rather than in a bed that takes an
hour. Image studies contribute no indicator, so the silence is counted over the coordinate
studies alone, and the magnitude is averaged over the image studies.

The competing explanation is worth stating because it fails on direction. Scrambling values among
voxels destroys spatial smoothness, which gives the permuted field *more* effectively independent
tests and therefore a *larger* maximum -- that makes a test conservative, not liberal, so it
cannot produce a familywise rate seven times nominal.
"""
import os, sys, tempfile, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import logging; logging.getLogger("nimare").setLevel(logging.ERROR)
import numpy as np
import nibabel as nib
from scipy import spatial, stats
from nimare.generate import create_effect_size_coordinate_studyset

N_SIMS = int(os.environ.get("NSIMS", 8))
N_STUDIES = int(os.environ.get("NSTUDIES", 20))
N_IMAGES = int(os.environ.get("NIMAGES", 2))
RADIUS = float(os.environ.get("RADIUS", 8.0))
ZOOMS, EXTENT = 4.0, 36.0
HALF = int(np.ceil(EXTENT / ZOOMS))
SHAPE = (2 * HALF + 1,) * 3
AFFINE = np.array([[ZOOMS, 0, 0, -ZOOMS * HALF], [0, ZOOMS, 0, -ZOOMS * HALF],
                   [0, 0, ZOOMS, -ZOOMS * HALF], [0, 0, 0, 1.0]])
MASK = nib.Nifti1Image(np.ones(SHAPE, np.int32), AFFINE)
THRESHOLDS = [2.3263, 3.0902, 3.2905, 4.2649]


def voxel_centres():
    grid = np.array(np.meshgrid(*[np.arange(s) for s in SHAPE], indexing="ij"))
    return nib.affines.apply_affine(AFFINE, grid.reshape(3, -1).T.astype(float))


def silent_counts(coordinates, image_ids, centres):
    """Per voxel, how many coordinate studies reported nothing within RADIUS."""
    counts = np.zeros(len(centres))
    per_study = {sid: g[["x", "y", "z"]].astype(float).to_numpy()
                 for sid, g in coordinates.groupby("study_id") if sid not in image_ids}
    for points in per_study.values():
        if not len(points):
            counts += 1.0            # reported nothing anywhere: silent everywhere
            continue
        near = spatial.cKDTree(points).query_ball_point(centres, RADIUS, return_length=True)
        counts += np.asarray(near) == 0
    return counts, per_study


if __name__ == "__main__":
    centres = voxel_centres()
    print(f"global null, {N_STUDIES} studies, {N_IMAGES} images, {len(centres)} voxels, "
          f"silence radius {RADIUS:.0f} mm\n", flush=True)
    print(f"{'sim':>4} {'coord st':>9} {'peaks':>6} {'silent range':>13} "
          f"{'rho(|g|, silent)':>17} {'p':>10}")
    rhos = []
    for sim in range(N_SIMS):
        collection = create_effect_size_coordinate_studyset(
            [(0, 0, 0)], effect_sizes=0.0, n_studies=N_STUDIES, sample_size=(20, 40),
            seed=41000 + sim, simulate_field=True, n_image_studies=N_IMAGES,
            image_dir=tempfile.mkdtemp(), noise_extent=EXTENT, field_zooms=ZOOMS,
            blob_fwhm=10.0,
            threshold_z=[THRESHOLDS[i % len(THRESHOLDS)] for i in range(N_STUDIES)],
        )
        # A study with no image carries NaN in the `g` column, not a path.
        image_rows = [r for r in collection.images.itertuples() if isinstance(r.g, str)]
        image_ids = {r.study_id for r in image_rows}
        counts, per_study = silent_counts(collection.coordinates, image_ids, centres)
        values = [np.abs(np.asarray(nib.load(r.g).dataobj, dtype=float).ravel())
                  for r in image_rows]
        if not values:
            print(f"{sim:4d}   no image maps in this realisation")
            continue
        magnitude = np.mean(values, axis=0)
        use = np.isfinite(magnitude) & np.isfinite(counts)
        # Spearman: the question is monotone association, and g has heavy tails at n = 20.
        rho, p = stats.spearmanr(magnitude[use], counts[use])
        rhos.append(float(rho))
        print(f"{sim:4d} {len(per_study):9d} {sum(len(v) for v in per_study.values()):6d} "
              f"{f'{counts.min():.0f}-{counts.max():.0f}':>13} {rho:+17.4f} {p:10.2e}",
              flush=True)

    if rhos:
        arr = np.array(rhos)
        print(f"\nmean rho {arr.mean():+.4f}, all {'negative' if (arr < 0).all() else 'mixed'}, "
              f"over {len(arr)} realisations")
        print("A non-zero rho means the fixed silence pattern and the moving image values were "
              "never\nindependent, so the observed map is scored on a pairing the null never "
              "offers it. The null\nis not exchangeable, and the familywise rate is the "
              "statistic that notices.")
