"""Does `coordinate_share` predict where the interval fails? If so it is actionable, not decorative.

The `se` runs 1.55 to 3.74 times the estimator's own spread across replications, and the excess is
located in the censoring term -- switching the silence off drops it from 1.78 to 1.29 where the
effect is. That was measured as a whole-map figure, which gives a reader nothing to act on.

`coordinate_share` is the per-voxel fraction of the information about `g` that came from the
coordinate indicators rather than the images' values. If the excess really lives in the censoring
term, then **se/sd should be worse where the share is high and near 1 where it is low** -- and a
reader gets a rule: trust the interval where the share is small.

That is a falsifiable prediction about a quantity the estimator now reports, so it is worth the
run either way. If se/sd is flat in the share, the diagnostic does not speak to the interval and
the docstring should not imply it does.

Scored against the field simulator's own signal, known exactly. `se/sd` needs the spread across
replications at each voxel, so the arms are per-voxel across seeds rather than per-seed.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tempfile
import numpy as np
import nibabel as nib
from nilearn.maskers import NiftiMasker
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma.effectsize import CBES

ZOOMS, EXTENT, BLOB = 4.0, 48.0, 10.0
HALF = int(np.ceil(EXTENT / ZOOMS))
SHAPE = (2 * HALF + 1,) * 3
AFFINE = np.array([[ZOOMS, 0, 0, -ZOOMS * HALF], [0, ZOOMS, 0, -ZOOMS * HALF],
                   [0, 0, ZOOMS, -ZOOMS * HALF], [0, 0, 0, 1.0]])
MASK = nib.Nifti1Image(np.ones(SHAPE, np.int32), AFFINE)
MASKER = NiftiMasker(MASK).fit()
FOCI = [(-28, -28, 0), (28, -28, 0), (-28, 28, 0), (28, 28, 0)]
EFFECTS = [0.2, 0.4, 0.6, 0.8]
N_REPS = int(os.environ.get("NREPS", 20))

grid = np.stack(np.indices(SHAPE), -1) * ZOOMS - ZOOMS * HALF
sd = BLOB / (2 * np.sqrt(2 * np.log(2)))
volume = np.zeros(SHAPE)
for focus, effect in zip(FOCI, EFFECTS):
    volume += effect * np.exp(
        -((grid - np.asarray(focus, float)) ** 2).sum(-1) / (2 * sd ** 2))
TRUTH = MASKER.transform(nib.Nifti1Image(volume.astype(np.float32), AFFINE)).ravel()

if __name__ == "__main__":
    gs, ses, shares = [], [], []
    for rep in range(N_REPS):
        ss = create_effect_size_coordinate_studyset(
            FOCI, effect_sizes=EFFECTS, n_studies=20, sample_size=(20, 40), tau=0.1,
            seed=2000 + rep, simulate_field=True, n_image_studies=2,
            image_dir=tempfile.mkdtemp(), noise_extent=EXTENT, field_zooms=ZOOMS,
            blob_fwhm=BLOB)
        res = CBES(mask=MASK, null_method="none", threshold="reporting_threshold").fit(ss)
        gs.append(res.get_map("g", return_type="array").ravel())
        ses.append(res.get_map("se", return_type="array").ravel())
        shares.append(res.get_map("coordinate_share", return_type="array").ravel())
    G, SE, SHARE = np.array(gs), np.array(ses), np.array(shares)

    share = SHARE.mean(0)
    spread = G.std(0, ddof=1)
    reported = np.where(np.isfinite(SE) & (SE > 0), SE, np.nan)
    reported = np.nanmean(reported, axis=0)
    usable = np.isfinite(reported) & (spread > 1e-6)

    print(f"four foci at {EFFECTS}, {N_REPS} replications, 20 studies, 2 images")
    print(f"coordinate_share over usable voxels: median {np.median(share[usable]):.3f}, "
          f"p95 {np.percentile(share[usable], 95):.3f}\n")
    edges = [0.0, 0.05, 0.10, 0.25, 0.50, 1.01]
    print(f"{'coordinate_share':>18} {'voxels':>8} {'se/sd':>7} {'mean se':>9} "
          f"{'sd of g':>9} {'bias':>8}")
    for lo, hi in zip(edges[:-1], edges[1:]):
        band = usable & (share >= lo) & (share < hi)
        if band.sum() < 20:
            continue
        ratio = float(np.mean(reported[band]) / np.mean(spread[band]))
        bias = float(np.mean(G.mean(0)[band] - TRUTH[band]))
        print(f"  [{lo:.2f}, {hi:.2f}){'':>4} {int(band.sum()):8d} {ratio:7.2f} "
              f"{np.mean(reported[band]):9.4f} {np.mean(spread[band]):9.4f} {bias:+8.4f}")

    top = usable & (share >= np.percentile(share[usable], 90))
    bottom = usable & (share <= np.percentile(share[usable], 10))
    print(f"\ntop decile of share:    se/sd {np.mean(reported[top]) / np.mean(spread[top]):.2f}")
    print(f"bottom decile of share: se/sd "
          f"{np.mean(reported[bottom]) / np.mean(spread[bottom]):.2f}")
