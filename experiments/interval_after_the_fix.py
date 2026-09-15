"""Re-measure the interval now that the coordinate channel actually contributes.

Every se/sd figure on record was taken when the reporting threshold was compared against
effect sizes on the *z* scale -- eighteen sampling standard deviations out, so every silence
was certain whatever the effect and the censoring term said nothing -- and before the reporting
limb of the indicator existed. So the whole interval section of the CBES docstring describes a
model that no longer exists.

Scored against the field simulator's own signal, which is known exactly and independent of the
peaks. Reported per stratum, because a whole-map figure is dominated by the 99% of voxels whose
truth is near zero and rewards any estimator that shrinks.

  se/sd   reported error over the spread of the estimate across replications. Below 1 is what a
          calibrated *conditional* error should give, the likelihood conditioning on where the
          foci fell while the replication spread is marginal over that.
  cov(t)  coverage of g +/- t(dof) * se, the documented recipe.
  width   half-width of that interval as a fraction of the truth, because coverage without
          width is not a measurement.
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
N_REPS = int(os.environ.get("NREPS", 40))

grid = np.stack(np.indices(SHAPE), -1) * ZOOMS - EXTENT
sd = BLOB / (2 * np.sqrt(2 * np.log(2)))
TRUTH = MASKER.transform(nib.Nifti1Image(
    (TRUE_G * np.exp(-(grid ** 2).sum(-1) / (2 * sd ** 2))).astype(np.float32), AFFINE)).ravel()
FOCUS = int(np.argmax(TRUTH))
STRATA = [("quiet", TRUTH < 0.05), ("middle", (TRUTH >= 0.05) & (TRUTH < 0.25)),
          ("effect", TRUTH >= 0.25)]

ARMS = [
    ("20 studies,  1 image", dict(n_studies=20, n_images=1), {}),
    ("20 studies,  2 images", dict(n_studies=20, n_images=2), {}),
    ("20 studies,  5 images", dict(n_studies=20, n_images=5), {}),
    ("20 studies, 20 images", dict(n_studies=20, n_images=20), {}),
    ("20 studies,  2 images, silence off", dict(n_studies=20, n_images=2),
     dict(selection_model="none")),
    ("20 studies,  2 images, tau 0.3", dict(n_studies=20, n_images=2, tau=0.3), {}),
]

if __name__ == "__main__":
    print(f"true g = {TRUE_G} at the focus, {N_REPS} replications per arm, "
          f"truth known exactly\n")
    header = (f"{'arm':>36} {'stratum':>8} {'bias':>7} {'se/sd':>6} {'dof':>5} "
              f"{'cov(t)':>7} {'width':>6}")
    print(header)
    for label, build, options in ARMS:
        gs, ses, dofs = [], [], []
        for rep in range(N_REPS):
            tmp = tempfile.mkdtemp()
            ss = create_effect_size_coordinate_studyset(
                [(0, 0, 0)], effect_sizes=TRUE_G, sample_size=(20, 40),
                seed=1000 + rep, simulate_field=True, image_dir=tmp,
                noise_extent=EXTENT, field_zooms=ZOOMS, blob_fwhm=BLOB,
                n_studies=build["n_studies"], n_image_studies=build["n_images"],
                tau=build.get("tau", 0.1),
            )
            try:
                res = CBES(mask=MASK, null_method="none",
                           threshold="reporting_threshold", **options).fit(ss)
            except ValueError:
                continue
            gs.append(np.abs(res.get_map("g", return_type="array").ravel()))
            ses.append(res.get_map("se", return_type="array").ravel())
            dofs.append(res.get_map("dof", return_type="array").ravel())
        if len(gs) < 5:
            print(f"{label:>36} {'--':>8}  too few usable replications")
            continue
        G, SE, DOF = np.array(gs), np.array(ses), np.array(dofs)
        for name, sel in STRATA:
            g, se, dof = G[:, sel], SE[:, sel], DOF[:, sel]
            truth = TRUTH[sel]
            bias = float(np.mean(g - truth))
            spread = float(np.mean(np.std(g, axis=0)))
            reported = float(np.mean(se[np.isfinite(se) & (se > 0)])) if np.any(se > 0) else np.nan
            median_dof = float(np.median(dof))
            crit = stats.t.isf(0.025, np.clip(dof, 1e-6, None))
            half = crit * se
            covered = np.abs(g - truth) <= half
            cov = float(np.mean(covered[np.isfinite(half)]))
            width = float(np.mean(half[np.isfinite(half)]) / max(np.mean(truth), 1e-6))
            ratio = reported / spread if spread > 0 else np.nan
            print(f"{label:>36} {name:>8} {bias:+7.3f} {ratio:6.2f} {median_dof:5.1f} "
                  f"{cov:7.2f} {width:6.2f}")
        print()
