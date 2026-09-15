"""Is the `se` over-stated, or is the point estimate rigid? Every se/sd on record used sd(|g|).

`g` is a signed inverse-variance mean, so `np.abs` before taking a spread across replications is
a transformation, not a no-op: at a null voxel sd(|g|) is smaller than sd(g), which inflates the
ratio. Every se/sd figure on record -- 1.55 to 3.74, and the flat-in-the-share table -- was
computed that way.

Two things are therefore untested:

  * **sd(|g|) against sd(g).** If the corrected spread is larger, the over-statement shrinks.
  * **variance against total error.** se/sd asks whether the interval matches the estimator's
    *spread*. A reader asks whether it covers the *truth*. Where bias exceeds the spread -- it
    did in the top share band, +0.19 against 0.16 -- those are different questions, and se/rmse
    is the one that bears on coverage.

A rigid over-shrunk estimator has a small spread by construction, so a large se/sd can mean the
point estimate is not moving rather than that the interval is too wide. se/rmse separates them.

**Calibration arm first.** Arm A gives every study an image and switches the selection model off,
which is a textbook local inverse-variance random-effects fit: se/sd must come out near 1. If it
does not, the harness is wrong and nothing else here can be believed.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tempfile
import numpy as np
import nibabel as nib
from scipy import stats
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

# Strata by the truth itself, so bias and spread are read where the signal is known.
STRATA = [("null", TRUTH < 0.02)]
for focus, effect in zip(FOCI, EFFECTS):
    ijk = np.round(np.linalg.inv(AFFINE) @ np.array(list(focus) + [1.0]))[:3].astype(int)
    hit = np.zeros(SHAPE, bool)
    hit[tuple(ijk)] = True
    STRATA.append((f"g={effect}", MASKER.transform(
        nib.Nifti1Image(hit.astype(np.int16), AFFINE)).ravel() > 0))

ARMS = [
    ("A textbook IVW (20 images, no selection)", dict(n_images=20), dict(selection_model="none")),
    ("B images only (2 images, no selection)", dict(n_images=2), dict(selection_model="none")),
    ("C shipped (2 images, zero-inflated)", dict(n_images=2), dict()),
    # D isolates the cost of the mixture itself. With every study carrying an image there are no
    # coordinate-only studies, so the reporting indicator is structurally empty and the only
    # difference from A is that the prevalence is fitted and profiled out.
    ("D 20 images, zero-inflated", dict(n_images=20), dict()),
]

if __name__ == "__main__":
    print(f"four foci at {EFFECTS}, {N_REPS} replications, 20 studies, tau=0.1\n")
    print("  se/sd(g)   reported se over the spread of the SIGNED estimate across replications")
    print("  se/sd(|g|) the same with the absolute value taken first -- the figure on record")
    print("  se/rmse    reported se over total error (bias and spread together)")
    print("  cov        share of replications whose t interval covers the truth\n")
    head = (f"{'arm':>42} {'stratum':>8} {'bias':>7} {'sd(g)':>7} {'sd|g|':>7} "
            f"{'se':>7} {'se/sd':>6} {'se/sd|':>7} {'se/rmse':>8} {'cov':>5}")
    print(head)
    for label, build, options in ARMS:
        gs, ses, dofs = [], [], []
        for rep in range(N_REPS):
            ss = create_effect_size_coordinate_studyset(
                FOCI, effect_sizes=EFFECTS, n_studies=20, sample_size=(20, 40), tau=0.1,
                seed=3000 + rep, simulate_field=True, n_image_studies=build["n_images"],
                image_dir=tempfile.mkdtemp(), noise_extent=EXTENT, field_zooms=ZOOMS,
                blob_fwhm=BLOB)
            try:
                res = CBES(mask=MASK, null_method="none",
                           threshold="reporting_threshold", **options).fit(ss)
            except ValueError:
                continue
            gs.append(res.get_map("g", return_type="array").ravel())
            ses.append(res.get_map("se", return_type="array").ravel())
            dofs.append(res.get_map("dof", return_type="array").ravel())
        if len(gs) < min(5, N_REPS):
            print(f"{label:>42} {'--':>8}  too few usable replications")
            continue
        G, SE, DOF = np.array(gs), np.array(ses), np.array(dofs)
        for name, sel in STRATA:
            g, se, dof = G[:, sel], SE[:, sel], DOF[:, sel]
            truth = TRUTH[sel]
            ok = np.isfinite(se) & (se > 0)
            bias = float(np.mean(g - truth))
            spread = float(np.mean(np.std(g, axis=0, ddof=1)))
            spread_abs = float(np.mean(np.std(np.abs(g), axis=0, ddof=1)))
            reported = float(np.mean(se[ok])) if ok.any() else np.nan
            rmse = float(np.sqrt(np.mean((g - truth) ** 2)))
            crit = stats.t.isf(0.025, np.clip(dof, 1e-6, None))
            half = crit * se
            fine = np.isfinite(half) & ok
            cov = float(np.mean((np.abs(g - truth) <= half)[fine])) if fine.any() else np.nan
            print(f"{label:>42} {name:>8} {bias:+7.3f} {spread:7.4f} {spread_abs:7.4f} "
                  f"{reported:7.4f} {reported / spread:6.2f} "
                  f"{reported / spread_abs:7.2f} {reported / rmse:8.2f} {cov:5.2f}")
        print()
