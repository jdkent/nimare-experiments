"""If the `se` inflation is the pi/mu trade-off, the interval on `pi*mu` should be the sound one.

The suspicion, from the arm table: a textbook all-image inverse-variance fit reads se/sd 1.04 to
1.24, but the same 20 images under the zero-inflated mixture read far worse -- with the reporting
indicator structurally empty, because every study has an image. If that holds, the inflation is
not the censoring term at all. It is the cost of profiling out the prevalence, and it should be
worst exactly where `pi` and `mu` are least separable: at a voxel with little or no effect, where
"a small effect in every study" and "a large effect in a few" fit equally well.

That has a consequence worth more than the diagnosis. The ridge runs along curves of roughly
constant `pi*mu`, so the product should be far better determined than either factor. `g_marginal`
and `se_marginal` are already emitted, so this is a direct test:

  * if se/sd for `g` is inflated where the truth is weak while se/sd for `g_marginal` is near 1
    there, the ridge is the mechanism and `g_marginal` is the map to put an interval on;
  * if both are inflated equally, the ridge is not the mechanism and something else is.

Scored against the field simulator's own signal. Note the two maps have *different* estimands:
`g` targets `mu` and `g_marginal` targets `pi*mu`, so each is scored against its own truth --
`pi` is 1 by construction here, since every study carries the effect, which makes them coincide
and keeps the comparison honest.
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
N_REPS = int(os.environ.get("NREPS", 32))

grid = np.stack(np.indices(SHAPE), -1) * ZOOMS - ZOOMS * HALF
sd = BLOB / (2 * np.sqrt(2 * np.log(2)))
volume = np.zeros(SHAPE)
for focus, effect in zip(FOCI, EFFECTS):
    volume += effect * np.exp(
        -((grid - np.asarray(focus, float)) ** 2).sum(-1) / (2 * sd ** 2))
TRUTH = MASKER.transform(nib.Nifti1Image(volume.astype(np.float32), AFFINE)).ravel()

STRATA = [("null", TRUTH < 0.02)]
for focus, effect in zip(FOCI, EFFECTS):
    ijk = np.round(np.linalg.inv(AFFINE) @ np.array(list(focus) + [1.0]))[:3].astype(int)
    hit = np.zeros(SHAPE, bool)
    hit[tuple(ijk)] = True
    STRATA.append((f"g={effect}", MASKER.transform(
        nib.Nifti1Image(hit.astype(np.int16), AFFINE)).ravel() > 0))

ARMS = [("2 images, 18 tables", 2), ("20 images, indicator empty", 20)]


def calibration(values, errors, truth):
    """se/sd, se/rmse and the pieces, for one map over one stratum."""
    ok = np.isfinite(errors) & (errors > 0)
    spread = float(np.mean(np.std(values, axis=0, ddof=1)))
    reported = float(np.mean(errors[ok])) if ok.any() else np.nan
    rmse = float(np.sqrt(np.mean((values - truth) ** 2)))
    return spread, reported, reported / spread, reported / rmse


if __name__ == "__main__":
    print(f"four foci at {EFFECTS}, {N_REPS} replications, 20 studies, tau=0.1, true pi = 1\n")
    head = (f"{'arm':>28} {'stratum':>8} | {'sd':>7} {'se':>7} {'se/sd':>6} | "
            f"{'sd_m':>7} {'se_m':>7} {'se/sd_m':>8} | {'pi':>5}")
    print(head)
    for label, n_images in ARMS:
        gs, ses, ms, sems, pis = [], [], [], [], []
        for rep in range(N_REPS):
            ss = create_effect_size_coordinate_studyset(
                FOCI, effect_sizes=EFFECTS, n_studies=20, sample_size=(20, 40), tau=0.1,
                seed=7000 + rep, simulate_field=True, n_image_studies=n_images,
                image_dir=tempfile.mkdtemp(), noise_extent=EXTENT, field_zooms=ZOOMS,
                blob_fwhm=BLOB)
            try:
                res = CBES(mask=MASK, null_method="none",
                           threshold="reporting_threshold").fit(ss)
            except ValueError:
                continue
            gs.append(res.get_map("g", return_type="array").ravel())
            ses.append(res.get_map("se", return_type="array").ravel())
            ms.append(res.get_map("g_marginal", return_type="array").ravel())
            sems.append(res.get_map("se_marginal", return_type="array").ravel())
            pis.append(res.get_map("prevalence", return_type="array").ravel())
        if len(gs) < min(5, N_REPS):
            print(f"{label:>28} {'--':>8}  too few usable replications")
            continue
        G, SE, M, SEM, PI = (np.array(a) for a in (gs, ses, ms, sems, pis))
        for name, sel in STRATA:
            s1, r1, k1, _ = calibration(G[:, sel], SE[:, sel], TRUTH[sel])
            s2, r2, k2, _ = calibration(M[:, sel], SEM[:, sel], TRUTH[sel])
            print(f"{label:>28} {name:>8} | {s1:7.4f} {r1:7.4f} {k1:6.2f} | "
                  f"{s2:7.4f} {r2:7.4f} {k2:8.2f} | {np.mean(PI[:, sel]):5.3f}")
        print()
