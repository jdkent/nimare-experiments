"""An algebraic prediction, tested: heterogeneous studies should separate `pi` from `mu`.

From silences alone the likelihood depends on `(pi, mu)` only through the probability of the
event each pair observed,

    P(pi, mu) = pi S(mu; sigma, c) + (1 - pi) S(0; sigma, c),   S = P(|g| < c | mu),

which is ONE equation in two unknowns -- for studies that share `(sigma, c)`. Verified directly:
along the profile ridge at a quiet voxel `P` is constant to 0.4% while `pi*mu` swings by a factor
of several hundred and changes sign.

But `S` depends on the study's own sampling sd and reporting cut. So `k` *distinct* `(sigma, c)`
configurations give `k` equations in the same two unknowns, and the ridge should collapse once
`k >= 2` with sufficiently different curvature. The prediction is therefore:

    **A collection of studies with widely spread sample sizes and thresholds identifies the
    magnitude; a collection of near-identical studies does not, however many studies it has.**

That is the opposite of the usual intuition, in which heterogeneity is a nuisance.

`interval="profile"` is the probe, because the profile interval on `mu` is bounded exactly when
the likelihood pins `mu` down -- so the *fraction of voxels with a bounded interval* is a direct
identifiability measure rather than a proxy. Also reported: the error in the fitted prevalence,
which is what the extra equations are supposed to buy.

Study count is held fixed at 20 across arms, so nothing here is about having more data. No cap
on peaks anywhere; each study reports whatever clears its own threshold.
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
N_REPS = int(os.environ.get("NREPS", 8))
N_STUDIES = 20

grid = np.stack(np.indices(SHAPE), -1) * ZOOMS - ZOOMS * HALF
sd = BLOB / (2 * np.sqrt(2 * np.log(2)))
volume = np.zeros(SHAPE)
for focus, effect in zip(FOCI, EFFECTS):
    volume += effect * np.exp(
        -((grid - np.asarray(focus, float)) ** 2).sum(-1) / (2 * sd ** 2))
TRUTH = MASKER.transform(nib.Nifti1Image(volume.astype(np.float32), AFFINE)).ravel()
NEAR = TRUTH >= 0.10

# Same number of studies in every arm; only how alike they are changes.
TIGHT_Z = [3.1] * N_STUDIES
SPREAD_Z = list(np.linspace(2.3, 4.5, N_STUDIES))
ARMS = [
    ("alike: n 28-32, one cut", (28, 32), TIGHT_Z),
    ("spread n only: 12-120", (12, 120), TIGHT_Z),
    ("spread cut only: 2.3-4.5", (28, 32), SPREAD_Z),
    ("both spread", (12, 120), SPREAD_Z),
]

if __name__ == "__main__":
    print(f"{N_STUDIES} studies, 2 with images, true prevalence 1.0, {N_REPS} replications")
    print("bounded = fraction of voxels whose profile interval on mu is finite both sides,")
    print("which is where the likelihood identifies the magnitude at all.\n")
    head = (f"{'arm':>27} {'bounded all':>12} {'bounded near':>13} {'pi near':>8} "
            f"{'|g| near':>9} {'truth':>7}")
    print(head)
    for label, sample_size, thresholds in ARMS:
        frac_all, frac_near, pis, mags = [], [], [], []
        for rep in range(N_REPS):
            ss = create_effect_size_coordinate_studyset(
                FOCI, effect_sizes=EFFECTS, n_studies=N_STUDIES, sample_size=sample_size,
                tau=0.1, seed=11000 + rep, simulate_field=True, n_image_studies=2,
                image_dir=tempfile.mkdtemp(), noise_extent=EXTENT, field_zooms=ZOOMS,
                blob_fwhm=BLOB, threshold_z=thresholds)
            try:
                res = CBES(mask=MASK, null_method="none", threshold="reporting_threshold",
                           interval="profile").fit(ss)
            except ValueError:
                continue
            lo = res.get_map("g_lower", return_type="array").ravel()
            hi = res.get_map("g_upper", return_type="array").ravel()
            bounded = np.isfinite(lo) & np.isfinite(hi)
            frac_all.append(bounded.mean())
            frac_near.append(bounded[NEAR].mean())
            pis.append(res.get_map("prevalence", return_type="array").ravel()[NEAR].mean())
            mags.append(
                np.abs(res.get_map("g", return_type="array").ravel())[NEAR].mean())
        if not frac_all:
            print(f"{label:>27}  no usable replications")
            continue
        print(f"{label:>27} {np.mean(frac_all):12.3f} {np.mean(frac_near):13.3f} "
              f"{np.mean(pis):8.3f} {np.mean(mags):9.3f} {TRUTH[NEAR].mean():7.3f}")
