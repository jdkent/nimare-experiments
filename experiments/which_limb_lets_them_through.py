"""The guard mitigates but does not control. Which of its two conditions is in the wrong place?

`fwe_after_the_guard.py` measured the two-foci cell after the max-statistic guard was installed:

    foci/study  guard fired  voxel FWE  rejected|kept  n kept
             2        0.783      0.050          0.231      13

The aggregate looks perfect and is not error control. A withheld fit cannot reject, so the
overall rate is just `(1 - 0.783) * 0.231 = 0.050` -- it lands on nominal by arithmetic. Among the
fits the guard *passes*, rejection is 0.231 against a nominal 0.05, exact binomial p = 0.0245 with
a 95% lower bound of 0.066. So the mechanism is right and the threshold is not.

The guard refuses only when a null is **both** sparse (fewer distinct attained maxima than
`0.10 * n_iters`, i.e. under 20 of 200) **and** narrow (coefficient of variation under 0.05). The
docstring's reason for the conjunction is that either alone can be low benignly -- "a coarse but
wide distribution still separates the observed value from the bulk, and a fine but narrow one can
arise when every study reports many foci of similar size".

That reasoning is testable and this is the test. For every fit, record both statistics and whether
it rejected. Then the passed fits that rejected can be attributed:

  * if they are sparse-but-wide, the "coarse but wide still separates" argument is wrong and
    sparseness should refuse on its own -- a defensible position on resolution grounds, since a
    p-value at alpha cannot be read off a distribution attaining fewer than 1/alpha values;
  * if they are fine-but-narrow, the CV limb is the one that is too permissive;
  * if they are neither, the two statistics do not capture what makes these nulls degenerate and
    a different diagnostic is needed rather than a different threshold.

Guessing between those three and tuning a constant to a thirteen-fit measurement is how a guard
ends up fitted to its own test set. This records the joint distribution instead.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from joblib import Parallel, delayed
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.meta.cbma.effectsize import (
    _null_maxima_diagnostics,
    _MIN_NULL_MAXIMA_DISTINCT_FRACTION,
    _MIN_NULL_MAXIMA_CV,
)

N_SIMS = int(os.environ.get("NSIMS", 80))
N_ITERS = int(os.environ.get("NITERS", 200))
ALPHA = 0.05
TRUTH = (0.0, 0.0, 0.0)


def build_mask():
    shape, step = (25, 25, 25), 4.0
    affine = np.eye(4)
    affine[:3, :3] *= step
    affine[:3, 3] = -step * (np.array(shape) - 1) / 2
    return nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)


def one(seed, n_noise_foci):
    mask = build_mask()
    studyset = create_effect_size_coordinate_studyset(
        [TRUTH], effect_sizes=0.0, n_studies=20, sample_size=(10, 1000),
        tau=0.0, prevalence=0.0, seed=seed, n_noise_foci=n_noise_foci, noise_extent=40.0)
    est = CBES(fwhm=10.0, mask=mask, peak_bias="per-study", n_iters=N_ITERS,
               seed=seed, cluster_threshold=None)
    result = est.fit(studyset)
    if np.all(result.get_map("p", return_type="array").ravel() == 1.0):
        return None
    maps, _, _ = est.correct_fwe_montecarlo(result, vfwe_only=True)
    logp = maps["logp_level-voxel"].ravel()
    covered = result.get_map("n_studies", return_type="array").ravel() > 0
    withheld = bool(np.all(logp[covered] == 0.0))
    rejected = bool((logp[covered] >= -np.log10(ALPHA)).any())

    null = est.null_distributions_.get("values_level-voxel")
    if null is None:
        for key, value in est.null_distributions_.items():
            if "voxel" in key and np.ndim(value) == 1 and np.size(value) > 10:
                null = value
                break
    if null is None:
        return None
    usable, n_distinct, cv = _null_maxima_diagnostics(null)
    return withheld, rejected, int(n_distinct), float(cv), int(np.size(null))


if __name__ == "__main__":
    print(f"{N_SIMS} simulations, 20 studies, N 10-1000, global null, {N_ITERS} permutations")
    print(f"guard refuses when distinct < {_MIN_NULL_MAXIMA_DISTINCT_FRACTION:.2f} * n_iters "
          f"AND cv < {_MIN_NULL_MAXIMA_CV}\n")
    for n_foci in (2, 6):
        rows = [r for r in Parallel(n_jobs=4)(
            delayed(one)(s, n_foci) for s in range(N_SIMS)) if r is not None]
        wh = np.array([r[0] for r in rows], dtype=bool)
        rej = np.array([r[1] for r in rows], dtype=bool)
        nd = np.array([r[2] for r in rows], dtype=float)
        cv = np.array([r[3] for r in rows], dtype=float)
        size = rows[0][4]
        limit = max(2, int(np.ceil(_MIN_NULL_MAXIMA_DISTINCT_FRACTION * size)))
        sparse = nd < limit
        narrow = cv < _MIN_NULL_MAXIMA_CV
        kept = ~wh
        print(f"--- {n_foci} foci per study ({len(rows)} fits, null size {size}, "
              f"sparse means distinct < {limit}) ---")
        print(f"  guard fired {wh.mean():.3f}   overall rejection {rej.mean():.3f}")
        print(f"  distinct maxima: median {np.median(nd):.0f}  range {nd.min():.0f}-{nd.max():.0f}")
        print(f"  cv:              median {np.median(cv):.3f}  range {cv.min():.3f}-{cv.max():.3f}")
        for name, mask_ in (("sparse & narrow (refused)", sparse & narrow),
                            ("sparse, not narrow", sparse & ~narrow),
                            ("narrow, not sparse", ~sparse & narrow),
                            ("neither", ~sparse & ~narrow)):
            n = int(mask_.sum())
            if not n:
                print(f"  {name:28s} n=0")
                continue
            print(f"  {name:28s} n={n:3d}  rejected {rej[mask_].mean():.3f}  "
                  f"median distinct {np.median(nd[mask_]):.0f}  median cv {np.median(cv[mask_]):.3f}")
        if kept.any():
            print(f"  among KEPT fits: n={int(kept.sum())}  rejected {rej[kept].mean():.3f}  "
                  f"median distinct {np.median(nd[kept]):.0f}  median cv {np.median(cv[kept]):.3f}")
        print(flush=True)
    print("The cell that both passes the guard and rejects is the one to read. Its median")
    print("distinct count and cv say which limb to move, and the six-foci rows say what moving")
    print("it would cost in legitimate refusals.")
