"""If a paper reports several contrasts, does the interval know they share subjects? (task #16)

Real papers report two to four contrasts, each with its own coordinate table, and a collection
assembled from them has several analyses per study. Those analyses are *not* independent: they
come from the same subjects, so their sampling errors are correlated. If the estimator treats
them as independent observations then the effective sample size is inflated and the standard
error is too small -- by roughly the square root of the number of contrasts per study if they
were perfectly correlated, less if partially.

Nothing in CBES appears to account for study identity in its weighting. `n_eff` is Kish's
effective count computed from kernel weights, which measures how unevenly the *spatial* kernel
distributes weight and knows nothing about which analyses came from the same scanner session.
The permutation null shuffles within an analysis, which is also blind to it.

So the question is quantitative: how much does the interval shrink, and how much does coverage
fall, as the same subjects are split into more reported contrasts? The truth is held fixed and
the number of *subjects* is held fixed; only the number of tables they are reported in changes.
A collection of 12 studies reporting one contrast each is compared against 12 studies reporting
two, three and four -- same subjects, same effect, more tables.

If the standard error falls as contrasts are added while the bias stays put, the interval is
counting correlated observations as independent, and coverage will drop the way it does under the
kernel-width sweep and under adding studies: a fixed bias against a shrinking interval.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from joblib import Parallel, delayed
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.utils import mm2vox

SITE = (0.0, 0.0, 0.0)
TRUE_MU = 0.60
TRUE_PI = 1.00
N_SIMS = int(os.environ.get("NSIMS", 40))
N_STUDIES = 12
REPORTING_P = 2.0 * stats.norm.sf(3.2905)
LOCALISATION_SD = 4.0
#: How strongly two contrasts from one study agree. 1.0 would make extra contrasts pure
#: duplication; 0.0 would make them genuinely independent replications of the same subjects,
#: which is impossible. Real contrasts from one session sit somewhere between, and the sweep
#: brackets it rather than guessing.
WITHIN_STUDY_R = (0.9, 0.5)


def build_mask():
    shape, step = (21, 21, 21), 4.0
    affine = np.eye(4)
    affine[:3, :3] *= step
    affine[:3, 3] = -step * (np.array(shape) - 1) / 2
    return nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)


def one(seed, n_contrasts, rho):
    """12 studies of the same size, each reporting `n_contrasts` correlated tables."""
    rng = np.random.default_rng(seed)
    mask = build_mask()
    studies = []
    for k in range(N_STUDIES):
        n = int(rng.integers(20, 41))
        cut = float(stats.t.isf(REPORTING_P / 2.0, n - 1))
        meta = {"sample_sizes": [n], "reporting_threshold": cut}
        # One study-level draw shared by every contrast, plus a contrast-specific part. The
        # subjects are the same, so the shared component is what makes the tables correlated.
        shared = rng.standard_normal()
        analyses = []
        for c in range(n_contrasts):
            own = rng.standard_normal()
            deviate = np.sqrt(rho) * shared + np.sqrt(1.0 - rho) * own
            t = TRUE_MU * np.sqrt(n) + deviate
            points = []
            if rng.random() < TRUE_PI and abs(t) >= cut:
                loc = np.asarray(SITE) + rng.normal(0, LOCALISATION_SD, 3)
                points.append((loc, t))
            loc = rng.uniform(-36, 36, 3)
            points.append((loc, (cut + rng.exponential(1.0 / cut)) * rng.choice([-1.0, 1.0])))
            analyses.append({"id": f"s{k}-{c}", "name": f"{c}", "metadata": meta, "points": [
                {"space": "MNI", "coordinates": [float(v) for v in loc],
                 "values": [{"kind": "T", "value": float(tv)}]} for loc, tv in points]})
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": analyses})
    est = CBES(fwhm=10.0, mask=mask, peak_bias=None, null_method="none",
               threshold="reporting_threshold")
    res = est.fit(Studyset({"id": "c", "name": "c", "studies": studies},
                           target=None, mask=mask))
    at = mm2vox(np.asarray([SITE]), mask.affine)[0]
    pos = int(np.ravel_multi_index(tuple(at), mask.shape))
    return (float(res.get_map("g", return_type="array").ravel()[pos]),
            float(res.get_map("se", return_type="array").ravel()[pos]),
            float(res.get_map("n_eff", return_type="array").ravel()[pos]))


if __name__ == "__main__":
    print(f"{N_STUDIES} studies, same subjects throughout, true mu {TRUE_MU}, "
          f"{N_SIMS} replications")
    print("Only the number of tables the same subjects are reported in changes.\n")
    for rho in WITHIN_STUDY_R:
        print(f"--- contrasts within a study correlated at {rho} ---")
        print(f"{'contrasts':>10s} {'mean g':>8s} {'bias':>8s} {'mean se':>9s} "
              f"{'sd of g':>9s} {'se/sd':>7s} {'n_eff':>7s} {'cover':>7s}")
        for n_contrasts in (1, 2, 3, 4):
            rows = [r for r in Parallel(n_jobs=6)(
                delayed(one)(s, n_contrasts, rho) for s in range(N_SIMS)) if r is not None]
            a = np.array(rows, dtype=float)
            ok = np.isfinite(a).all(axis=1) & (a[:, 1] > 0)
            g, se, neff = a[ok, 0], a[ok, 1], a[ok, 2]
            cover = np.mean((g - 1.96 * se <= TRUE_MU) & (TRUE_MU <= g + 1.96 * se))
            print(f"{n_contrasts:10d} {g.mean():8.3f} {g.mean()-TRUE_MU:+8.3f} "
                  f"{se.mean():9.3f} {g.std(ddof=1):9.3f} "
                  f"{se.mean()/max(g.std(ddof=1),1e-9):7.2f} {neff.mean():7.1f} "
                  f"{cover:7.2f}", flush=True)
        print()
    print("If mean se falls with the contrast count while the bias holds, the interval is")
    print("counting correlated tables as independent studies. n_eff is the place a correction")
    print("would show up if one existed -- it should stop rising once the subjects run out.")
