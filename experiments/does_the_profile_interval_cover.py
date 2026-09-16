"""Does the profile interval cover where the Wald interval does not?

The Wald interval is wrong in both directions and no iteration count fixes it: se/sd is 0.80 at
the default and 2.55 at convergence, for the same data at prevalence 0.6. The arithmetic is exact
-- the observed information matches a finite-difference Hessian of the same likelihood to four
decimals -- so the fault is that curvature at one point does not describe how far mu moves along a
near-flat ridge from one collection to the next.

``interval="profile"`` is already implemented and does not use the curvature at all. It inverts
the likelihood ratio, which is parameterisation-invariant and has better higher-order properties
than a Wald interval. Nobody has measured its coverage in the regime the estimator exists for. If
it is calibrated where Wald is not, the fix is a default change rather than new theory.

**Kill conditions, stated before running.**

  1. If the profile interval's coverage is no closer to 0.95 than Wald's at prevalence 0.6, it is
     not the answer and the interval problem needs the harder route.
  2. If the profile interval is calibrated only by being enormously wide -- median width more than
     three times Wald's while covering -- it is honest but useless, and that has to be said rather
     than reported as a win. Coverage without width is on the metrics-that-lie list for a reason.

Environment knobs: REPS, VOXELS, NIMG, NTAB, MU, PREV, TAU, MAXITER.
"""
import os
import sys

import numpy as np

sys.path.insert(0, "/home/user/NiMARE")
from nimare.meta.cbma.effectsize import CBES  # noqa: E402

REPS = int(os.environ.get("REPS", "200"))
VOXELS = int(os.environ.get("VOXELS", "32"))
N_IMG = int(os.environ.get("NIMG", "2"))
N_TAB = int(os.environ.get("NTAB", "20"))
MU = float(os.environ.get("MU", "0.5"))
PREV = float(os.environ.get("PREV", "0.6"))
TAU = float(os.environ.get("TAU", "0.15"))
MAX_ITER = int(os.environ.get("MAXITER", "25"))

N_SUBJ = 20
VAR_G = 1.0 / N_SUBJ
CUT = 3.09 / np.sqrt(N_SUBJ)
SD_A = np.sqrt(VAR_G + TAU**2)
SD_0 = np.sqrt(VAR_G)


def draw(rng, width):
    active_img = rng.random((N_IMG, width)) < PREV
    g_img = np.where(
        active_img,
        MU + SD_A * rng.standard_normal((N_IMG, width)),
        SD_0 * rng.standard_normal((N_IMG, width)),
    )
    active_tab = rng.random((N_TAB, width)) < PREV
    latent = np.where(
        active_tab,
        MU + SD_A * rng.standard_normal((N_TAB, width)),
        SD_0 * rng.standard_normal((N_TAB, width)),
    )
    return g_img, np.where(np.abs(latent) >= CUT, -1, 1)


def fit(g_img, indicator, interval):
    width = g_img.shape[1]
    n = N_IMG + N_TAB
    weights = np.zeros((n, width))
    weights[:N_IMG] = 1.0
    values = np.zeros((n, width))
    values[:N_IMG] = g_img
    variances = np.full((n, width), VAR_G)
    full_indicator = np.zeros((n, width))
    full_indicator[N_IMG:] = indicator

    est = CBES(selection_model="zero-inflated", interval=interval, tau2_method="none",
               max_iter=MAX_ITER)
    mu, _pi, se, _sem, _share, lower, upper = est._fit_chunk(
        weights=weights,
        g_obs=values,
        var_obs=variances,
        indicator=full_indicator,
        tau2=np.full(width, TAU**2),
        null_var=np.full(n, VAR_G),
        cutoffs=np.full(n, CUT),
        start=g_img.mean(axis=0),
    )
    return mu, se, lower, upper


def main():
    rng = np.random.default_rng(0)
    mus, ses, lowers, uppers = [], [], [], []
    for _ in range(REPS):
        g_img, indicator = draw(rng, VOXELS)
        mu, se, lower, upper = fit(g_img, indicator, "profile")
        mus.append(mu)
        ses.append(se)
        lowers.append(lower)
        uppers.append(upper)

    mu = np.concatenate(mus)
    se = np.concatenate(ses)
    lower = np.concatenate(lowers)
    upper = np.concatenate(uppers)
    err = mu - MU

    # Compare on the *same* voxels. A profile bound is infinite wherever the data do not reject
    # pi = 0, so scoring profile on its 55% and Wald on its 93% would compare two different
    # populations of voxels -- and the profile's are selected for being informative, which is
    # exactly the direction that would flatter it.
    bounded = np.isfinite(lower) & np.isfinite(upper) & np.isfinite(se)
    wald_all = np.isfinite(se)
    wald_cov_all = np.abs(err[wald_all]) <= 1.96 * se[wald_all]
    wald_cov = np.abs(err[bounded]) <= 1.96 * se[bounded]
    wald_width = 2 * 1.96 * se[bounded]
    prof_cov = (lower[bounded] <= MU) & (MU <= upper[bounded])
    prof_width = (upper - lower)[bounded]

    print(f"true mu {MU}, prevalence {PREV}, tau {TAU}, {N_IMG} images, {N_TAB} tables, "
          f"max_iter {MAX_ITER}")
    print(f"{mu.size} fits; sd of the estimate {err.std():.4f}")
    print()
    print(f"{'interval':10s} {'usable':>8s} {'coverage':>9s} {'median width':>13s} "
          f"{'mean width':>11s}")
    print(f"{'wald (all)':10s} {wald_all.mean():8.3f} {wald_cov_all.mean():9.3f} "
          f"{'':13s} {'':11s}")
    print(f"{'wald':10s} {bounded.mean():8.3f} {wald_cov.mean():9.3f} "
          f"{np.median(wald_width):13.4f} {wald_width.mean():11.4f}")
    if bounded.any():
        print(f"{'profile':10s} {bounded.mean():8.3f} {prof_cov.mean():9.3f} "
              f"{np.median(prof_width):13.4f} {prof_width.mean():11.4f}")
    else:
        print(f"{'profile':10s} {0.0:8.3f} {'--':>9s} {'--':>13s} {'--':>11s}")

    print()
    print("Both rows after the first are scored on the voxels where the profile bound is finite,")
    print("so they compare like with like. 'wald (all)' is Wald over every voxel it can score, to")
    print("show how much that subset differs from the whole.")
    if not bounded.any():
        print()
        print("KILL CONDITION: the profile interval is unbounded everywhere here, so it cannot")
        print("be the fix. The interval problem needs the harder route.")
        return
    gap_wald = abs(wald_cov.mean() - 0.95)
    gap_prof = abs(prof_cov.mean() - 0.95)
    print()
    if gap_prof >= gap_wald:
        print(f"KILL CONDITION 1 MET: profile is no closer to 0.95 ({gap_prof:.3f} against "
              f"{gap_wald:.3f}).")
    elif np.median(prof_width) > 3 * np.median(wald_width):
        print(f"KILL CONDITION 2 MET: calibrated only by being {np.median(prof_width) / np.median(wald_width):.1f} "
              "times wider.")
    else:
        print(f"Profile is closer to nominal ({gap_prof:.3f} against {gap_wald:.3f}) at "
              f"{np.median(prof_width) / np.median(wald_width):.2f} times the width.")


if __name__ == "__main__":
    main()
