"""Is the zero-inflated fit the MLE, and is its interval the right one?

The estimator's Notes say ``se/sd`` runs 1.2 to 2.2 against a known truth and attribute it to
"the cost of profiling out a prevalence the indicator only partly identifies". That explanation
is not obviously right. Profiling out a nuisance parameter is supposed to *give* the correct
standard error for the parameter of interest, not inflate it by a factor of two. Either the
explanation is wrong or the Schur complement is.

A first attempt compared the fit against the argmax of a 2-D grid over (mu, pi). That oracle was
worse than the estimator it was judging -- bias +0.19 against -0.01, sd 0.48 against 0.40 -- which
is the bed's fault, not the estimator's. The mixture likelihood has a near-flat ridge: many
(mu, pi) pairs explain the data about equally well, so a global argmax wanders while an EM started
at the pooled estimate settles on a stable point nearby. "The MLE" is not a well-defined target
here the way it is without the mixture.

So ask the two questions separately, neither of which needs a unique maximum.

  1. *Is the fit at a maximum?* Compare the log-likelihood the estimator achieves against the best
     the grid can find. A fit that is not below the grid is converged, whichever optimum it chose.
  2. *Is the interval right?* Evaluate the exact 2-D Hessian at the estimator's own point, by
     finite differences on the analytic log-likelihood rather than on the coarse grid, and take
     the (mu, mu) element of its inverse. That is what the Schur complement claims to equal, and
     comparing them there is independent of which optimum was found.

The same approach found two bugs in the non-mixture path (see censoring_versus_imputation.py), so
the calibration check is that this bed reproduces those on ``selection_model="none"`` before being
believed on the mixture: with CALIBRATE=1 it refits the non-mixture case and asserts agreement.

Environment knobs: REPS, VOXELS, NIMG, NTAB, MU, PREV, TAU, CALIBRATE, MAXITER.
"""
import os
import sys

import numpy as np
from scipy.special import ndtr

sys.path.insert(0, "/home/user/NiMARE")
from nimare.meta.cbma.effectsize import CBES  # noqa: E402

REPS = int(os.environ.get("REPS", "200"))
VOXELS = int(os.environ.get("VOXELS", "32"))
N_IMG = int(os.environ.get("NIMG", "2"))
N_TAB = int(os.environ.get("NTAB", "20"))
MU = float(os.environ.get("MU", "0.5"))
PREV = float(os.environ.get("PREV", "0.6"))
TAU = float(os.environ.get("TAU", "0.15"))
CALIBRATE = os.environ.get("CALIBRATE", "1") == "1"
MAX_ITER = int(os.environ.get("MAXITER", "25"))

N_SUBJ = 20
VAR_G = 1.0 / N_SUBJ
CUTOFF_G = 3.09 / np.sqrt(N_SUBJ)

MU_GRID = np.linspace(-0.5, 1.6, 211)
PI_GRID = np.linspace(0.02, 0.98, 97)


def draw(rng, width, prevalence):
    """Image values and reporting indicators, with a fraction ``1 - prevalence`` of null studies.

    A null study has no effect at this voxel at all: its image value is centred on zero and its
    table reports only if noise alone cleared the cut.
    """
    sd_active = np.sqrt(VAR_G + TAU**2)
    sd_null = np.sqrt(VAR_G)

    active_img = rng.random((N_IMG, width)) < prevalence
    g_img = np.where(
        active_img,
        rng.normal(MU, sd_active, size=(N_IMG, width)),
        rng.normal(0.0, sd_null, size=(N_IMG, width)),
    )

    active_tab = rng.random((N_TAB, width)) < prevalence
    latent = np.where(
        active_tab,
        rng.normal(MU, sd_active, size=(N_TAB, width)),
        rng.normal(0.0, sd_null, size=(N_TAB, width)),
    )
    indicator = np.where(np.abs(latent) >= CUTOFF_G, -1, 1).astype(np.int8)
    return g_img, indicator


def cbes_arm(g_img, indicator, selection_model):
    """``_fit_chunk`` on hand-built arrays, so what ships is what is measured."""
    width = g_img.shape[1]
    est = CBES(
        selection_model=selection_model, interval="wald", tau2_method="none", max_iter=MAX_ITER
    )
    n = N_IMG + N_TAB
    weights = np.zeros((n, width))
    weights[:N_IMG] = 1.0                       # "Every image contributes weight 1"
    values = np.zeros((n, width))
    values[:N_IMG] = g_img
    variances = np.full((n, width), VAR_G)
    full_indicator = np.zeros((n, width), dtype=np.int8)
    full_indicator[N_IMG:] = indicator

    mu, pi, se, _se_marginal, _share, surviving, *_ = est._fit_chunk(
        weights=weights,
        g_obs=values,
        var_obs=variances,
        indicator=full_indicator,
        tau2=np.full(width, TAU**2),
        null_var=np.full(n, VAR_G),
        cutoffs=np.full(n, CUTOFF_G),
        start=g_img.mean(axis=0),
    )
    return mu, pi, se, surviving


def _log_likelihood_surface(g_img, indicator):
    """``ll[i, j, voxel]`` over the (mu, pi) grid, for the zero-inflated mixture."""
    sd_active = np.sqrt(VAR_G + TAU**2)
    sd_null = np.sqrt(VAR_G)
    mu = MU_GRID[:, None, None]
    pi = PI_GRID[None, :, None]

    silent_active = np.clip(
        ndtr((CUTOFF_G - mu) / sd_active) - ndtr((-CUTOFF_G - mu) / sd_active), 1e-300, None
    )
    silent_null = float(ndtr(CUTOFF_G / sd_null) - ndtr(-CUTOFF_G / sd_null))

    n_sil = (indicator > 0).sum(axis=0)[None, None, :]
    n_rep = (indicator < 0).sum(axis=0)[None, None, :]
    mix_silent = pi * silent_active + (1.0 - pi) * silent_null
    ll = n_sil * np.log(np.clip(mix_silent, 1e-300, None))
    ll = ll + n_rep * np.log(np.clip(1.0 - mix_silent, 1e-300, None))

    for study in range(N_IMG):
        value = g_img[study][None, None, :]
        density_active = np.exp(-0.5 * ((value - mu) / sd_active) ** 2) / sd_active
        density_null = np.exp(-0.5 * (value / sd_null) ** 2) / sd_null
        ll = ll + np.log(
            np.clip(pi * density_active + (1.0 - pi) * density_null, 1e-300, None)
        )
    return ll


def log_likelihood_at(mu, pi, g_img, indicator):
    """Analytic log-likelihood of the mixture, elementwise over voxels."""
    sd_active = np.sqrt(VAR_G + TAU**2)
    sd_null = np.sqrt(VAR_G)
    pi = np.clip(pi, 1e-9, 1.0 - 1e-9)

    silent_active = np.clip(
        ndtr((CUTOFF_G - mu) / sd_active) - ndtr((-CUTOFF_G - mu) / sd_active), 1e-300, None
    )
    silent_null = float(ndtr(CUTOFF_G / sd_null) - ndtr(-CUTOFF_G / sd_null))
    mix_silent = pi * silent_active + (1.0 - pi) * silent_null

    total = (indicator > 0).sum(axis=0) * np.log(np.clip(mix_silent, 1e-300, None))
    total = total + (indicator < 0).sum(axis=0) * np.log(
        np.clip(1.0 - mix_silent, 1e-300, None)
    )
    for study in range(N_IMG):
        value = g_img[study]
        density_active = np.exp(-0.5 * ((value - mu) / sd_active) ** 2) / sd_active
        density_null = np.exp(-0.5 * (value / sd_null) ** 2) / sd_null
        total = total + np.log(
            np.clip(pi * density_active + (1.0 - pi) * density_null, 1e-300, None)
        )
    return total


def exact_se_at(mu, pi, g_img, indicator, h_mu=1e-3, h_pi=1e-3):
    """Schur-complement standard error from the exact 2-D Hessian at a given (mu, pi).

    Independent of which optimum the estimator chose, so it isolates the information from the
    optimisation. Pi is nudged off the boundary before differencing, since the Hessian does not
    exist at 0 or 1.
    """
    pi = np.clip(pi, 2 * h_pi, 1.0 - 2 * h_pi)

    def ll(dm, dp):
        return log_likelihood_at(mu + dm * h_mu, pi + dp * h_pi, g_img, indicator)

    centre = ll(0, 0)
    h_mm = (ll(1, 0) - 2 * centre + ll(-1, 0)) / h_mu**2
    h_pp = (ll(0, 1) - 2 * centre + ll(0, -1)) / h_pi**2
    h_mp = (ll(1, 1) - ll(1, -1) - ll(-1, 1) + ll(-1, -1)) / (4 * h_mu * h_pi)

    schur = np.where(h_pp < 0, -h_mm + h_mp**2 / np.where(h_pp < 0, h_pp, -1.0), -h_mm)
    return np.where(schur > 0, 1.0 / np.sqrt(np.maximum(schur, 1e-12)), np.inf)


def best_grid_log_likelihood(g_img, indicator):
    """The highest log-likelihood a 2-D grid can find, per voxel."""
    return _log_likelihood_surface(g_img, indicator).reshape(-1, g_img.shape[1]).max(axis=0)


def summarise(label, mu, se, truth):
    err = (mu - truth).ravel()
    finite = np.isfinite(se.ravel())
    width = 2 * 1.96 * se.ravel()[finite]
    covered = np.abs(err[finite]) <= 1.96 * se.ravel()[finite]
    return (
        f"{label:26s} {err.mean():+8.4f} {err.std():8.4f} "
        f"{np.sqrt((err**2).mean()):8.4f} {width.mean():9.4f} {covered.mean():9.3f} "
        f"{se.ravel()[finite].mean() / err.std():8.2f}"
    )


def main():
    if CALIBRATE:
        rng = np.random.default_rng(11)
        g_img, indicator = draw(rng, 64, 1.0)
        mu_c, _pi, se_c, _surviving = cbes_arm(g_img, indicator, "none")
        # At prevalence 1 the mixture is inert, so the non-mixture fit must equal the 1-D MLE.
        sd = np.sqrt(VAR_G + TAU**2)
        grid = MU_GRID[:, None]
        silent = np.clip(
            ndtr((CUTOFF_G - grid) / sd) - ndtr((-CUTOFF_G - grid) / sd), 1e-300, None
        )
        ll = (indicator > 0).sum(axis=0)[None, :] * np.log(silent)
        ll = ll + (indicator < 0).sum(axis=0)[None, :] * np.log(
            np.clip(1.0 - silent, 1e-300, None)
        )
        ll = ll - 0.5 * (((g_img[:, None, :] - grid[None]) / sd) ** 2).sum(axis=0)
        gap = np.abs(mu_c - MU_GRID[ll.argmax(axis=0)]).max()
        print(f"calibration: non-mixture fit against the 1-D grid MLE, worst gap {gap:.4f}")
        if gap > 2 * (MU_GRID[1] - MU_GRID[0]):
            print("  FAILED -- stop, this bed cannot be trusted on the mixture either.")
            return

    rng = np.random.default_rng(0)
    errors, reported_se, exact_se, ll_gap, fitted_pi = [], [], [], [], []
    survivors = []
    for _ in range(REPS):
        g_img, indicator = draw(rng, VOXELS, PREV)
        mu_c, pi_c, se_c, surviving_c = cbes_arm(g_img, indicator, "zero-inflated")
        survivors.append(surviving_c)
        errors.append(mu_c - MU)
        reported_se.append(se_c)
        exact_se.append(exact_se_at(mu_c, pi_c, g_img, indicator))
        ll_gap.append(
            best_grid_log_likelihood(g_img, indicator)
            - log_likelihood_at(mu_c, pi_c, g_img, indicator)
        )
        fitted_pi.append(pi_c)

    err = np.concatenate(errors)
    reported = np.concatenate(reported_se)
    exact = np.concatenate(exact_se)
    gap = np.concatenate(ll_gap)
    pi_hat = np.concatenate(fitted_pi)
    usable = np.isfinite(reported) & np.isfinite(exact) & (exact > 0)

    print()
    print(f"true mu {MU}, prevalence {PREV}, tau {TAU}, {N_IMG} images, {N_TAB} tables")
    print(f"{err.size} fits; {usable.mean():.3f} have a finite se from both routes")
    print()
    print(f"  bias {err.mean():+.4f}   sd {err.std():.4f}   "
          f"rmse {np.sqrt((err**2).mean()):.4f}   "
          f"fitted prevalence median {np.median(pi_hat):.3f} (true {PREV})")
    print()
    print("1. is the fit at a maximum? log-likelihood below the best the grid found:")
    print(f"     median {np.median(gap):+.5f}, 90th pct {np.percentile(gap, 90):+.5f}, "
          f"max {gap.max():+.5f}")
    print(f"     voxels more than 0.01 below the grid: {(gap > 0.01).mean():.4f}")
    print()
    print(f"   se/sd {reported[usable].mean() / err.std():.2f}")
    print("2. is the interval right? reported se against the exact Schur se at the same point:")
    ratio = reported[usable] / exact[usable]
    print(f"     ratio median {np.median(ratio):.4f}, mean {ratio.mean():.4f}, "
          f"10th-90th {np.percentile(ratio, 10):.4f}-{np.percentile(ratio, 90):.4f}")
    print()
    print(f"     se/sd reported {reported[usable].mean() / err.std():.2f}, "
          f"exact {exact[usable].mean() / err.std():.2f}")
    print("     se/sd above 1 with the exact Hessian means the width is the likelihood's,")
    print("     not the estimator's arithmetic.")

    # The diagnostic the algebra predicts: identified_share is 1 - I_mupi^2 / (I_mumu I_pipi),
    # the fraction of the information about mu that survives estimating the prevalence. The Wald
    # interval should fail exactly where it is small. If coverage does not track it, the
    # diagnostic does not earn its place in the estimator.
    share = np.concatenate(survivors)
    covered = np.abs(err) <= 1.96 * reported
    print()
    print("3. does identified_share predict where the Wald interval fails?")
    edges = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0001]
    print(f"     {'identified_share':>18s} {'voxels':>8s} {'coverage':>9s} {'mean se':>9s} "
          f"{'sd of err':>10s}")
    for lo, hi in zip(edges[:-1], edges[1:]):
        band = (share >= lo) & (share < hi) & np.isfinite(reported)
        if band.sum() < 20:
            continue
        print(f"     {f'{lo:.1f}-{hi:.1f}':>18s} {int(band.sum()):8d} "
              f"{covered[band].mean():9.3f} {reported[band].mean():9.4f} "
              f"{err[band].std():10.4f}")


if __name__ == "__main__":
    main()
