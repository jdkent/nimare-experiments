"""Is "censoring, not imputation" a real advantage, or a slogan?

The estimator's docstring asserts that non-reporting is censoring and that nothing should be
imputed. That is a claim about a competitor, and until now nothing had measured the competitor.

MetaNSUE (Albajes-Eizagirre, Solanes & Radua, Stat Methods Med Res 2018,
doi:10.1177/0962280218811349) states the same premise CBES rests on -- a study that reports "not
significant" without an effect size cannot be dropped and cannot be entered as zero -- and
resolves the resulting bound by *multiple imputation*: draw plausible values inside the bound,
meta-analyse each completed dataset, pool with Rubin's rules. SDM-PSI carries that algorithm into
neuroimaging. CBES instead writes the bound into the likelihood and maximises it by EM.

Both are principled. Which is better is an empirical question with four possible answers, and
three of them require the docstring to change:

  * the censored likelihood wins on both point estimate and interval -- keep the line
  * they tie on the point estimate and the censored one gives tighter honest intervals -- narrow
    the claim to the interval
  * they tie outright -- the line is rhetoric and should go
  * imputation wins -- the model should change, not the docstring

The censored arm here is not a reimplementation: it calls ``CBES._fit_chunk`` on hand-built
arrays, so whatever ships is what is measured. Only the imputation arm is written here.

Calibration, per the protocol. The full-scale published-pain bed (8 splits, rr_pain.log) says
that switching the silences off raises the mean error by 0.067 and the error at the strongest
decile by 0.132, both p < 0.001 -- that is, dropping the silent studies *overestimates*. This bed
must reproduce that sign, with the drop-silent arm biased high and the censored arm closer to the
truth. It is a 1-D-free scalar bed, so absolute magnitudes are not expected to transfer; #71
records that the scalar bed reads se/sd 0.74 where the field bed reads 1.83, so intervals here
are read relatively (arm against arm in the same bed) and never as absolute coverage.

Environment knobs: REPS, VOXELS, NIMG, NTAB, MU, TAU, MIMP (imputations), TITER (imputation
outer iterations).
"""
import os
import sys

import numpy as np
from scipy.special import ndtr, ndtri

sys.path.insert(0, "/home/user/NiMARE")
from nimare.meta.cbma.effectsize import CBES  # noqa: E402

REPS = int(os.environ.get("REPS", "200"))
VOXELS = int(os.environ.get("VOXELS", "64"))
N_IMG = int(os.environ.get("NIMG", "2"))
N_TAB = int(os.environ.get("NTAB", "20"))
MU = float(os.environ.get("MU", "0.5"))
TAU = float(os.environ.get("TAU", "0.15"))
M_IMPUTE = int(os.environ.get("MIMP", "20"))
T_ITER = int(os.environ.get("TITER", "6"))

#: Hedges' g sampling variance for a one-sample design at this many subjects, and the reporting
#: cut in the same units. Fixed across studies so that the only thing separating the arms is how
#: each treats the bound -- a spread of cutoffs would also separate them, but through
#: identifiability rather than through the censoring machinery.
N_SUBJ = 20
VAR_G = 1.0 / N_SUBJ
CUTOFF_G = 3.09 / np.sqrt(N_SUBJ)


def draw(rng, width):
    """One synthetic voxel set: image values, and reporting indicators from the table studies."""
    sd_img = np.sqrt(VAR_G + TAU**2)
    g_img = rng.normal(MU, sd_img, size=(N_IMG, width))

    sd_tab = np.sqrt(VAR_G + TAU**2)
    latent = rng.normal(MU, sd_tab, size=(N_TAB, width))
    reported = np.abs(latent) >= CUTOFF_G
    indicator = np.where(reported, -1, 1).astype(np.int8)
    return g_img, indicator


def censored_arm(g_img, indicator):
    """CBES's own EM, called directly on the arrays. Returns (mu, se)."""
    width = g_img.shape[1]
    est = CBES(selection_model="none", interval="wald", tau2_method="none")
    # Unit weights, because that is what the estimator supplies: "Every image contributes weight
    # 1" (_accumulate). Passing inverse-variance weights here leaves the point estimate alone --
    # they scale the whole score, so the root does not move -- but it multiplies the information
    # and makes the interval look far worse than it is.
    weights = np.zeros((N_IMG + N_TAB, width))
    weights[:N_IMG] = 1.0
    values = np.zeros((N_IMG + N_TAB, width))
    values[:N_IMG] = g_img
    variances = np.full((N_IMG + N_TAB, width), VAR_G)
    full_indicator = np.zeros((N_IMG + N_TAB, width), dtype=np.int8)
    full_indicator[N_IMG:] = indicator

    mu, _pi, se, _sem, _share, _lo, _hi = est._fit_chunk(
        weights=weights,
        g_obs=values,
        var_obs=variances,
        indicator=full_indicator,
        tau2=np.full(width, TAU**2),
        null_var=np.full(N_IMG + N_TAB, VAR_G),
        cutoffs=np.full(N_IMG + N_TAB, CUTOFF_G),
        # The shipped code starts the EM at the naive pooled estimate (``start=fit["g"]``), not
        # at zero. It matters: where a study reported, log(1 - P(silent)) is convex, so at mu = 0
        # the total curvature can be positive, the Newton step is refused, and the fit never
        # leaves the start value. Starting anywhere sensible is enough; starting at zero is not.
        start=g_img.mean(axis=0),
    )
    return mu, se


def _truncated_normal(rng, mean, sd, low, high):
    """Inverse-CDF draw from N(mean, sd^2) restricted to [low, high], elementwise."""
    alpha = ndtr((low - mean) / sd)
    beta = ndtr((high - mean) / sd)
    u = alpha + rng.random(mean.shape) * np.clip(beta - alpha, 1e-12, None)
    return mean + sd * ndtri(np.clip(u, 1e-12, 1.0 - 1e-12))


def _impute_outside(rng, mean, sd, cut):
    """Draw from N(mean, sd^2) restricted to |x| >= cut: pick a tail by its mass, then draw."""
    mass_low = ndtr((-cut - mean) / sd)
    mass_high = 1.0 - ndtr((cut - mean) / sd)
    total = np.clip(mass_low + mass_high, 1e-12, None)
    take_high = rng.random(mean.shape) < (mass_high / total)
    big = np.full_like(mean, np.inf)
    low = np.where(take_high, cut, -big)
    high = np.where(take_high, big, -cut)
    return _truncated_normal(rng, mean, sd, low, high)


def imputation_arm(rng, g_img, indicator):
    """MetaNSUE-style multiple imputation within the bounds. Returns (mu, se).

    Imputes every table study, not only the silent ones. That is the honest translation of
    MetaNSUE to this setting: CBES discards reported peak heights, so a study that reported is
    censored too, just on the other side of the cut. Giving the imputation arm the heights it
    would have in a table-based meta-analysis would be comparing two different inputs.
    """
    width = g_img.shape[1]
    sd = np.sqrt(VAR_G + TAU**2)
    precision_img = N_IMG / (VAR_G + TAU**2)
    mu_hat = g_img.mean(axis=0)

    silent = indicator > 0
    for _ in range(T_ITER):
        draws = np.empty((M_IMPUTE, width))
        for m in range(M_IMPUTE):
            mean = np.broadcast_to(mu_hat, indicator.shape)
            inside = _truncated_normal(
                rng, mean, sd, np.full(indicator.shape, -CUTOFF_G),
                np.full(indicator.shape, CUTOFF_G),
            )
            outside = _impute_outside(rng, mean, sd, CUTOFF_G)
            completed = np.where(silent, inside, outside)
            # Fixed-variance inverse-variance pool over images plus completed table studies.
            total = g_img.sum(axis=0) + completed.sum(axis=0)
            draws[m] = total / (N_IMG + N_TAB)
        mu_hat = draws.mean(axis=0)

    within = 1.0 / (precision_img + N_TAB / (VAR_G + TAU**2))
    between = draws.var(axis=0, ddof=1)
    total_var = within + (1.0 + 1.0 / M_IMPUTE) * between
    return mu_hat, np.sqrt(total_var)


def oracle_arm(g_img, indicator):
    """The censored-normal MLE by brute-force grid search. Returns (mu, se).

    Both real arms target this. It exists so that a disagreement between them can be attributed
    rather than argued about: whichever arm is further from the grid is the one with the defect.
    Variance is known here, so the likelihood is the image values' normal density times the
    indicator probabilities, and a grid fine enough to beat the Monte Carlo noise is cheap.
    """
    width = g_img.shape[1]
    sd = np.sqrt(VAR_G + TAU**2)
    grid = np.linspace(-1.5, 2.0, 1401)[:, None]           # (grid, 1)
    n_sil = (indicator > 0).sum(axis=0)[None, :]
    n_rep = (indicator < 0).sum(axis=0)[None, :]

    silent_prob = np.clip(
        ndtr((CUTOFF_G - grid) / sd) - ndtr((-CUTOFF_G - grid) / sd), 1e-300, None
    )
    ll = n_sil * np.log(silent_prob) + n_rep * np.log(np.clip(1.0 - silent_prob, 1e-300, None))
    ll = ll - 0.5 * (((g_img[:, None, :] - grid[None]) / sd) ** 2).sum(axis=0)

    best = ll.argmax(axis=0)
    mu = grid[:, 0][best]
    # Curvature of the profile at the maximum, by second difference on the grid.
    step = grid[1, 0] - grid[0, 0]
    inner = np.clip(best, 1, grid.size - 2)
    second = (ll[inner + 1, np.arange(width)] - 2 * ll[inner, np.arange(width)]
              + ll[inner - 1, np.arange(width)]) / step**2
    se = np.where(second < 0, 1.0 / np.sqrt(-np.minimum(second, -1e-12)), np.inf)
    return mu, se


def drop_silent_arm(g_img, indicator):
    """The estimator MetaNSUE's abstract warns about: keep the reporters, drop the silent."""
    width = g_img.shape[1]
    sd = np.sqrt(VAR_G + TAU**2)
    reported = indicator < 0
    # A reported study contributes the conditional mean of |g| given it cleared the cut, which is
    # the best a table gives you without heights.
    mass = np.clip(1.0 - ndtr((CUTOFF_G - MU) / sd) + ndtr((-CUTOFF_G - MU) / sd), 1e-12, None)
    conditional = MU + sd * (
        np.exp(-0.5 * ((CUTOFF_G - MU) / sd) ** 2) - np.exp(-0.5 * ((-CUTOFF_G - MU) / sd) ** 2)
    ) / np.sqrt(2 * np.pi) / mass
    n_rep = reported.sum(axis=0)
    total = g_img.sum(axis=0) + conditional * n_rep
    return total / (N_IMG + n_rep), np.full(width, np.nan)


def zero_fill_arm(g_img, indicator):
    """The other estimator the abstract warns about: enter every silence as a zero."""
    width = g_img.shape[1]
    sd = np.sqrt(VAR_G + TAU**2)
    mass = np.clip(1.0 - ndtr((CUTOFF_G - MU) / sd) + ndtr((-CUTOFF_G - MU) / sd), 1e-12, None)
    conditional = MU + sd * (
        np.exp(-0.5 * ((CUTOFF_G - MU) / sd) ** 2) - np.exp(-0.5 * ((-CUTOFF_G - MU) / sd) ** 2)
    ) / np.sqrt(2 * np.pi) / mass
    reported = (indicator < 0).sum(axis=0)
    total = g_img.sum(axis=0) + conditional * reported
    return total / (N_IMG + N_TAB), np.full(width, np.nan)


def main():
    rng = np.random.default_rng(0)
    arms = {"grid MLE (oracle)": [], "censored (CBES)": [],
            "imputation (MetaNSUE-style)": [], "drop silent": [], "zero fill": []}
    errors = {name: [] for name in arms}
    widths = {name: [] for name in arms}
    covered = {name: [] for name in arms}

    for _ in range(REPS):
        g_img, indicator = draw(rng, VOXELS)
        results = {
            "grid MLE (oracle)": oracle_arm(g_img, indicator),
            "censored (CBES)": censored_arm(g_img, indicator),
            "imputation (MetaNSUE-style)": imputation_arm(rng, g_img, indicator),
            "drop silent": drop_silent_arm(g_img, indicator),
            "zero fill": zero_fill_arm(g_img, indicator),
        }
        for name, (mu, se) in results.items():
            errors[name].append(mu - MU)
            widths[name].append(2 * 1.96 * se)
            covered[name].append(np.abs(mu - MU) <= 1.96 * se)

    print(f"true mu {MU}, tau {TAU}, {N_IMG} image studies, {N_TAB} table studies")
    print(f"cut {CUTOFF_G:.3f} in g units; {(REPS * VOXELS)} fits per arm")
    reported_share = None
    print()
    print(f"{'arm':30s} {'bias':>8s} {'sd':>8s} {'rmse':>8s} {'CI width':>10s} {'coverage':>9s}")
    for name in arms:
        err = np.concatenate([np.asarray(e).ravel() for e in errors[name]])
        wid = np.concatenate([np.asarray(w).ravel() for w in widths[name]])
        cov = np.concatenate([np.asarray(c).ravel() for c in covered[name]])
        finite = np.isfinite(wid)
        width_text = f"{np.nanmean(wid[finite]):10.4f}" if finite.any() else f"{'--':>10s}"
        cov_text = f"{np.nanmean(cov[finite]):9.3f}" if finite.any() else f"{'--':>9s}"
        print(
            f"{name:30s} {err.mean():+8.4f} {err.std():8.4f} "
            f"{np.sqrt((err**2).mean()):8.4f} {width_text} {cov_text}"
        )

    print()
    print("calibration: the full-scale pain bed says dropping the silences overestimates.")
    drop_bias = np.concatenate([np.asarray(e).ravel() for e in errors["drop silent"]]).mean()
    cen_bias = np.concatenate([np.asarray(e).ravel() for e in errors["censored (CBES)"]]).mean()
    verdict = "reproduced" if drop_bias > cen_bias and drop_bias > 0 else "FAILED"
    print(f"  drop-silent bias {drop_bias:+.4f} against censored {cen_bias:+.4f}: {verdict}")
    if verdict == "FAILED":
        print("  stop: this bed cannot be trusted on anything else until that is understood.")
    return reported_share


if __name__ == "__main__":
    main()
