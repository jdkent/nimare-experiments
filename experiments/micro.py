"""A one-dimensional testbed for intuitions about thresholded, selectively reported data.

Every claim this project has argued over -- whether peak height or peak count carries the signal,
whether a truncated normal is the right selection model, whether the effect scale and the
reporting threshold are separately identifiable -- is a claim about the statistics of a smooth
random field that has been thresholded and reduced to a few reported maxima. None of it is about
brains. In one dimension, on a few thousand grid points, each of those questions answers in
milliseconds instead of the tens of minutes a whole-brain fit costs, and the truth is known
exactly rather than estimated from held-out subjects.

The trade is real and worth stating: a 1-D field has different topology from a 3-D one, so the
*numbers* here will not match a brain. Resel counts, peak densities and the Euler characteristic
all differ. What transfers is the *direction and mechanism* of an effect. So the rule for using
this file is: form the intuition here, then confirm the one that survives on real data. To keep
that honest, `calibration()` reproduces a result already established in 3-D, and anything this
testbed says should be discounted if that check drifts.

Run with no arguments to execute every test.
"""
import sys
import numpy as np
from scipy import optimize, stats
from scipy.ndimage import gaussian_filter1d

LENGTH = 4096
SMOOTH = 4.0          # grid points; the field's FWHM is about 2.35 times this
rng = np.random.default_rng(0)


#: Standard deviation of the smoothing filter's output for unit-variance white input. Computed
#: once from the filter's own coefficients rather than from each realisation, because dividing a
#: realisation by its *own* empirical standard deviation gives the field a Student-t marginal on
#: roughly the effective degrees of freedom. That is a real defect -- the tail is what every
#: question here depends on -- though it turned out not to be the one that broke the calibration.
_impulse = np.zeros(LENGTH)
_impulse[LENGTH // 2] = 1.0
NOISE_SD = float(np.sqrt((gaussian_filter1d(_impulse, SMOOTH, mode="wrap") ** 2).sum()))


def noise(n=1):
    """Unit-variance smooth Gaussian noise, `n` independent realisations."""
    raw = gaussian_filter1d(rng.standard_normal((n, LENGTH)), SMOOTH, axis=-1, mode="wrap")
    return raw / NOISE_SD


def truth_field(n_bumps=12, width=30.0):
    """A latent effect field: bumps of widely varying height, so strata exist to compare."""
    x = np.arange(LENGTH)
    field = np.zeros(LENGTH)
    centres = rng.choice(LENGTH, n_bumps, replace=False)
    heights = np.exp(rng.uniform(np.log(0.1), np.log(1.2), n_bumps))
    for c, h in zip(centres, heights):
        d = np.minimum(np.abs(x - c), LENGTH - np.abs(x - c))
        field += h * np.exp(-0.5 * (d / width) ** 2)
    return field


def peaks_above(z, cut):
    """Indices of local maxima of |z| at or above `cut`."""
    mag = np.abs(z)
    interior = (mag[1:-1] >= mag[:-2]) & (mag[1:-1] >= mag[2:]) & (mag[1:-1] >= cut)
    return np.flatnonzero(interior) + 1


def report(z, cut, per_cluster=False):
    """What a paper prints: every local maximum, or one representative per supra-threshold run.

    The distinction matters more than the threshold does. Reporting every maximum thins the
    pattern as the cut rises but leaves each reported value free to vary. Reporting one
    representative per cluster pins every value near the top of whatever cluster it came from,
    which is what destroys the magnitude channel -- and it is what cluster-extent reporting does.
    """
    if not per_cluster:
        return peaks_above(z, cut)
    above = np.abs(z) >= cut
    if not above.any():
        return np.array([], dtype=int)
    edges = np.flatnonzero(np.diff(above.astype(np.int8)))
    starts = np.r_[0 if above[0] else [], edges[::2] + 1 if above[0] else edges[0::2] + 1]
    bounds = np.flatnonzero(np.diff(np.r_[False, above, False].astype(np.int8)))
    out = []
    for lo, hi in zip(bounds[0::2], bounds[1::2]):
        out.append(lo + int(np.argmax(np.abs(z[lo:hi]))))
    return np.array(out, dtype=int)


# --------------------------------------------------------------- calibration against 3-D
def calibration():
    """Reproduce, in 1-D, a result already established in 3-D: the truncated-normal MLE

    applied to *peaks* under-recovers a small mean and over-recovers a large one, while the same
    MLE applied to ordinary supra-threshold *draws* is roughly right. In 3-D at a cut of 3.29 the
    peak MLE returned 0.26 for a true 0.5. If this testbed cannot show the same direction, nothing
    else it says should be believed.
    """
    print("calibration: truncated-normal MLE, peaks versus draws (3-D gave 0.26 for a true 0.5)")
    print(f"  {'true mean':>10} {'cut':>5} {'mean peak':>10} {'MLE on peaks':>13} "
          f"{'MLE on draws':>13} {'eff. draws':>11}")
    for mu in (0.5, 2.0):
        for cut in (2.5, 3.29):
            found, drawn = [], []
            for _ in range(2500):
                z = mu + noise()[0]
                # Signed, not absolute. The likelihood below is for values drawn from
                # N(m, 1) and kept when |x| >= cut, so folding the negative exceedances onto
                # the positive side inflates the mean -- by about 0.19 at a true 0.5 and a cut
                # of 3.29, where 2.8% of exceedances are negative, and by nothing at a true 2.0
                # where essentially none are. That asymmetry was the whole discrepancy.
                found.append(z[peaks_above(z, cut)])
                drawn.append(z[np.abs(z) >= cut])
            found = np.concatenate(found)
            drawn = np.concatenate(drawn)

            def mle(values):
                def neg(m):
                    surv = max(stats.norm.sf(cut - m) + stats.norm.cdf(-cut - m), 1e-12)
                    return -(stats.norm.logpdf(values, m, 1.0).sum()
                             - values.size * np.log(surv))
                return optimize.minimize_scalar(neg, bounds=(-4, 6), method="bounded").x

            sub = lambda a: a if a.size < 8000 else rng.choice(a, 8000, replace=False)
            effective = drawn.size / (SMOOTH * np.sqrt(4 * np.pi))
            print(f"  {mu:10.1f} {cut:5.2f} {found.mean():10.3f} "
                  f"{mle(sub(found)):13.3f} {mle(sub(drawn)):13.3f} {effective:11.0f}")
    print()


# --------------------------------------------------------- T1: which channel carries the signal
def channel_sweep():
    """Height versus count as the reporting threshold slides, so density varies continuously.

    The 3-D runs sampled two regimes and found the answer flipped between them. Here the cut moves
    smoothly, so the crossover itself is visible and can be located rather than bracketed.
    """
    print("T1: which channel predicts the truth, by reporting rule and density")
    print(f"  {'rule':>14} {'cut':>5} {'foci/study':>11} {'r, height':>10} {'r, count':>9} "
          f"{'winner':>8}")
    mu = truth_field()
    n_studies, n_subj = 20, 25
    smooth_kernel = 12.0
    for per_cluster in (False, True):
      for cut in (2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0):
        height_num = np.zeros(LENGTH)
        height_den = np.zeros(LENGTH)
        count = np.zeros(LENGTH)
        total = 0
        for _ in range(n_studies):
            z = mu * np.sqrt(n_subj) + noise()[0]
            idx = report(z, cut, per_cluster)
            total += idx.size
            if not idx.size:
                continue
            spike = np.zeros(LENGTH)
            spike[idx] = np.abs(z[idx]) / np.sqrt(n_subj)       # reported effect size
            reach = np.zeros(LENGTH)
            reach[idx] = 1.0
            height_num += gaussian_filter1d(spike, smooth_kernel, mode="wrap")
            height_den += gaussian_filter1d(reach, smooth_kernel, mode="wrap")
            near = np.zeros(LENGTH)
            near[idx] = 1.0
            count += (gaussian_filter1d(near, smooth_kernel, mode="wrap") > 1e-6)
        label = "one per cluster" if per_cluster else "every maximum"
        use = height_den > 1e-9
        if use.sum() < 100:
            print(f"  {label:>14} {cut:5.2f} {total / n_studies:11.1f}   too few reported")
            continue
        height = np.zeros(LENGTH)
        height[use] = height_num[use] / height_den[use]
        r_h = stats.pearsonr(height[use], mu[use])[0]
        r_c = stats.pearsonr(count[use], mu[use])[0]
        print(f"  {label:>14} {cut:5.2f} {total / n_studies:11.1f} {r_h:10.3f} {r_c:9.3f} "
              f"{('count' if r_c > r_h else 'height'):>8}")
      print()


# ------------------------------------------- T3: is the effect scale identifiable from the cut?
def identifiability():
    """Can the effect scale and the reporting threshold be told apart from coordinates alone?

    The claim borrowed from ecology was that they cannot -- that intensity and detectability are
    confounded, and unbiased data is needed to separate them. That claim was overstated when
    repeated from secondary sources, so it is tested here directly.

    The shape of the effect is treated as known and only its scale `a` is unknown, alongside the
    cut `u`. If the profile of the log-likelihood in `a` is flat, the two are confounded; if it
    has a well-separated minimum, they are not. Then one "image" study -- a direct noisy
    observation of the field, which is what an unthresholded map is -- is added, and the profile
    re-examined.
    """
    print("T3: is the effect scale identifiable from coordinates alone?")
    shape = truth_field()
    shape = shape / shape.max()
    a_true, u_true, n_subj, n_studies = 0.8, 3.0, 25, 30

    tables = []
    for _ in range(n_studies):
        z = a_true * shape * np.sqrt(n_subj) + noise()[0]
        tables.append(peaks_above(z, u_true))

    def neg_loglik(a, u):
        """Poisson likelihood for the peak pattern, intensity from the expected-maxima rate.

        The rate of local maxima above `u` in a smooth Gaussian field with local mean `m` is
        proportional to the zero-mean rate evaluated at `u - m`, which is the approximation the
        random field literature uses at high thresholds.
        """
        m = a * shape * np.sqrt(n_subj)
        rate = np.exp(-0.5 * np.maximum(u - m, 0.0) ** 2) + 1e-12
        rate = rate / rate.sum()
        total = sum(t.size for t in tables)
        ll = 0.0
        for t in tables:
            ll += np.log(rate[t]).sum()
        return -(ll - total * np.log(1.0))

    print(f"  {'':>16} {'a = 0.4':>9} {'a = 0.6':>9} {'a = 0.8':>9} {'a = 1.0':>9} "
          f"{'a = 1.5':>9}")
    grid_a = (0.4, 0.6, 0.8, 1.0, 1.5)
    row = []
    for a in grid_a:
        # Profile over the cut: the best the model can do at this scale.
        best = min(neg_loglik(a, u) for u in np.linspace(1.5, 5.0, 36))
        row.append(best)
    base = min(row)
    print(f"  {'coordinates only':>16} " + " ".join(f"{v - base:9.1f}" for v in row))

    # One image study: a direct noisy observation of the same field.
    image = a_true * shape + noise()[0] / np.sqrt(n_subj)
    sigma = 1.0 / np.sqrt(n_subj)
    row_img = []
    for a in grid_a:
        best = min(neg_loglik(a, u) for u in np.linspace(1.5, 5.0, 36))
        gauss = 0.5 * (((image - a * shape) / sigma) ** 2).sum()
        row_img.append(best + gauss)
    base_img = min(row_img)
    print(f"  {'plus one image':>16} " + " ".join(f"{v - base_img:9.1f}" for v in row_img))
    print("\n  Values are the increase in negative log-likelihood from the best scale; the true"
          "\n  scale is 0.8. A flat row means the scale is not identifiable from that data.\n")


if __name__ == "__main__":
    which = sys.argv[1:] or ["calibration", "channel_sweep", "identifiability"]
    for name in which:
        globals()[name]()
