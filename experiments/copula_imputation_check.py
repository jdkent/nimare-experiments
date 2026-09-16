"""Check two claims made about how a completed dataset would have to be built.

Both were asserted in conversation before being measured, which is the wrong order.

  1. *Smoothing after imputation breaks the bounds.* Draw independent values inside each voxel's
     censoring interval, then apply a spatial kernel to make the map look like a brain map. The
     claim is that the smoothed values no longer satisfy the intervals they were drawn from, and
     that their variance collapses.
  2. *A Gaussian copula fixes it.* Draw a smooth Gaussian field with the target covariance, map it
     through Phi to uniforms, then through each voxel's truncated-normal inverse CDF. The claim is
     that this satisfies every bound exactly while keeping most of the spatial correlation.

One dimension is enough for both: they are claims about marginals and about a covariance, neither
of which needs three dimensions. What 1-D cannot speak to is how much correlation survives in a
volume, where a kernel of a given FWHM couples many more neighbours -- so read the retained
correlation here as a direction, not a number.
"""
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.special import ndtr, ndtri

LENGTH = 4096
SMOOTH = 4.0
CUT = 0.69
MU = 0.3
SD = 0.27
rng = np.random.default_rng(0)


def smooth_field(n):
    """Unit-variance smooth Gaussian noise."""
    impulse = np.zeros(LENGTH)
    impulse[LENGTH // 2] = 1.0
    scale = np.sqrt((gaussian_filter1d(impulse, SMOOTH, mode="wrap") ** 2).sum())
    raw = gaussian_filter1d(rng.standard_normal((n, LENGTH)), SMOOTH, axis=-1, mode="wrap")
    return raw / scale


def bounds():
    """A silent/reported pattern, and the interval each voxel's completed value must satisfy.

    Silent voxels are bounded to (-CUT, CUT); reported ones to |x| >= CUT, taking the side the
    latent value fell on, which is what a table's sign tells you.
    """
    latent = MU + SD * smooth_field(1)[0]
    silent = np.abs(latent) < CUT
    low = np.where(silent, -CUT, np.where(latent > 0, CUT, -np.inf))
    high = np.where(silent, CUT, np.where(latent > 0, np.inf, -CUT))
    return silent, low, high


def independent_then_smooth(low, high, fwhm_points):
    u = rng.random(LENGTH)
    alpha, beta = ndtr((low - MU) / SD), ndtr((high - MU) / SD)
    draw = MU + SD * ndtri(np.clip(alpha + u * (beta - alpha), 1e-12, 1 - 1e-12))
    return draw, gaussian_filter1d(draw, fwhm_points / 2.355, mode="wrap")


def copula(low, high):
    field = smooth_field(1)[0]
    u = ndtr(field)
    alpha, beta = ndtr((low - MU) / SD), ndtr((high - MU) / SD)
    return MU + SD * ndtri(np.clip(alpha + u * (beta - alpha), 1e-12, 1 - 1e-12))


def violations(x, low, high):
    return float(np.mean((x < low - 1e-9) | (x > high + 1e-9)))


def autocorrelation(x, lag):
    a, b = x[:-lag] - x.mean(), x[lag:] - x.mean()
    return float((a * b).mean() / x.var())


def main():
    silent, low, high = bounds()
    target = smooth_field(1)[0]
    lag = int(round(SMOOTH))
    print(f"{silent.mean():.3f} of voxels silent; target field autocorrelation at lag {lag}: "
          f"{autocorrelation(target, lag):+.3f}")
    print()
    print(f"{'construction':34s} {'bound violations':>17s} {'sd':>8s} {'autocorr':>9s}")

    raw, smoothed = independent_then_smooth(low, high, fwhm_points=SMOOTH * 2.355)
    print(f"{'independent draw, unsmoothed':34s} {violations(raw, low, high):17.4f} "
          f"{raw.std():8.4f} {autocorrelation(raw, lag):9.3f}")
    print(f"{'independent draw, then smoothed':34s} {violations(smoothed, low, high):17.4f} "
          f"{smoothed.std():8.4f} {autocorrelation(smoothed, lag):9.3f}")

    drawn = copula(low, high)
    print(f"{'Gaussian copula':34s} {violations(drawn, low, high):17.4f} "
          f"{drawn.std():8.4f} {autocorrelation(drawn, lag):9.3f}")

    print()
    print("Claim 1 holds if smoothing shows violations above 0 and a smaller sd.")
    print("Claim 2 holds if the copula shows exactly 0 violations and an autocorrelation")
    print("well above the independent draw's.")


if __name__ == "__main__":
    main()
