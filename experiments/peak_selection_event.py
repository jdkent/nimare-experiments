"""Is a reported peak's height a truncated normal, or a peak height?

The truncation term added to the censored likelihood asks: "what mean explains a value that
had to exceed c?" and answers it with the truncated-normal MLE. That is the right question
only if the selection event were a single draw clearing a threshold. It was not: the value is
the height of a *local maximum of a smooth field*, selected as the largest thing in its
neighbourhood, so it clears c far more often than one draw would and its excess over c is
governed by the field's smoothness rather than by the mean.

If that is what is happening, the truncated-normal MLE will sit well below the true mean, and
it will do so by more the higher the threshold -- because the peak-height distribution tends
to ``u + Exp(u)``, which forgets the signal.

This decides whether the fix is to re-locate the truncation (evaluate it at the peak's own
mean rather than the kernel voxel's) or to replace the selection event entirely.
"""
import numpy as np
from scipy import stats, optimize
from scipy.ndimage import gaussian_filter, maximum_filter

rng = np.random.default_rng(0)
SHAPE = (48, 48, 48)
SMOOTH = 2.0           # voxels; ~5 voxel FWHM
N_FIELDS = 300


def field(mu, sd=1.0):
    raw = gaussian_filter(rng.standard_normal(SHAPE), SMOOTH)
    raw /= raw.std()
    return mu + sd * raw


def peaks(volume, cut):
    mag = np.abs(volume)
    hit = (mag == maximum_filter(mag, size=3)) & (mag >= cut)
    return volume[hit]


def truncated_normal_mle(values, cut, sd):
    """MLE of the mean of N(mu, sd^2) observed only where |x| > cut."""
    def neg_ll(mu):
        surv = 1.0 - (stats.norm.cdf((cut - mu) / sd) - stats.norm.cdf((-cut - mu) / sd))
        surv = max(surv, 1e-12)
        return -(stats.norm.logpdf(values, mu, sd).sum() - len(values) * np.log(surv))
    return optimize.minimize_scalar(neg_ll, bounds=(-6.0, 6.0), method="bounded").x


print(f"{'mu':>5} {'cut':>5} {'peaks':>7} {'mean peak':>10} {'draw-trunc MLE':>15} "
      f"{'single-draw mean':>17} {'single-draw MLE':>16}")
for mu in (0.0, 0.5, 1.0, 2.0):
    for cut in (2.5, 3.29, 4.0):
        found, drawn = [], []
        for _ in range(N_FIELDS):
            vol = field(mu)
            found.append(peaks(vol, cut))
            flat = vol.ravel()
            drawn.append(flat[np.abs(flat) >= cut])
        found = np.concatenate(found)
        drawn = np.concatenate(drawn)
        if found.size < 30 or drawn.size < 30:
            print(f"{mu:5.1f} {cut:5.2f} {found.size:7d}   too few")
            continue
        sub = found if found.size < 20000 else rng.choice(found, 20000, replace=False)
        sub_d = drawn if drawn.size < 20000 else rng.choice(drawn, 20000, replace=False)
        print(f"{mu:5.1f} {cut:5.2f} {found.size:7d} {found.mean():10.3f} "
              f"{truncated_normal_mle(sub, cut, 1.0):15.3f} "
              f"{drawn.mean():17.3f} {truncated_normal_mle(sub_d, cut, 1.0):16.3f}")
