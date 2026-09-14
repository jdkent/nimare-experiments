"""What does the likelihood at a one-reporting-study voxel actually look like?

The EM does not converge at these voxels -- g moved 0.24 between max_iter 25 and 200, and 67%
of voxels on the pain collection have exactly one reporting study. Two very different diagnoses
fit that: (a) the joint surface has a real peak but the EM crawls toward it along a ridge, in
which case a better optimizer fixes it, or (b) the surface is genuinely flat, in which case no
optimizer fixes anything and the only honest output is an interval.

Built directly from the model: one reporting study with g_1 at the reporting threshold, m silent
studies, and the exact zero-inflated censored likelihood on a grid.
"""
import numpy as np
from scipy.special import ndtr

SIGMA = 0.25          # within-study sd of g
CUTOFF = 0.55         # |g| below which a study stays silent
G_OBS = 0.70          # the one reported value, just over the threshold


def log_likelihood(mu, pi, m_silent, g_obs=G_OBS):
    """Zero-inflated: with prob pi the effect is present at mu, otherwise it is absent."""
    present = np.exp(-0.5 * ((g_obs - mu) / SIGMA) ** 2) / (SIGMA * np.sqrt(2 * np.pi))
    absent = np.exp(-0.5 * (g_obs / SIGMA) ** 2) / (SIGMA * np.sqrt(2 * np.pi))
    reported = pi * present + (1 - pi) * absent

    silent_present = ndtr((CUTOFF - mu) / SIGMA) - ndtr((-CUTOFF - mu) / SIGMA)
    silent_absent = ndtr(CUTOFF / SIGMA) - ndtr(-CUTOFF / SIGMA)
    silent = pi * silent_present + (1 - pi) * silent_absent
    return np.log(np.maximum(reported, 1e-300)) + m_silent * np.log(np.maximum(silent, 1e-300))


mu_grid = np.linspace(0.0, 2.0, 801)
pi_grid = np.linspace(0.01, 0.99, 197)

print(f"one reported g = {G_OBS} at threshold {CUTOFF}, within-study sd {SIGMA}\n")
print(f"{'m silent':>9s} {'profile argmax':>15s} {'profile range':>26s} "
      f"{'within 2 logL of peak':>24s}")
for m_silent in (1, 2, 5, 10, 20):
    surface = log_likelihood(mu_grid[:, None], pi_grid[None, :], m_silent)
    profile = surface.max(axis=1)
    peak = profile.max()
    best = mu_grid[int(np.argmax(profile))]
    flat = mu_grid[profile >= peak - 2.0]
    print(f"{m_silent:9d} {best:15.3f} {f'{profile.min():.2f} to {peak:.2f}':>26s} "
          f"{f'[{flat.min():.2f}, {flat.max():.2f}]':>24s}")

print("\nFor contrast, the same voxel with two reporting studies (both at g = 0.70):")


def log_likelihood_two(mu, pi, m_silent):
    return log_likelihood(mu, pi, m_silent) + log_likelihood(mu, pi, 0)


for m_silent in (5, 20):
    surface = log_likelihood_two(mu_grid[:, None], pi_grid[None, :], m_silent)
    profile = surface.max(axis=1)
    peak = profile.max()
    flat = mu_grid[profile >= peak - 2.0]
    print(f"  m silent {m_silent:3d}   argmax {mu_grid[int(np.argmax(profile))]:.3f}   "
          f"within 2 logL [{flat.min():.2f}, {flat.max():.2f}]")
