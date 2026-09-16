"""Does the boundary test for pi = 1 have its nominal size, and where does its power come from?

`proofs/boundary_test_for_full_prevalence.py` derives the test: pi = 1 sits on the edge of the
parameter space, so 2 log Lambda follows Chernoff's half-and-half mixture of chi^2_0 and chi^2_1
rather than chi^2_1, and the level-0.05 critical value is 2.7055 rather than 3.8415. The proof
also predicts that power rises with the number of coordinate tables, through I_pipi = n_t v^2 / D,
even though the same tables contribute nothing to separating mu from pi.

Both are claims about a specific likelihood, and the boundary result is asserted in the proof
rather than derived, so neither is established until measured.

**Kill conditions, stated before running.**

  1. If the empirical size at 2.7055 is outside 0.05 +- 3 binomial standard errors, the boundary
     distribution does not describe this likelihood and the critical value cannot be used as
     derived. Nothing downstream survives that.
  2. If the size at 3.8415 is not close to half the size at 2.7055, the factor-of-two claim --
     the whole practical content of the boundary correction -- is wrong.
  3. If power does not rise with the number of tables, the proof's I_pipi argument is wrong and
     the test cannot be sold as something a coordinate corpus helps with.

The statistic is computed on a grid rather than by the estimator's EM, deliberately: the question
here is whether the *test* is calibrated, not whether the optimiser finds the maximum. Confounding
those two is how the max_iter result stayed hidden for so long.

Environment knobs: REPS, NIMG, NTAB, MU, TAU, PREVS, NTABS, NSUBJ.
"""
import os
import sys

import numpy as np
from scipy.special import ndtr

REPS = int(os.environ.get("REPS", "4000"))
N_IMG = int(os.environ.get("NIMG", "2"))
N_TAB = int(os.environ.get("NTAB", "20"))
MU = float(os.environ.get("MU", "0.5"))
TAU = float(os.environ.get("TAU", "0.15"))
PREVS = [float(x) for x in os.environ.get("PREVS", "1.0,0.9,0.8,0.6,0.4").split(",")]
NTABS = [int(x) for x in os.environ.get("NTABS", "5,20,80").split(",")]

#: Sample size of the coordinate studies. The per-table information about pi is v^2 / D with
#: v = s_a - s_0, and s_0 is fixed at 2*Phi(z) - 1 whatever the sample size, because c / sigma_0
#: is just z. So the whole of v is carried by s_a, which falls steeply as studies grow: a large
#: study with the effect reports it, a small one stays silent either way and its silence says
#: nothing about which component it came from. I_pipi runs 0.124 to 21.2 from n=12 to n=200.
N_SUBJ = float(os.environ.get("NSUBJ", "20"))
VAR_G = 1.0 / N_SUBJ
CUT = 3.09 / np.sqrt(N_SUBJ)
SD_A = np.sqrt(VAR_G + TAU**2)
SD_0 = np.sqrt(VAR_G)

BOUNDARY_CUT = 2.7055      # P(chi^2_1 > t) = 0.10, the half-mixture's 0.05 point
NAIVE_CUT = 3.8415         # P(chi^2_1 > t) = 0.05, what the naive recipe would use

#: Dense near the boundary, because that is where the constrained and unconstrained maxima part
#: company and a coarse grid there would manufacture zeros and a conservative test.
PI_GRID = np.unique(np.concatenate([[1.0], 1.0 - np.logspace(-3.0, 0.0, 160)[:-1]]))
MU_GRID = np.linspace(-0.6, 2.0, 521)

SILENT_NULL = float(ndtr(CUT / SD_0) - ndtr(-CUT / SD_0))


def draw(rng, reps, prevalence, n_tab):
    """One replicate per column: image values and reporting indicators."""
    active_img = rng.random((N_IMG, reps)) < prevalence
    g_img = np.where(
        active_img,
        MU + SD_A * rng.standard_normal((N_IMG, reps)),
        SD_0 * rng.standard_normal((N_IMG, reps)),
    )
    active_tab = rng.random((n_tab, reps)) < prevalence
    latent = np.where(
        active_tab,
        MU + SD_A * rng.standard_normal((n_tab, reps)),
        SD_0 * rng.standard_normal((n_tab, reps)),
    )
    reported = np.abs(latent) >= CUT
    return g_img, reported.sum(axis=0), n_tab - reported.sum(axis=0)


def log_likelihood(g_img, n_rep, n_sil):
    """``ll[i_mu, i_pi, replicate]`` for the zero-inflated mixture."""
    mu = MU_GRID[:, None, None]
    pi = PI_GRID[None, :, None]

    silent_active = np.clip(
        ndtr((CUT - mu) / SD_A) - ndtr((-CUT - mu) / SD_A), 1e-300, None
    )
    mix_silent = pi * silent_active + (1.0 - pi) * SILENT_NULL
    total = n_sil[None, None, :] * np.log(np.clip(mix_silent, 1e-300, None))
    total = total + n_rep[None, None, :] * np.log(np.clip(1.0 - mix_silent, 1e-300, None))

    for study in range(N_IMG):
        value = g_img[study][None, None, :]
        density_active = np.exp(-0.5 * ((value - mu) / SD_A) ** 2) / SD_A
        density_null = np.exp(-0.5 * (value / SD_0) ** 2) / SD_0
        total = total + np.log(
            np.clip(pi * density_active + (1.0 - pi) * density_null, 1e-300, None)
        )
    return total


def statistic(g_img, n_rep, n_sil, chunk=250):
    """2 log Lambda for H0: pi = 1 against the one-sided alternative pi < 1."""
    out = np.empty(n_rep.size)
    boundary = int(np.flatnonzero(PI_GRID == 1.0)[0])
    for lo in range(0, n_rep.size, chunk):
        hi = min(lo + chunk, n_rep.size)
        ll = log_likelihood(g_img[:, lo:hi], n_rep[lo:hi], n_sil[lo:hi])
        unconstrained = ll.max(axis=(0, 1))
        constrained = ll[:, boundary, :].max(axis=0)
        out[lo:hi] = 2.0 * (unconstrained - constrained)
    # A maximum on the boundary gives exactly zero, which is the chi^2_0 atom.
    return np.maximum(out, 0.0)


def main():
    rng = np.random.default_rng(0)
    se = np.sqrt(0.05 * 0.95 / REPS)
    print(f"true mu {MU}, tau {TAU}, {N_IMG} images, cut z=3.09 at n={N_SUBJ}")
    print(f"{REPS} replicates per cell; binomial se at 0.05 is {se:.4f}")
    print(f"boundary critical value {BOUNDARY_CUT}, naive chi^2_1 value {NAIVE_CUT}")
    print()
    print(f"{'tables':>7} {'prevalence':>11} {'atom at 0':>10} {'reject@2.71':>12} "
          f"{'reject@3.84':>12}")

    size_at_boundary = {}
    power = {}
    for n_tab in NTABS:
        for prevalence in PREVS:
            g_img, n_rep, n_sil = draw(rng, REPS, prevalence, n_tab)
            stat = statistic(g_img, n_rep, n_sil)
            atom = float(np.mean(stat <= 1e-9))
            r_boundary = float(np.mean(stat > BOUNDARY_CUT))
            r_naive = float(np.mean(stat > NAIVE_CUT))
            if prevalence == 1.0:
                size_at_boundary[n_tab] = (r_boundary, r_naive)
            else:
                power[(n_tab, prevalence)] = r_boundary
            print(f"{n_tab:7d} {prevalence:11.2f} {atom:10.3f} {r_boundary:12.4f} "
                  f"{r_naive:12.4f}")

    print()
    verdicts = []
    for n_tab, (r_boundary, r_naive) in sorted(size_at_boundary.items()):
        ok = abs(r_boundary - 0.05) <= 3 * se
        verdicts.append(ok)
        print(f"  size at {n_tab} tables: {r_boundary:.4f} "
              f"{'within' if ok else 'OUTSIDE'} 0.05 +- {3 * se:.4f}")
        ratio = r_naive / r_boundary if r_boundary > 0 else float("nan")
        print(f"    naive cut rejects {r_naive:.4f}, a ratio of {ratio:.2f} "
              f"(the factor-of-two claim wants about 0.5)")

    print()
    for prevalence in PREVS:
        if prevalence == 1.0:
            continue
        row = [f"{power.get((n, prevalence), float('nan')):.3f}" for n in NTABS]
        print(f"  power at prevalence {prevalence:.1f}: " +
              "  ".join(f"{n} tables {p}" for n, p in zip(NTABS, row)))

    rising = all(
        power.get((NTABS[0], p), 0) <= power.get((NTABS[-1], p), 0)
        for p in PREVS if p != 1.0
    )
    print()
    print(f"  power rises with the number of tables: {rising}")
    if not all(verdicts):
        print("  KILL CONDITION 1 MET: the boundary distribution does not describe this")
        print("  likelihood, and the derived critical value cannot be used.")


if __name__ == "__main__":
    main()
