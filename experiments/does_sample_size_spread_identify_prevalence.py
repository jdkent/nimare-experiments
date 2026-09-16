"""Does the coordinate channel need a spread of sample sizes to say anything about prevalence?

`proofs/mixture_identification.py` proves that the reporting indicator's information matrix has
determinant

    det I = pi^2 * sum_{j<k} w_j w_k (u_j v_k - u_k v_j)^2

which is *exactly* zero when every study shares the same (cutoff, sigma). Sigma is set by the
sample size. So a collection whose studies all carry the same sample size and the same threshold
has a coordinate channel that cannot separate mu from pi at all -- it contributes a rank-1
constraint and nothing more. The images keep the model identified, so this is not a claim that the
prevalence becomes meaningless; it is a claim about what the tables add.

The step from the theorem to "therefore requiring a sample size matters" is an inference the
theorem does not make, and an inference of exactly the shape that was wrong an hour earlier
(Schur -> 0 was proved, "therefore the interval is bad" was assumed, and simulation said the
opposite). So it gets measured.

**Kill condition, stated before running.** If the prevalence error at no spread is within 10% of
the error at a 12-to-120 spread, the theorem is still true and the inference from it is wrong: the
tables were never adding separating information that the spread could enable, and requiring a
sample size rests only on the threshold conversion and on the silent dropping of studies that
lack one -- both real, both smaller.

Environment knobs: REPS, VOXELS, NIMG, NTAB, MU, PREV, TAU, MAXITER.
"""
import os
import sys

import numpy as np

sys.path.insert(0, "/home/user/NiMARE")
from nimare.meta.cbma.effectsize import CBES  # noqa: E402

REPS = int(os.environ.get("REPS", "300"))
VOXELS = int(os.environ.get("VOXELS", "32"))
N_IMG = int(os.environ.get("NIMG", "2"))
N_TAB = int(os.environ.get("NTAB", "20"))
MU = float(os.environ.get("MU", "0.5"))
PREV = float(os.environ.get("PREV", "0.6"))
TAU = float(os.environ.get("TAU", "0.15"))
MAX_ITER = int(os.environ.get("MAXITER", "25"))

#: The images are held at a fixed sample size throughout. Only the coordinate studies vary, so
#: any change in the fit is attributable to the channel the theorem is about.
IMAGE_N = 20.0
IMAGE_VAR = 1.0 / IMAGE_N
IMAGE_CUT = 3.09 / np.sqrt(IMAGE_N)


def design(rng, size_range, z_range):
    """Per-study sample sizes and thresholds, log-uniform and uniform respectively."""
    lo, hi = size_range
    sizes = np.full(N_TAB, lo) if lo == hi else np.exp(
        rng.uniform(np.log(lo), np.log(hi), N_TAB)
    )
    lo_z, hi_z = z_range
    zs = np.full(N_TAB, lo_z) if lo_z == hi_z else rng.uniform(lo_z, hi_z, N_TAB)
    return sizes, zs


def cross_products(sizes, zs):
    """The determinant factor the proof identifies, evaluated at the true mu."""
    from scipy.stats import norm

    sd_a, sd_0 = np.sqrt(1.0 / sizes + TAU**2), np.sqrt(1.0 / sizes)
    cut = zs / np.sqrt(sizes)
    hi, lo = (cut - MU) / sd_a, (-cut - MU) / sd_a
    s_a = norm.cdf(hi) - norm.cdf(lo)
    s_0 = norm.cdf(cut / sd_0) - norm.cdf(-cut / sd_0)
    u = -(norm.pdf(hi) - norm.pdf(lo)) / sd_a
    v = s_a - s_0
    return float(sum(
        (u[j] * v[k] - u[k] * v[j]) ** 2
        for j in range(len(u))
        for k in range(j + 1, len(u))
    ))


def one_fit(rng, sizes, zs, width):
    """Generate a voxel set under this design and fit it. Returns (mu_hat, pi_hat)."""
    var_tab, cut_tab = 1.0 / sizes, zs / np.sqrt(sizes)
    sd_a_tab = np.sqrt(var_tab + TAU**2)[:, None]
    sd_0_tab = np.sqrt(var_tab)[:, None]
    sd_a_img, sd_0_img = np.sqrt(IMAGE_VAR + TAU**2), np.sqrt(IMAGE_VAR)

    active_img = rng.random((N_IMG, width)) < PREV
    g_img = np.where(
        active_img,
        MU + sd_a_img * rng.standard_normal((N_IMG, width)),
        sd_0_img * rng.standard_normal((N_IMG, width)),
    )
    active_tab = rng.random((N_TAB, width)) < PREV
    latent = np.where(
        active_tab,
        MU + sd_a_tab * rng.standard_normal((N_TAB, width)),
        sd_0_tab * rng.standard_normal((N_TAB, width)),
    )
    indicator = np.where(np.abs(latent) >= cut_tab[:, None], -1, 1)

    n = N_IMG + N_TAB
    weights = np.zeros((n, width))
    weights[:N_IMG] = 1.0
    values = np.zeros((n, width))
    values[:N_IMG] = g_img
    variances = np.full((n, width), IMAGE_VAR)
    full_indicator = np.zeros((n, width))
    full_indicator[N_IMG:] = indicator

    est = CBES(selection_model="zero-inflated", interval="wald", tau2_method="none",
               max_iter=MAX_ITER)
    mu, pi, *_ = est._fit_chunk(
        weights=weights,
        g_obs=values,
        var_obs=variances,
        indicator=full_indicator,
        tau2=np.full(width, TAU**2),
        null_var=np.concatenate([np.full(N_IMG, IMAGE_VAR), var_tab]),
        cutoffs=np.concatenate([np.full(N_IMG, IMAGE_CUT), cut_tab]),
        start=g_img.mean(axis=0),
    )
    return mu, pi


ARMS = {
    "no spread at all": ((20, 20), (3.09, 3.09)),
    "sample sizes 12-40": ((12, 40), (3.09, 3.09)),
    "sample sizes 12-120": ((12, 120), (3.09, 3.09)),
    "thresholds 2.3-4.5": ((20, 20), (2.3, 4.5)),
    "both spread": ((12, 120), (2.3, 4.5)),
}


def main():
    print(f"true mu {MU}, prevalence {PREV}, tau {TAU}, {N_IMG} images at n={IMAGE_N:.0f}, "
          f"{N_TAB} tables, max_iter {MAX_ITER}")
    print(f"{REPS * VOXELS} fits per arm")
    print()
    print(f"{'arm':22s} {'cross-prod':>11s} {'pi bias':>9s} {'pi rmse':>9s} "
          f"{'mu bias':>9s} {'mu rmse':>9s}")
    results = {}
    for label, (size_range, z_range) in ARMS.items():
        # Two streams, so the arms are paired. Drawing the design from the same generator as the
        # data would desynchronise them -- the no-spread arm consumes no random numbers while
        # choosing its design and every other arm consumes 20 or 40 -- and the comparison would
        # then be partly between different noise realisations.
        design_rng = np.random.default_rng(1)
        data_rng = np.random.default_rng(0)
        pis, mus, dets = [], [], []
        for _ in range(REPS):
            sizes, zs = design(design_rng, size_range, z_range)
            dets.append(cross_products(sizes, zs))
            mu, pi = one_fit(data_rng, sizes, zs, VOXELS)
            mus.append(mu)
            pis.append(pi)
        pi_err = np.concatenate(pis) - PREV
        mu_err = np.concatenate(mus) - MU
        results[label] = np.sqrt((pi_err**2).mean())
        print(f"{label:22s} {np.mean(dets):11.4f} {pi_err.mean():+9.4f} "
              f"{np.sqrt((pi_err**2).mean()):9.4f} {mu_err.mean():+9.4f} "
              f"{np.sqrt((mu_err**2).mean()):9.4f}")

    flat = results["no spread at all"]
    wide = results["sample sizes 12-120"]
    print()
    print(f"prevalence rmse, no spread {flat:.4f} against widest spread {wide:.4f}: "
          f"{100 * (flat - wide) / flat:+.1f}%")
    if flat - wide < 0.10 * flat:
        print("KILL CONDITION MET: the spread buys under 10%. The theorem stands and the")
        print("inference from it does not -- requiring a sample size rests on the threshold")
        print("conversion and the silent drop, not on identification.")
    else:
        print("The spread buys more than 10%, so sample-size heterogeneity is doing the")
        print("identifying work the theorem says it must, and a defaulted sample size removes it.")


if __name__ == "__main__":
    main()
