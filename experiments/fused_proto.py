"""Prototype: is a fused numba EM sweep worth building?

The numpy path makes ~25 passes over the per-pair arrays for one EM iteration and every one
of them is memory-bound. A fused kernel keeps each pair in registers and touches memory once.
Benchmarked on the silent half only, at the 0.45M pairs a real fit sees and at 8M.
"""
import math
import time
import numpy as np
from numba import njit
from scipy.special import ndtr

INV_SQRT_2PI = 1.0 / math.sqrt(2.0 * math.pi)
INV_SQRT_2 = 1.0 / math.sqrt(2.0)
PROB_FLOOR = 1e-12
LOGP_FLOOR = 1e-300


@njit(cache=True, inline="always")
def _ndtr_diff(upper, lower):
    """P(lower < Z < upper), accurate in either tail.

    ``0.5 * (erf(u/r2) - erf(l/r2))`` loses every digit when both limits are far out on the
    same side, because erf has saturated at +-1 there; ``erfc`` of the reflected limit has not.
    """
    if upper <= 0.0:
        return 0.5 * (math.erfc(-upper * INV_SQRT_2) - math.erfc(-lower * INV_SQRT_2))
    if lower >= 0.0:
        return 0.5 * (math.erfc(lower * INV_SQRT_2) - math.erfc(upper * INV_SQRT_2))
    return 0.5 * (math.erf(upper * INV_SQRT_2) - math.erf(lower * INV_SQRT_2))


# nogil, and serial: the permutation null already runs whole fits in parallel, so threading
# inside the kernel would oversubscribe. The gain here is fusion, not cores.
@njit(cache=True, nogil=True)
def _silent_sweep(voxel, weight, inv_sigma, inv_sigma_sq, cutoff_scaled, twice_cutoff_scaled,
                  prob_silent_null, mu, pi, width):
    loglik = np.zeros(width)
    claimed = np.zeros(width)
    score = np.zeros(width)
    curvature = np.zeros(width)
    if True:
        for i in range(voxel.size):
            v = voxel[i]
            inv_s = inv_sigma[i]
            upper = cutoff_scaled[i] - mu[v] * inv_s
            lower = upper - twice_cutoff_scaled[i]
            prob = _ndtr_diff(upper, lower)
            if prob < PROB_FLOOR:
                prob = PROB_FLOOR
            pdf_u = INV_SQRT_2PI * math.exp(-0.5 * upper * upper)
            pdf_l = INV_SQRT_2PI * math.exp(-0.5 * lower * lower)
            censor_score = -(pdf_u - pdf_l) * inv_s / prob
            d2 = -(upper * pdf_u - lower * pdf_l) * inv_sigma_sq[i] / prob

            p = pi[v]
            resp = p * prob
            mixture = resp + (1.0 - p) * prob_silent_null[i] + LOGP_FLOOR
            resp = resp / mixture
            w = weight[i]
            loglik[v] += w * math.log(mixture)
            claimed[v] += w * resp
            wr = w * resp
            score[v] += wr * censor_score
            curvature[v] += wr * (d2 - censor_score * censor_score)
    return loglik, claimed, score, curvature


def numpy_sweep(voxel, weight, inv_sigma, inv_sigma_sq, cutoff_scaled, twice_cutoff_scaled,
                prob_silent_null, mu, pi, width):
    upper = mu[voxel] * -inv_sigma
    upper += cutoff_scaled
    lower = upper - twice_cutoff_scaled
    prob = ndtr(upper)
    prob -= ndtr(lower)
    np.clip(prob, PROB_FLOOR, None, out=prob)

    def pdf(x):
        out = x * x
        out *= -0.5
        np.exp(out, out=out)
        out *= INV_SQRT_2PI
        return out

    pdf_u, pdf_l = pdf(upper), pdf(lower)
    censor_score = (pdf_u - pdf_l) * -inv_sigma / prob
    pdf_u *= upper
    pdf_l *= lower
    d2 = (pdf_u - pdf_l) * -inv_sigma_sq / prob

    p = pi[voxel]
    resp = p * prob
    mixture = resp + (1.0 - p) * prob_silent_null + LOGP_FLOOR
    resp = resp / mixture
    loglik = np.bincount(voxel, weights=weight * np.log(mixture), minlength=width)
    claimed = np.bincount(voxel, weights=weight * resp, minlength=width)
    wr = weight * resp
    score = np.bincount(voxel, weights=wr * censor_score, minlength=width)
    curvature = np.bincount(voxel, weights=wr * (d2 - censor_score**2), minlength=width)
    return loglik, claimed, score, curvature


for n, width in ((450_000, 20_000), (8_000_000, 300_000)):
    rng = np.random.default_rng(0)
    voxel = np.sort(rng.integers(0, width, n))
    sigma = np.abs(rng.normal(0.25, 0.05, n)) + 0.05
    inv_sigma = 1.0 / sigma
    cutoff = np.full(n, 0.55)
    args = (voxel, rng.random(n) + 0.5, inv_sigma, inv_sigma**2,
            cutoff * inv_sigma, 2.0 * cutoff * inv_sigma,
            np.clip(rng.random(n), 0.01, 0.99),
            np.abs(rng.normal(0.3, 0.2, width)), np.clip(rng.random(width), 0.01, 0.99), width)
    _silent_sweep(*args)

    def bench(fn, repeats=7):
        ts = []
        for _ in range(repeats):
            t = time.perf_counter()
            out = fn(*args)
            ts.append(time.perf_counter() - t)
        return min(ts), out

    t_np, o_np = bench(numpy_sweep)
    t_nb, o_nb = bench(_silent_sweep)
    print(f"{n/1e6:5.2f}M pairs, {width} voxels:  numpy {t_np*1e3:7.1f} ms   "
          f"numba {t_nb*1e3:7.1f} ms   {t_np/t_nb:5.2f}x")
    for name, a, b in zip(("loglik", "claimed", "score", "curvature"), o_np, o_nb):
        d = np.abs(a - b) / np.maximum(np.abs(a), 1e-12)
        print(f"    {name:10s} max rel {np.max(d):.3e}")
