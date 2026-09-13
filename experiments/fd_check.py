"""Finite-difference check of CBES's analytic score and curvature in mu."""
import numpy as np
from scipy.special import ndtr
from nimare.meta.cbma.effectsize import _censoring_terms, _mu_derivatives

def loglik(mu, *, g_rep, sigma_rep, w_rep, rep_voxel, cutoff_sil, sigma_sil, w_sil,
           sil_voxel, width):
    """The mu-dependent part of the weighted log-likelihood the derivatives come from."""
    mu_rep = mu[rep_voxel]
    dens = -0.5 * ((g_rep - mu_rep) / sigma_rep) ** 2
    upper = (cutoff_sil - mu[sil_voxel]) / sigma_sil
    lower = (-cutoff_sil - mu[sil_voxel]) / sigma_sil
    cens = np.log(np.clip(ndtr(upper) - ndtr(lower), 1e-300, None))
    return (np.bincount(rep_voxel, weights=w_rep * dens, minlength=width)
            + np.bincount(sil_voxel, weights=w_sil * cens, minlength=width))

rng = np.random.default_rng(0)
worst_score = worst_curv = 0.0
for trial in range(6):
    width = 40
    n_rep, n_sil = 120, 150
    rep_voxel = rng.integers(0, width, n_rep)
    sil_voxel = rng.integers(0, width, n_sil)
    g_rep = rng.normal(0.5, 0.4, n_rep)
    sigma_rep = rng.uniform(0.15, 0.45, n_rep)
    w_rep = rng.uniform(0.1, 1.0, n_rep)
    cutoff_sil = rng.uniform(0.3, 0.9, n_sil)
    sigma_sil = rng.uniform(0.15, 0.45, n_sil)
    w_sil = rng.uniform(0.1, 1.0, n_sil)
    mu = rng.normal(0.4, 0.3, width)

    kw = dict(g_rep=g_rep, sigma_rep=sigma_rep, w_rep=w_rep, rep_voxel=rep_voxel,
              cutoff_sil=cutoff_sil, sigma_sil=sigma_sil, w_sil=w_sil,
              sil_voxel=sil_voxel, width=width)

    censoring = _censoring_terms(mu[sil_voxel], cutoff_sil, sigma_sil)
    score, curvature = _mu_derivatives(
        width=width, mu_rep=mu[rep_voxel], g_rep=g_rep,
        precision_rep=1.0 / sigma_rep**2, rep_voxel=rep_voxel, weight_rep=w_rep,
        sil_voxel=sil_voxel, weight_sil=w_sil, censoring=censoring)

    h = 1e-5
    fd_score = np.zeros(width); fd_curv = np.zeros(width)
    for v in range(width):
        up, dn = mu.copy(), mu.copy()
        up[v] += h; dn[v] -= h
        lu, ld, l0 = loglik(up, **kw)[v], loglik(dn, **kw)[v], loglik(mu, **kw)[v]
        fd_score[v] = (lu - ld) / (2 * h)
        fd_curv[v] = (lu - 2 * l0 + ld) / h**2

    worst_score = max(worst_score, np.max(np.abs(score - fd_score) / (np.abs(fd_score) + 1e-6)))
    worst_curv = max(worst_curv, np.max(np.abs(curvature - fd_curv) / (np.abs(fd_curv) + 1e-6)))

print(f"max relative error vs finite differences over 6 random problems:")
print(f"  score     {worst_score:.2e}")
print(f"  curvature {worst_curv:.2e}")
print("  -> analytic derivatives agree" if max(worst_score, worst_curv) < 1e-4
      else "  -> MISMATCH: analytic derivatives are wrong")
