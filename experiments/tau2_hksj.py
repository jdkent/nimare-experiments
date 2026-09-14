"""Test 4: is kernel-weighted DerSimonian-Laird good enough, and does HKSJ help the SE?

DL is a one-step moment estimator and is known to be downward-biased with few studies and
heterogeneous within-study variances -- which is CBES's situation, since the effective number of
studies per voxel is a handful. Paule-Mandel and REML are the usual recommendations. And the
pooled SE, sqrt(sum w^2/V)/sum(w/V), treats tau2-hat as known, which is exactly the case
Hartung-Knapp-Sidik-Jonkman exists to fix.

Three questions, in order of what would change the code:

  1. Do PM and REML recover a known tau2 better than DL at CBES's study counts?
  2. Does HKSJ give better coverage of the true pooled effect than the model-based SE?
  3. Do either matter on a real fit, or is the difference lost in the other errors?

The kernel-weighted forms are checked against PyMARE with every weight set to 1 first -- if they
do not reduce to the textbook estimators there, nothing downstream is worth reading.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from scipy import stats
from scipy.optimize import brentq
from pymare import Dataset as PyMAREDataset
from pymare.estimators import DerSimonianLaird as PyMARE_DL
from pymare.estimators import VarianceBasedLikelihoodEstimator

from nimare.meta.cbma.effectsize import _local_dersimonian_laird

MAX_TAU2 = 100.0


def dl(g, s2, w):
    """The estimator CBES uses, called through its own moment sums."""
    a = w / s2
    return float(_local_dersimonian_laird(
        np.array([w.sum()]), np.array([a.sum()]), np.array([(a**2).sum()]),
        np.array([(a * g).sum()]), np.array([(a * g * g).sum()]),
        np.array([(w**2 / s2).sum()]), np.array([len(g)]))[0])


def _weighted_q_minus_expectation(tau2, g, s2, w):
    """Kernel-weighted generalized Q at tau2, minus its own expectation at that tau2.

    Paule-Mandel solves Q(tau2) = k - 1. The kernel-weighted expectation is not ``sum(w) - 1``:
    taking the DL derivation's E[Q] with ``s2 + tau2`` in place of ``s2`` gives

        E[Q] = sum(w) - sum(a^2 (s2 + tau2)) / sum(a),    a = w / (s2 + tau2)

    which collapses to ``k - 1`` when every w is 1, and stays positive when the weights are
    small. ``sum(w) - 1`` does not: at a voxel reached only by distant foci the weights sum to
    less than one, the target goes negative, Q can never fall to it and the solver returns the
    bracket's upper end. That drove the mean estimate to 2.1 against a true 0.0 at k = 3.
    """
    a = w / (s2 + tau2)
    mean = (a * g).sum() / a.sum()
    q = (a * (g - mean) ** 2).sum()
    expected = w.sum() - (a**2 * (s2 + tau2)).sum() / a.sum()
    return q - expected


def paule_mandel(g, s2, w):
    if len(g) < 2:
        return 0.0
    if _weighted_q_minus_expectation(0.0, g, s2, w) <= 0:
        return 0.0
    if _weighted_q_minus_expectation(MAX_TAU2, g, s2, w) > 0:
        return MAX_TAU2
    return float(
        brentq(_weighted_q_minus_expectation, 0.0, MAX_TAU2, args=(g, s2, w), xtol=1e-10)
    )


def _neg_reml(tau2, g, s2, w):
    a = w / (s2 + tau2)
    mean = (a * g).sum() / a.sum()
    # Weighted REML log-likelihood, up to a constant; w enters as a case weight.
    return 0.5 * (
        (w * np.log(s2 + tau2)).sum()
        + (a * (g - mean) ** 2).sum()
        + np.log(a.sum())
    )


def reml(g, s2, w, rescale=True):
    """Weighted REML. ``rescale`` normalizes the weights to sum to k, and is not cosmetic.

    As tau2 grows the objective tends to ``0.5 (sum(w) - 1) log tau2``, so with weights summing
    to less than one it is unbounded below and the estimate runs to whatever the grid allows --
    2.1 against a true 0.0 at k = 3 before this. Rescaling to sum to k keeps the relative
    weighting, restores properness, and leaves the every-w-is-1 case untouched.
    """
    if len(g) < 2:
        return 0.0
    if rescale and w.sum() > 0:
        w = w * (len(g) / w.sum())
    grid = np.concatenate([[0.0], np.geomspace(1e-6, MAX_TAU2, 400)])
    values = [_neg_reml(t, g, s2, w) for t in grid]
    return float(grid[int(np.argmin(values))])


def pooled(g, s2, w, tau2):
    a = w / (s2 + tau2)
    est = (a * g).sum() / a.sum()
    se_model = np.sqrt((w**2 / (s2 + tau2)).sum()) / a.sum()
    # HKSJ: replace the model SE by a weighted residual variance, on k_eff - 1 df.
    #
    # k_eff is Kish's (sum w)^2 / sum w^2, not sum(w). Both give k when every weight is 1, but
    # sum(w) is not scale-invariant in the weights: at a voxel reached only by distant foci it
    # falls below 1, the df goes to zero and the t critical value explodes -- which read as
    # HKSJ covering 99.3% everywhere before this, an artefact of the df and not of HKSJ. Kish's
    # form is invariant to rescaling the weights, and is what the estimator already reports as
    # its "n_eff" map.
    k_eff = (w.sum() ** 2) / (w**2).sum() if (w**2).sum() > 0 else 0.0
    if len(g) > 1 and k_eff > 1:
        se_hksj = np.sqrt((a * (g - est) ** 2).sum() / ((k_eff - 1.0) * a.sum()))
    else:
        se_hksj = se_model
    return est, se_model, max(se_hksj, 1e-12), k_eff


print("=" * 74)
print("oracle: with every weight 1, do these reduce to the textbook estimators?")
print("=" * 74)
rng = np.random.default_rng(0)
for k in (3, 8, 20):
    g = rng.normal(0.4, 0.3, k)
    s2 = np.abs(rng.normal(0.08, 0.02, k)) + 0.01
    w = np.ones(k)
    ds = PyMAREDataset(y=g[:, None], v=s2[:, None])
    ref_dl = float(np.asarray(PyMARE_DL().fit_dataset(ds).params_["tau2"]).ravel()[0])
    ref_reml = float(np.asarray(
        VarianceBasedLikelihoodEstimator(method="reml").fit_dataset(ds).params_["tau2"]
    ).ravel()[0])
    print(f"  k={k:3d}  DL ours {dl(g, s2, w):.6f} vs pymare {ref_dl:.6f}   "
          f"REML ours {reml(g, s2, w):.4f} vs pymare {ref_reml:.4f}   "
          f"PM ours {paule_mandel(g, s2, w):.4f}")

print()
print("=" * 74)
print("1. recovery of a known tau2, 4000 replicates per cell")
print("=" * 74)
print(f"{'k':>4s} {'true tau2':>10s} {'DL':>18s} {'Paule-Mandel':>18s} {'REML':>18s}")
for k in (3, 5, 10):
    for true_tau2 in (0.0, 0.05, 0.20):
        rng = np.random.default_rng(1)
        got = {"dl": [], "pm": [], "reml": []}
        for _ in range(4000):
            s2 = np.abs(rng.normal(0.08, 0.03, k)) + 0.01
            w = rng.uniform(0.2, 1.0, k)  # kernel weights, as a real voxel sees them
            g = rng.normal(0.4, np.sqrt(s2 + true_tau2))
            got["dl"].append(dl(g, s2, w))
            got["pm"].append(paule_mandel(g, s2, w))
            got["reml"].append(reml(g, s2, w))
        cells = []
        for name in ("dl", "pm", "reml"):
            v = np.array(got[name])
            cells.append(f"{v.mean():.3f} ({v.mean() - true_tau2:+.3f})")
        print(f"{k:4d} {true_tau2:10.2f} {cells[0]:>18s} {cells[1]:>18s} {cells[2]:>18s}",
              flush=True)

print()
print("=" * 74)
print("2. coverage of the true pooled effect by a nominal 95% interval")
print("=" * 74)
print(f"{'k':>4s} {'true tau2':>10s} {'model SE (z)':>14s} {'HKSJ (t)':>10s} {'model SE (t)':>14s}")
for k in (3, 5, 10):
    for true_tau2 in (0.0, 0.05, 0.20):
        rng = np.random.default_rng(2)
        hit = {"model_z": 0, "hksj_t": 0, "model_t": 0}
        trials = 4000
        for _ in range(trials):
            s2 = np.abs(rng.normal(0.08, 0.03, k)) + 0.01
            w = rng.uniform(0.2, 1.0, k)
            g = rng.normal(0.4, np.sqrt(s2 + true_tau2))
            tau2 = paule_mandel(g, s2, w)
            est, se_model, se_hksj, k_eff = pooled(g, s2, w, tau2)
            t_crit = stats.t.ppf(0.975, max(k_eff - 1.0, 1.0))
            hit["model_z"] += abs(est - 0.4) <= 1.959964 * se_model
            hit["hksj_t"] += abs(est - 0.4) <= t_crit * se_hksj
            hit["model_t"] += abs(est - 0.4) <= t_crit * se_model
        print(f"{k:4d} {true_tau2:10.2f} {100 * hit['model_z'] / trials:13.1f}% "
              f"{100 * hit['hksj_t'] / trials:9.1f}% {100 * hit['model_t'] / trials:13.1f}%",
              flush=True)
