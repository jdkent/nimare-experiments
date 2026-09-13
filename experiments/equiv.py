"""Check the sparse-pair EM against a literal transcription of the dense version."""
import numpy as np
from scipy.special import ndtr
from nimare.meta.cbma import effectsize as E
from nimare.meta.cbma import CBES

def dense_fit_chunk(self, *, weights, g_obs, var_obs, covered, tau2, null_var, cutoffs, start):
    zero_inflated = self.selection_model == "zero-inflated"
    reports = weights > 0
    silent = ~covered
    sigma_reported = np.sqrt(var_obs + tau2)
    sigma_silent = np.sqrt(null_var + tau2)
    sigma_null = np.sqrt(np.broadcast_to(null_var, silent.shape))
    weight_reported = np.where(reports, weights, 0.0)
    n_reporting = reports.sum(axis=0)
    reporter_scale = np.divide(weight_reported.sum(axis=0), n_reporting,
                               out=np.ones(weights.shape[1]), where=n_reporting > 0)
    weight_silent = silent.astype(float) * reporter_scale
    prob_silent_null = np.clip(ndtr(cutoffs / sigma_null) - ndtr(-cutoffs / sigma_null), 1e-12, None)
    density_null = E._normal_pdf(g_obs / np.sqrt(var_obs)) / np.sqrt(var_obs)

    def derivs(mu, wr, ws):
        precision = 1.0 / sigma_reported**2
        score = (wr * (g_obs - mu) * precision).sum(axis=0)
        curvature = -(wr * precision).sum(axis=0)
        upper = (cutoffs - mu) / sigma_silent
        lower = (-cutoffs - mu) / sigma_silent
        prob = np.clip(ndtr(upper) - ndtr(lower), 1e-12, None)
        pu, pl = E._normal_pdf(upper), E._normal_pdf(lower)
        d_prob = -(pu - pl) / sigma_silent
        d2_prob = -(upper * pu - lower * pl) / sigma_silent**2
        cs = d_prob / prob
        return score + (ws * cs).sum(axis=0), curvature + (ws * (d2_prob / prob - cs**2)).sum(axis=0)

    mu = start.copy()
    pi = np.full(mu.shape, 0.5 if zero_inflated else 1.0)
    total_weight = weight_reported.sum(axis=0) + weight_silent.sum(axis=0)
    resp_rep = np.ones_like(g_obs); resp_sil = np.ones_like(g_obs)
    for _ in range(self.max_iter):
        if zero_inflated:
            de = E._normal_pdf((g_obs - mu) / sigma_reported) / sigma_reported
            resp_rep = pi * de
            resp_rep = resp_rep / (resp_rep + (1.0 - pi) * density_null + 1e-300)
            upper = (cutoffs - mu) / sigma_silent; lower = (-cutoffs - mu) / sigma_silent
            pse = np.clip(ndtr(upper) - ndtr(lower), 1e-12, None)
            resp_sil = pi * pse
            resp_sil = resp_sil / (resp_sil + (1.0 - pi) * prob_silent_null + 1e-300)
            num = (weight_reported * resp_rep).sum(axis=0) + (weight_silent * resp_sil).sum(axis=0)
            pi = np.clip(np.divide(num, total_weight, out=np.zeros_like(mu), where=total_weight > 0),
                         1e-4, 1.0 - 1e-4)
        score, curv = derivs(mu, weight_reported * resp_rep, weight_silent * resp_sil)
        mu = mu + np.clip(np.where(curv < 0, -score / curv, 0.0), -1.0, 1.0)
    _, curv = derivs(mu, weight_reported * resp_rep, weight_silent * resp_sil)
    se = np.full(mu.shape, np.inf); ok = curv < 0
    se[ok] = 1.0 / np.sqrt(-curv[ok])
    return mu, pi, se

rng = np.random.default_rng(0)
E._EM_COMPACTION_FRACTION = 2.0      # never compact
E._EM_TOLERANCE = 0.0                # never retire early -> run all max_iter, like the dense loop

for model in ("zero-inflated", "tobit"):
    for trial in range(5):
        K, V = rng.integers(5, 25), rng.integers(50, 400)
        weights = np.where(rng.random((K, V)) < 0.4, rng.random((K, V)), 0.0)
        g_obs = np.where(weights > 0, rng.normal(0.5, 0.4, (K, V)), 0.0)
        var_obs = rng.uniform(0.02, 0.1, (K, V))
        covered = rng.random((K, V)) < 0.6
        tau2 = rng.uniform(0, 0.05, V)
        null_var = rng.uniform(0.02, 0.08, (K, 1))
        cutoffs = rng.uniform(0.3, 0.8, (K, 1))
        start = rng.normal(0.5, 0.2, V)
        est = CBES(selection_model=model, max_iter=12, null_method="parametric")
        kw = dict(weights=weights, g_obs=g_obs, var_obs=var_obs, covered=covered,
                  tau2=tau2, null_var=null_var, cutoffs=cutoffs, start=start)
        a = est._fit_chunk(**kw)
        b = dense_fit_chunk(est, **kw)
        for name, x, y in zip(("mu", "pi", "se"), a, b):
            finite = np.isfinite(x) & np.isfinite(y)
            assert np.array_equal(np.isfinite(x), np.isfinite(y)), f"{model} {name} inf mismatch"
            assert np.allclose(x[finite], y[finite], rtol=1e-10, atol=1e-12), (
                f"{model} {name} max diff {np.abs(x[finite]-y[finite]).max():.3e}")
    print(f"{model}: sparse EM matches the dense reference on 5 random problems")
