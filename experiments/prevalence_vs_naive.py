"""Does the fitted `prevalence` beat counting studies with a focus nearby?

`prevalence` is the one output no other method produces: an IBMA gives a magnitude, ALE gives
convergence, and neither claims to say what fraction of studies have an effect at a voxel. If
the feature earns its complexity anywhere it is here. But there is an obvious cheap rival:

    naive(v) = (studies with a focus within r mm of v) / (studies in the collection)

which needs no likelihood, no EM, no selection model and no threshold inference. It is biased
low by construction -- a study that has the effect but failed to clear its threshold is counted
as not having it -- and that is exactly the bias the censored likelihood exists to undo. So the
test is whether undoing it actually helps, measured against a prevalence that is known because
it was simulated.

Both estimators are read from the *same* simulated collections, so nothing differs but the
estimator. Three things are scored, because they can disagree:

  bias         mean estimate minus true prevalence, per site.
  RMSE         total error per site, which is what a user pays.
  ordering     within-map Spearman against the truth and the rate of exactly correct rankings,
               since the ordinal reading is the only one the docstring claims.

If naive wins on RMSE and matches on ordering, the machinery earns nothing for this output and
the honest thing is to say so.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from joblib import Parallel, delayed
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.utils import mm2vox

SITES = [(-32.0, 0.0, 0.0), (-10.0, 0.0, 0.0), (12.0, 0.0, 0.0), (34.0, 0.0, 0.0)]
SITE_PREVALENCE = [0.25, 0.50, 0.75, 1.00]
N_SIMS = int(os.environ.get("NSIMS", 40))
N_STUDIES = 24
U = 3.2905
LOCALISATION_SD = 4.0
N_NOISE = 1
#: Radius for the naive count. 10 mm is the scale of the estimator's own kernel; 15 mm is the
#: looser radius a reader eyeballing a table would use. Both are scored, because the naive
#: estimator's one free parameter should not be allowed to be the reason it loses.
NAIVE_RADII = (10.0, 15.0)


def build_mask():
    shape, step = (25, 21, 21), 4.0
    affine = np.eye(4)
    affine[:3, :3] *= step
    affine[:3, 3] = -step * (np.array(shape) - 1) / 2
    return nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)


def one(seed, effect):
    rng = np.random.default_rng(seed)
    mask = build_mask()
    studies, per_study_points = [], []
    for k in range(N_STUDIES):
        n = int(rng.integers(20, 41))
        meta = {"sample_sizes": [n]}
        points = []
        for site, prevalence in zip(SITES, SITE_PREVALENCE):
            if rng.random() >= prevalence:
                continue
            var = 1.0 / n + effect**2 / (2.0 * n)
            z = rng.normal(effect, np.sqrt(var)) * np.sqrt(n)
            if abs(z) < U:
                continue
            points.append((np.asarray(site) + rng.normal(0, LOCALISATION_SD, 3), z))
        for _ in range(N_NOISE):
            loc = rng.uniform(-44, 44, 3)
            points.append((loc, (U + rng.exponential(1.0 / U)) * rng.choice([-1.0, 1.0])))
        if not points:
            continue
        per_study_points.append(np.array([p[0] for p in points]))
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}", "name": "1", "metadata": meta, "points": [
                {"space": "MNI", "coordinates": [float(c) for c in loc],
                 "values": [{"kind": "Z", "value": float(z)}]} for loc, z in points]}]})
    if len(studies) < 8:
        return None

    # The naive count uses the same denominator the estimator does: every study in the
    # collection, including those whose only focus is noise. A study that dropped out entirely
    # is absent from both, which is the one thing neither estimator can help.
    naive = {}
    for radius in NAIVE_RADII:
        naive[radius] = [
            float(np.mean([np.any(np.linalg.norm(pts - np.asarray(site), axis=1) <= radius)
                           for pts in per_study_points]))
            for site in SITES]

    est = CBES(fwhm=10.0, mask=mask, peak_bias=None, null_method="none")
    result = est.fit(Studyset({"id": "w", "name": "w", "studies": studies},
                              target=None, mask=mask))
    pi_map = result.get_map("prevalence", return_type="array").ravel()
    fitted = []
    for site in SITES:
        at = mm2vox(np.asarray([site]), mask.affine)[0]
        fitted.append(float(pi_map[int(np.ravel_multi_index(tuple(at), mask.shape))]))
    return fitted, naive


def score(name, block, truth):
    bias = block.mean(0) - truth
    rmse = np.sqrt(((block - truth) ** 2).mean(0))
    rho = [stats.spearmanr(row, truth)[0] for row in block]
    exact = np.mean([list(np.argsort(row)) == list(np.argsort(truth)) for row in block])
    print(f"  {name:22s} " + " ".join(f"{v:8.3f}" for v in block.mean(0))
          + f"   mean|bias| {np.abs(bias).mean():.3f}  mean RMSE {rmse.mean():.3f}"
          + f"  rho {np.nanmean(rho):+.3f}  exact {exact:.0%}")


if __name__ == "__main__":
    truth = np.array(SITE_PREVALENCE)
    print(f"{N_STUDIES} coordinate studies, {N_SIMS} replications, U = {U}")
    print("Same collections scored by both estimators; lower RMSE and higher rho are better.\n")
    for effect in (0.8, 0.5, 0.4):
        rows = [r for r in Parallel(n_jobs=8)(delayed(one)(s, effect) for s in range(N_SIMS))
                if r is not None]
        if not rows:
            print(f"effect {effect}: no usable simulation\n")
            continue
        fitted = np.array([r[0] for r in rows], dtype=float)
        print(f"effect {effect}, {len(rows)} usable replications")
        print(f"  {'true prevalence':22s} " + " ".join(f"{v:8.2f}" for v in truth))
        score("CBES prevalence", fitted, truth)
        for radius in NAIVE_RADII:
            block = np.array([r[1][radius] for r in rows], dtype=float)
            score(f"naive count, {radius:.0f} mm", block, truth)
        print(flush=True)
    print("The naive count is biased low by construction -- studies that had the effect but")
    print("missed their threshold are counted as not having it. If it still wins on RMSE, the")
    print("censored likelihood is not buying back more than it costs in variance.")
