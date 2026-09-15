"""Does the window of detectability show up in CBES itself, or only in the idealised model?

An exact occupancy likelihood, fitted on exact detection records with the true detection
function, recovered both prevalence and magnitude well where the detection gradient `dD/dmu` was
large (1.9 to 2.2) and failed at both ends -- nothing to detect below, saturation above. That was
a ceiling on what any estimator could do, not a measurement of this one: CBES must infer or be
told a threshold, tolerate localisation error, work through a smoothing kernel, and fit a
censored mixture rather than a Bernoulli.

Task #49 proposes emitting a per-voxel diagnostic of where a voxel sits relative to that window.
That is only worth building if the window governs CBES's own behaviour, which is what this
measures: sites spanning the detectability range in one fit, with the recoverability of each
factor scored against the detection gradient at its truth.

The threshold is supplied rather than inferred, since threshold inference is separately known to
dominate the prevalence and would confound the comparison.
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

#: Sites spanning the window: undetectable, rising, in the window, saturating, saturated.
SITES = [(-36.0, 0.0, 0.0), (-18.0, 0.0, 0.0), (0.0, 0.0, 0.0),
         (18.0, 0.0, 0.0), (36.0, 0.0, 0.0)]
SITE_MU = [0.15, 0.30, 0.50, 0.75, 1.10]
TRUE_PI = 0.60                     # the same at every site, so only detectability varies
N_SIMS = int(os.environ.get("NSIMS", 40))
N_STUDIES = 24
REPORTING_P = 2.0 * stats.norm.sf(3.2905)
LOCALISATION_SD = 4.0


def build_mask():
    shape, step = (25, 21, 21), 4.0
    affine = np.eye(4)
    affine[:3, :3] *= step
    affine[:3, 3] = -step * (np.array(shape) - 1) / 2
    return nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)


def detection_gradient(mu, n, cut):
    """dD/dmu at the truth, on the t scale the studies actually report."""
    h = 0.02
    def D(m):
        return float(stats.nct.sf(cut, n - 1, m * np.sqrt(n))
                     + stats.nct.cdf(-cut, n - 1, m * np.sqrt(n)))
    return (D(mu + h) - D(mu - h)) / (2 * h)


def one(seed):
    rng = np.random.default_rng(seed)
    mask = build_mask()
    studies, reported = [], np.zeros(len(SITES))
    for k in range(N_STUDIES):
        n = int(rng.integers(20, 41))
        cut = float(stats.t.isf(REPORTING_P / 2.0, n - 1))
        meta = {"sample_sizes": [n], "reporting_threshold": cut}
        points = []
        for s, (site, mu) in enumerate(zip(SITES, SITE_MU)):
            if rng.random() >= TRUE_PI:
                continue
            t = float(stats.nct.rvs(df=n - 1, nc=mu * np.sqrt(n), random_state=rng))
            if abs(t) < cut:
                continue
            reported[s] += 1
            loc = np.asarray(site) + rng.normal(0, LOCALISATION_SD, 3)
            points.append((loc, t))
        loc = rng.uniform(-44, 44, 3)
        points.append((loc, (cut + rng.exponential(1.0 / cut)) * rng.choice([-1.0, 1.0])))
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}", "name": "1", "metadata": meta, "points": [
                {"space": "MNI", "coordinates": [float(c) for c in loc],
                 "values": [{"kind": "T", "value": float(t)}]} for loc, t in points]}]})
    est = CBES(fwhm=10.0, mask=mask, peak_bias=None, null_method="none",
               threshold="reporting_threshold")
    res = est.fit(Studyset({"id": "w", "name": "w", "studies": studies},
                           target=None, mask=mask))
    g = res.get_map("g", return_type="array").ravel()
    pi = res.get_map("prevalence", return_type="array").ravel()
    out = []
    for site in SITES:
        at = mm2vox(np.asarray([site]), mask.affine)[0]
        pos = int(np.ravel_multi_index(tuple(at), mask.shape))
        out.append((float(g[pos]), float(pi[pos])))
    return np.array(out), reported / N_STUDIES


if __name__ == "__main__":
    grads = [detection_gradient(mu, 30, float(stats.t.isf(REPORTING_P / 2.0, 29)))
             for mu in SITE_MU]
    print(f"{N_STUDIES} coordinate-only studies, {N_SIMS} replications, "
          f"true prevalence {TRUE_PI} at every site, threshold supplied")
    print("Only detectability varies across sites, so the prevalence should come back equal")
    print("everywhere if CBES is not affected by where a site sits in the window.\n")
    rows = [r for r in Parallel(n_jobs=6)(delayed(one)(s) for s in range(N_SIMS))
            if r is not None]
    est = np.array([r[0] for r in rows])
    dens = np.array([r[1] for r in rows])
    print(f"{'':22s}" + "".join(f"{f'mu={m:.2f}':>10s}" for m in SITE_MU))
    print(f"{'dD/dmu at the truth':22s}" + "".join(f"{v:10.3f}" for v in grads))
    print(f"{'fraction reporting':22s}" + "".join(f"{v:10.3f}" for v in dens.mean(0)))
    print(f"{'fitted g':22s}" + "".join(f"{v:10.3f}" for v in est[:, :, 0].mean(0)))
    print(f"{'  error vs truth':22s}"
          + "".join(f"{v:+10.3f}" for v in est[:, :, 0].mean(0) - np.array(SITE_MU)))
    print(f"{'fitted prevalence':22s}" + "".join(f"{v:10.3f}" for v in est[:, :, 1].mean(0)))
    print(f"{'  error vs truth':22s}"
          + "".join(f"{v:+10.3f}" for v in est[:, :, 1].mean(0) - TRUE_PI))
    print("\nIf the prevalence error is smallest where dD/dmu is largest, the window governs")
    print("CBES too and the per-voxel diagnostic in task #49 is worth building. If the error")
    print("is flat across sites, the window is a property of the idealised model only.")
