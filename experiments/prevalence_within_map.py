"""Does prevalence order *voxels within one map*, which is what the docstring claims?

`prevalence_calibration` varied prevalence across collections at a single voxel and found a floor
near 0.2 at a true zero, affine but unstable coefficients, and non-monotonicity at weak effects.
None of that bears on the claim the estimator actually makes: that "ordering survives, so
comparing voxels within one map is sound" -- an ordering across *voxels inside one fit*, not
across collections.

So one collection is built containing four well-separated sites whose prevalences differ by
construction, and the question is whether the fitted map ranks those four correctly.

The studyset is assembled directly rather than through the simulator, because the simulator
controls one prevalence for the whole collection and `Studyset.coordinates` is read-only. Each
study independently has each site's effect with that site's probability, draws its own observed
effect size, and reports a focus only if that clears the reporting threshold -- so the reporting
fraction per site is recorded alongside, and is lower than the prevalence by construction.
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
N_SIMS = 16
N_STUDIES = 24
#: Reporting threshold as a two-sided p, so each study's cut sits on its own t scale. A fixed
#: number on the z scale would hand the estimator a statistic whose convention it does not
#: assume -- it reads a reported value as a t on n - 1 degrees of freedom -- and that inflates
#: every recovered magnitude by about a third. See PROTOCOL.md.
REPORTING_P = 2.0 * stats.norm.sf(3.2905)
LOCALISATION_SD = 4.0      # mm of jitter between the true site and the reported focus
N_NOISE = 1


def build_mask():
    shape, step = (25, 21, 21), 4.0
    affine = np.eye(4)
    affine[:3, :3] *= step
    affine[:3, 3] = -step * (np.array(shape) - 1) / 2
    return nib.Nifti1Image(np.ones(shape, dtype=np.int32), affine)


def one(seed, effect, threshold="study-min"):
    rng = np.random.default_rng(seed)
    mask = build_mask()
    studies, reported = [], np.zeros(len(SITES))
    for k in range(N_STUDIES):
        n = int(rng.integers(20, 41))
        cut = float(stats.t.isf(REPORTING_P / 2.0, n - 1))
        meta = {"sample_sizes": [n]}
        points = []
        for s, (site, prevalence) in enumerate(zip(SITES, SITE_PREVALENCE)):
            if rng.random() >= prevalence:
                continue                      # this study has no effect at this site
            # A genuine noncentral t: the effect is the noncentrality and the denominator
            # carries its own n - 1 degrees of freedom, which is the statistic the estimator
            # assumes a reported value to be.
            z = float(stats.nct.rvs(df=n - 1, nc=effect * np.sqrt(n), random_state=rng))
            if abs(z) < cut:
                continue                      # had the effect, failed to clear the threshold
            reported[s] += 1
            loc = np.asarray(site) + rng.normal(0, LOCALISATION_SD, 3)
            points.append((loc, z))
        for _ in range(N_NOISE):
            loc = rng.uniform(-44, 44, 3)
            # A noise peak just clears the cut, with the exponential overshoot of a smooth field.
            points.append((loc, (cut + rng.exponential(1.0 / cut)) * rng.choice([-1.0, 1.0])))
        if len(points) < 1:
            continue
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}", "name": "1", "metadata": meta, "points": [
                {"space": "MNI", "coordinates": [float(c) for c in loc],
                 "values": [{"kind": "T", "value": float(z)}]} for loc, z in points]}]})
    if len(studies) < 8:
        return None
    est = CBES(fwhm=10.0, mask=mask, peak_bias=None, null_method="none",
               threshold=threshold)
    result = est.fit(Studyset({"id": "w", "name": "w", "studies": studies},
                              target=None, mask=mask))
    pi_map = result.get_map("prevalence", return_type="array").ravel()
    out = []
    for site in SITES:
        at = mm2vox(np.asarray([site]), mask.affine)[0]
        out.append(float(pi_map[int(np.ravel_multi_index(tuple(at), mask.shape))]))
    return out + list(reported / max(len(studies), 1))


for effect, threshold, label in ((0.8, "study-min", "study-min"),
                                 (0.8, 3.2905267314919255, "fixed 3.2905"),
                                 (0.4, "study-min", "study-min"),
                                 (0.4, 3.2905267314919255, "fixed 3.2905")):
    rows = [r for r in Parallel(n_jobs=6)(
        delayed(one)(seed, effect, threshold) for seed in range(N_SIMS)) if r is not None]
    if not rows:
        print(f"effect {effect}: no usable simulation\n", flush=True)
        continue
    block = np.array(rows, dtype=float)
    fitted, reporting = block[:, :len(SITES)], block[:, len(SITES):]
    print(f"effect {effect}, threshold={label}, {len(rows)} sims, {N_STUDIES} studies")
    print(f"  {'site prevalence':>18} " + " ".join(f"{p:8.2f}" for p in SITE_PREVALENCE))
    print(f"  {'fraction reporting':>18} " + " ".join(f"{v:8.3f}" for v in reporting.mean(0)))
    print(f"  {'fitted prevalence':>18} " + " ".join(f"{v:8.3f}" for v in fitted.mean(0)))
    print(f"  {'fitted sd':>18} " + " ".join(f"{v:8.3f}" for v in fitted.std(0)))
    rhos = [stats.spearmanr(row, SITE_PREVALENCE)[0] for row in fitted]
    exact = np.mean([list(np.argsort(row)) == [0, 1, 2, 3] for row in fitted])
    print(f"  rank correlation within a map: {np.mean(rhos):+.3f} (sd {np.std(rhos):.3f}); "
          f"exact ordering in {exact:.0%} of maps\n", flush=True)
print("The docstring claims ordering within one map survives. A rank correlation near 1 supports")
print("it; one near zero would mean even the ordinal reading is unsupported.")
