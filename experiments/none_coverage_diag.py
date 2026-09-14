"""Why is selection_model="none" + montecarlo inflated on the sparse-coverage generator?

Hypothesis: with "none" the statistic exists only where a kernel reaches, and every other voxel
enters both the observed map and the null histogram as |z| = 0. The null histogram therefore
carries an atom of mass (1 - c_null) at zero, and a covered voxel's p is

    p(z) = c_null * P(|Z| >= z | covered).

Rejections can only come from covered voxels, so the rate over the whole mask is

    c_obs * P(p < .05) = c_obs * 0.05 / c_null.

That is exactly nominal when c_obs == c_null and inflated by their ratio when it is not. The
relocation null spreads foci uniformly over the mask; the observed foci need not be uniform,
so the two coverage fractions can differ and nothing in the machinery notices.

This measures c_obs and c_null directly and checks the predicted ratio against the realised
rejection rate.
"""
import json, sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from pathlib import Path
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset

OUT = Path("/tmp/claude-0/nonecov"); OUT.mkdir(parents=True, exist_ok=True)
SHAPE = (12, 12, 12)
AFF = np.diag([4.0, 4.0, 4.0, 1.0]); AFF[:3, 3] = -22.0
MASK = str(OUT / "mask.nii.gz")
nib.save(nib.Nifti1Image(np.ones(SHAPE, np.int32), AFF), MASK)


def build(seed, n_coord=15, span=20.0):
    rng = np.random.default_rng(seed)
    studies = []
    for k in range(n_coord):
        n = int(rng.integers(20, 40))
        pts = [{"space": "MNI", "coordinates": [float(v) for v in rng.uniform(-span, span, 3)],
                "values": [{"kind": "Z", "value": float(rng.uniform(3.3, 5.0))}]}
               for _ in range(4)]
        studies.append({"id": f"c{k}", "name": f"c{k}", "metadata": {"sample_sizes": [n]},
            "analyses": [{"id": f"c{k}-1", "name": "1", "metadata": {"sample_sizes": [n]},
                "points": pts, "images": []}]})
    path = OUT / f"ss_{seed}_{span:.0f}.json"
    path.write_text(json.dumps({"id": "m", "name": "m", "studies": studies}))
    return Studyset(str(path))


N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
print(f"{N} global-null fits, coordinates only, fwhm 10 mm on a 12^3 / 4 mm mask.")
print("c_obs  = share of mask where the observed statistic exists")
print("c_null = mean share per relocation; predicted rate = 0.05 * c_obs / c_null\n", flush=True)
print(f"{'selection':>14s} {'span':>5s} {'c_obs':>7s} {'c_null':>7s} {'ratio':>6s} "
      f"{'pred':>7s} {'actual':>7s} {'covered-only':>13s}", flush=True)

for sel in ("none", "zero-inflated"):
    for span in (20.0, 12.0):
        co, cn, actual, cov_only = [], [], [], []
        for seed in range(N):
            est = CBES(fwhm=10.0, mask=MASK, null_method="montecarlo", selection_model=sel,
                       threshold="study-min", peak_bias="per-study", n_iters=100, seed=seed)
            res = est.fit(build(seed, span=span))
            z = res.get_map("z", return_type="array")
            p = res.get_map("p", return_type="array")
            covered = z != 0.0
            co.append(float(covered.mean()))
            # The null histogram's zero bin is the uncovered share, averaged over relocations.
            hist = est.null_distributions_["histweights_corr-none_method-montecarlo"]
            cn.append(1.0 - float(hist[0] / hist.sum()))
            actual.append(float(np.mean(p < 0.05)))
            cov_only.append(float(np.mean(p[covered] < 0.05)) if covered.any() else np.nan)
        r = np.mean(co) / np.mean(cn)
        print(f"{sel:>14s} {span:5.0f} {np.mean(co):7.3f} {np.mean(cn):7.3f} {r:6.2f} "
              f"{0.05 * r:7.4f} {np.mean(actual):7.4f} {np.nanmean(cov_only):13.4f}", flush=True)
