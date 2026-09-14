"""Is inference valid when only one or two studies contribute images?

Recovery at 1 and 2 images was measured earlier and was fine, but that used null_method="none",
so the null was never exercised there. Every null validation so far used 0, 5, 8 or 10 images.

Small N is the regime where the image side of the null can under-randomize: sign-flipping gives
2**N distinct configurations, so one image study contributes 2 states and two contribute 4. The
coordinate side still varies freely, but if the image contribution dominates the statistic --
images enter at weight 1 at every voxel, coordinates only near peaks -- a coarse image null
could leave the test miscalibrated.

No true effect anywhere; any departure from 0.05 is a false positive rate problem.
"""
import json, sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from pathlib import Path
from scipy import stats as st
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset

OUT = Path("/tmp/claude-0/fewimages"); OUT.mkdir(parents=True, exist_ok=True)
SHAPE = (12, 12, 12)
AFF = np.diag([4.0, 4.0, 4.0, 1.0]); AFF[:3, 3] = -22.0
MASK = str(OUT / "mask.nii.gz")
nib.save(nib.Nifti1Image(np.ones(SHAPE, np.int32), AFF), MASK)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 20


def build(seed, n_image, n_coord=20):
    rng = np.random.default_rng(seed)
    studies = []
    for k in range(n_image):
        n = int(rng.integers(20, 40))
        nib.save(nib.Nifti1Image(rng.normal(0, 1 / np.sqrt(n), SHAPE).astype(np.float32), AFF),
                 OUT / f"{seed}_{k}_g.nii.gz")
        nib.save(nib.Nifti1Image(np.full(SHAPE, 1.0 / n, np.float32), AFF),
                 OUT / f"{seed}_{k}_v.nii.gz")
        studies.append({"id": f"i{k}", "name": f"i{k}", "metadata": {"sample_sizes": [n]},
            "analyses": [{"id": f"i{k}-1", "name": "1", "metadata": {"sample_sizes": [n]},
                "points": [], "images": [
                    {"url": str(OUT / f"{seed}_{k}_g.nii.gz"), "filename": "g.nii.gz",
                     "space": "MNI", "value_type": "g"},
                    {"url": str(OUT / f"{seed}_{k}_v.nii.gz"), "filename": "v.nii.gz",
                     "space": "MNI", "value_type": "g_var"}]}]})
    for k in range(n_coord):
        n = int(rng.integers(20, 40))
        pts = [{"space": "MNI", "coordinates": [float(v) for v in rng.uniform(-20, 20, 3)],
                "values": [{"kind": "Z", "value": float(rng.uniform(3.3, 5.0))}]}
               for _ in range(4)]
        studies.append({"id": f"c{k}", "name": f"c{k}", "metadata": {"sample_sizes": [n]},
            "analyses": [{"id": f"c{k}-1", "name": "1", "metadata": {"sample_sizes": [n]},
                "points": pts, "images": []}]})
    path = OUT / f"ss_{seed}_{n_image}.json"
    path.write_text(json.dumps({"id": "m", "name": "m", "studies": studies}))
    return Studyset(str(path))


print(f"{N} simulations per cell; global null; nominal p<.05 = 0.050\n", flush=True)
print(f"{'null':>12s} {'images':>7s} {'sign states':>12s} {'unc p<.05':>10s} {'FWE .05':>8s} "
      f"{'verdict':>8s}", flush=True)
for null in ("permute-magnitudes",):
    for n_img in (1, 2, 3, 5):
        unc, fwe = [], []
        for seed in range(N):
            kw = dict(n_iters=200)
            est = CBES(fwhm=10.0, mask=MASK, use_images=True, null_method=null,
                       threshold="study-min", peak_bias="per-study", seed=seed, **kw)
            res = est.fit(build(seed, n_img))
            p = res.get_map("p", return_type="array")
            unc.append(float(np.mean(p < 0.05)))
            maps, _, _ = est.correct_fwe_montecarlo(res, vfwe_only=True)
            fwe.append(bool(np.any(maps["logp_level-voxel"] > -np.log10(0.05))))
        m = float(np.mean(unc))
        # two-sided exact binomial on the pooled voxel-level rate is over-powered by spatial
        # correlation, so judge on the simulation-level spread instead
        se = float(np.std(unc)) / np.sqrt(N)
        ok = abs(m - 0.05) <= max(3 * se, 0.015)
        f = f"{np.mean(fwe):8.2f}" if fwe else f"{'-':>8s}"
        print(f"{null:>12s} {n_img:7d} {2 ** n_img:12d} {m:10.4f} {f} "
              f"{'PASS' if ok else 'FAIL':>8s}", flush=True)
