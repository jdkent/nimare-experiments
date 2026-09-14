"""Does the approximate null -- the default -- account for image studies at all?

_approximate_null builds its draws from the focus table and the sample sizes. It never
references _image_studies_. But the observed statistic is fitted with the images included. If
that is right, the observed |z| is computed from more studies than the null ever sees, so it is
systematically larger and the p-values are anti-conservative -- in the default configuration,
on exactly the mixed collections the estimator is meant to support.

No true effect anywhere. Any rejection rate above the nominal 0.05 is a false positive rate.
"""
import json, sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from pathlib import Path
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset

OUT = Path("/tmp/claude-0/approxnull"); OUT.mkdir(parents=True, exist_ok=True)
SHAPE = (14, 14, 14)
AFF = np.diag([4.0, 4.0, 4.0, 1.0]); AFF[:3, 3] = -26.0
MASK = str(OUT / "mask.nii.gz")
nib.save(nib.Nifti1Image(np.ones(SHAPE, np.int32), AFF), MASK)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 8


def build(seed, n_image, n_coord=20):
    rng = np.random.default_rng(seed)
    studies = []
    for k in range(n_image):
        n = int(rng.integers(20, 40))
        g = rng.normal(0, 1 / np.sqrt(n), SHAPE).astype(np.float32)   # pure noise
        nib.save(nib.Nifti1Image(g, AFF), OUT / f"{seed}_{k}_g.nii.gz")
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
        pts = [{"space": "MNI", "coordinates": [float(v) for v in rng.uniform(-24, 24, 3)],
                "values": [{"kind": "Z", "value": float(rng.uniform(3.3, 5.0))}]}
               for _ in range(4)]
        studies.append({"id": f"c{k}", "name": f"c{k}", "metadata": {"sample_sizes": [n]},
            "analyses": [{"id": f"c{k}-1", "name": "1", "metadata": {"sample_sizes": [n]},
                "points": pts, "images": []}]})
    path = OUT / f"ss_{seed}_{n_image}.json"
    path.write_text(json.dumps({"id": "m", "name": "m", "studies": studies}))
    return Studyset(str(path))


print(f"{N} simulations; no true effect anywhere; nominal p<.05 rate = 0.050\n", flush=True)
print(f"{'null method':>13s} {'images':>7s} {'p<.05':>8s} {'p<.001':>9s} {'verdict':>16s}",
      flush=True)
for null in ("approximate", "montecarlo"):
    for n_img in (0, 10):
        r5, r1 = [], []
        for seed in range(N):
            kw = dict(n_iters=150) if null == "montecarlo" else {}
            res = CBES(fwhm=10.0, mask=MASK, use_images=True, null_method=null,
                       selection_model="none", seed=seed, **kw).fit(build(seed, n_img))
            p = res.get_map("p", return_type="array")
            r5.append(float(np.mean(p < 0.05)))
            r1.append(float(np.mean(p < 0.001)))
        m5 = float(np.mean(r5))
        verdict = "ok" if m5 < 0.08 else f"INFLATED {m5 / 0.05:.1f}x"
        print(f"{null:>13s} {n_img:7d} {m5:8.4f} {np.mean(r1):9.5f} {verdict:>16s}", flush=True)
