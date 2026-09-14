"""Fit every main configuration of CBES and check it runs and is self-consistent.

Not a correctness proof -- a sweep for combinations that crash, silently disagree with their
own null, or produce output that cannot be right. Each cell reports:

    run     did it fit at all
    null    is the null built from the same studies the observed statistic used
    fpr     rejection rate under a global null, which should sit near the nominal 0.05

The null column is the one that matters. A null built from a different set of studies than the
observed statistic is not a null for that statistic, however well calibrated it looks on
coordinate-only data.
"""
import json, itertools, sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from pathlib import Path
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset

OUT = Path("/tmp/claude-0/configmatrix"); OUT.mkdir(parents=True, exist_ok=True)
SHAPE = (12, 12, 12)
AFF = np.diag([4.0, 4.0, 4.0, 1.0]); AFF[:3, 3] = -22.0
MASK = str(OUT / "mask.nii.gz")
nib.save(nib.Nifti1Image(np.ones(SHAPE, np.int32), AFF), MASK)


def build(seed, n_image, n_coord=15):
    """Global null: nothing is true anywhere."""
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


N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
SELECTION = ("zero-inflated", "none")
NULLS = ("approximate", "montecarlo", "none")
DATA = ((0, "coords only"), (8, "mixed w/ images"))

print(f"{N} simulations per cell; global null throughout (nominal p<.05 = 0.050)\n", flush=True)
print(f"{'selection':>14s} {'null':>12s} {'data':>16s} {'peak_bias':>10s} {'result':>28s}",
      flush=True)
for sel, null, (n_img, label), bias in itertools.product(
        SELECTION, NULLS, DATA, (None, "per-study")):
    rates = []
    error = None
    for seed in range(N):
        try:
            kw = dict(n_iters=100) if null == "montecarlo" else {}
            res = CBES(fwhm=10.0, mask=MASK, use_images=True, null_method=null,
                       selection_model=sel, peak_bias=bias, threshold="study-min",
                       seed=seed, **kw).fit(build(seed, n_img))
            rates.append(float(np.mean(res.get_map("p", return_type="array") < 0.05)))
        except Exception as exc:
            error = f"{type(exc).__name__}: {str(exc)[:40]}"
            break
    if error:
        result = f"CRASH {error}"
    elif null == "none":
        result = "no null (p=1 by design)"
    else:
        m = float(np.mean(rates))
        flag = "" if m < 0.08 else f"  <-- INFLATED {m / 0.05:.1f}x"
        result = f"p<.05 = {m:.4f}{flag}"
    print(f"{sel:>14s} {null:>12s} {label:>16s} {str(bias):>10s} {result:>28s}", flush=True)
