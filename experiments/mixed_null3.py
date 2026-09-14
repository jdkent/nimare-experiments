"""When does holding the images fixed actually contaminate the null?

The null pools |z| over every voxel, so a few fixed signal voxels barely move its quantiles --
which is why focal image-only effects are still detected. The prediction is that this fails when
the image signal is *widespread*: once a large share of the brain carries it, the pooled null
picks up that signal in every permutation, its tail thickens, and the test loses power against
the very effect it is meant to find.

Signal extent is swept from focal to most-of-the-brain, with the effect size held constant.
Power is read at a location that always carries the effect.
"""
import json, sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from pathlib import Path
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset

OUT = Path("/tmp/claude-0/mixednull3"); OUT.mkdir(parents=True, exist_ok=True)
SHAPE = (14, 14, 14)
AFF = np.diag([4.0, 4.0, 4.0, 1.0]); AFF[:3, 3] = -26.0
MASK = str(OUT / "mask.nii.gz")
nib.save(nib.Nifti1Image(np.ones(SHAPE, np.int32), AFF), MASK)
grid = np.indices(SHAPE).astype(float)
CENTRE = tuple((np.array(SHAPE) - 1) // 2)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 6


def run(width):
    """width sets how much of the volume carries the effect."""
    truth = 0.8 * np.exp(-sum((grid[i] - CENTRE[i]) ** 2 for i in range(3)) / (2 * width**2))
    frac = float(np.mean(truth > 0.2))
    ps, gs = [], []
    for seed in range(N):
        rng = np.random.default_rng(seed)
        studies = []
        for k in range(10):
            n = int(rng.integers(20, 40))
            g = (truth + rng.normal(0, 1 / np.sqrt(n), SHAPE)).astype(np.float32)
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
        for k in range(20):
            n = int(rng.integers(20, 40))
            pts = []
            for _ in range(4):
                xyz = rng.uniform(-24, 24, 3)
                pts.append({"space": "MNI", "coordinates": [float(v) for v in xyz],
                            "values": [{"kind": "Z", "value": float(rng.uniform(3.3, 5.0))}]})
            studies.append({"id": f"c{k}", "name": f"c{k}", "metadata": {"sample_sizes": [n]},
                "analyses": [{"id": f"c{k}-1", "name": "1", "metadata": {"sample_sizes": [n]},
                    "points": pts, "images": []}]})
        path = OUT / f"ss_{seed}_{width}.json"
        path.write_text(json.dumps({"id": "m", "name": "m", "studies": studies}))
        res = CBES(fwhm=10.0, mask=MASK, use_images=True, null_method="montecarlo",
                   n_iters=150, selection_model="none", seed=seed).fit(Studyset(str(path)))
        ps.append(float(res.get_map("p").get_fdata()[CENTRE]))
        gs.append(float(res.get_map("g").get_fdata()[CENTRE]))
    return frac, float(np.median(gs)), float(np.median(ps)), float(np.mean(np.array(ps) < 0.05))


print(f"{N} simulations per row; effect size fixed at g = 0.8 at the centre,\n"
      f"only the spatial extent of the image signal changes.\n", flush=True)
print(f"{'signal extent':>14s} {'% of volume':>12s} {'median g':>9s} {'median p':>9s} "
      f"{'power':>7s}", flush=True)
for width in (1.5, 3.0, 5.0, 8.0):
    frac, g, p, pw = run(width)
    print(f"{width:11.1f}mm {100 * frac:11.0f}% {g:9.3f} {p:9.3f} {pw:7.2f}", flush=True)
