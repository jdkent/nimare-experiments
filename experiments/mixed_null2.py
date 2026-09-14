"""Is the mixed relocation null blind to an effect that only the images see?

Holding images fixed while relocating coordinates makes the test conditional: it asks whether
the coordinates add structure beyond the images, not whether there is an effect. If so, an
effect visible only in the images should be undetectable however strong it is -- because it is
present identically in the observed map and in every null iteration, so it cancels.

Two locations, one studyset:
  origin (0,0,0)  -- effect in the images AND reported as peaks by the coordinate studies
  offset (20,0,0) -- effect in the images only; no coordinate study reports there

A test that answers "is there an effect" should find both. A test conditional on the images
should find only the first.
"""
import json, sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from pathlib import Path
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import t_to_z

OUT = Path("/tmp/claude-0/mixednull2"); OUT.mkdir(parents=True, exist_ok=True)
SHAPE = (14, 14, 14)
AFF = np.diag([4.0, 4.0, 4.0, 1.0]); AFF[:3, 3] = -26.0
MASK = str(OUT / "mask.nii.gz")
nib.save(nib.Nifti1Image(np.ones(SHAPE, np.int32), AFF), MASK)
grid = np.indices(SHAPE).astype(float)
ORIGIN = tuple((np.array(SHAPE) - 1) // 2)
OFFSET = (ORIGIN[0] + 5, ORIGIN[1], ORIGIN[2])
N_SIMS = int(sys.argv[1]) if len(sys.argv) > 1 else 8


def blob(centre):
    return np.exp(-sum((grid[i] - centre[i]) ** 2 for i in range(3)) / 8.0)


def build(seed, n_image=10, n_coord=20):
    rng = np.random.default_rng(seed)
    truth = 0.8 * blob(ORIGIN) + 0.8 * blob(OFFSET)     # both locations real, in the images
    studies = []
    for k in range(n_image):
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
    for k in range(n_coord):
        n = int(rng.integers(20, 40))
        # reports the origin only; never the offset, which only the images know about
        pts = [{"space": "MNI", "coordinates": [0.0, 0.0, 0.0],
                "values": [{"kind": "Z", "value": float(max(t_to_z(
                    np.array([0.8 * np.sqrt(n)]), n - 1)[0], 3.3))}]}]
        for _ in range(4):
            xyz = rng.uniform(-24, 24, 3)
            if abs(xyz[0] - 20.0) < 10 and abs(xyz[1]) < 10 and abs(xyz[2]) < 10:
                continue                                  # keep the offset unreported
            pts.append({"space": "MNI", "coordinates": [float(v) for v in xyz],
                        "values": [{"kind": "Z", "value": float(rng.uniform(3.3, 5.0))}]})
        studies.append({"id": f"c{k}", "name": f"c{k}", "metadata": {"sample_sizes": [n]},
            "analyses": [{"id": f"c{k}-1", "name": "1", "metadata": {"sample_sizes": [n]},
                "points": pts, "images": []}]})
    path = OUT / f"ss_{seed}.json"
    path.write_text(json.dumps({"id": "m", "name": "m", "studies": studies}))
    return Studyset(str(path))


print(f"{N_SIMS} simulations; both locations carry a real effect in the images,\n"
      f"but only the origin is ever reported as a coordinate.\n", flush=True)
p_o, p_f, g_o, g_f = [], [], [], []
for seed in range(N_SIMS):
    res = CBES(fwhm=10.0, mask=MASK, use_images=True, null_method="montecarlo",
               n_iters=150, selection_model="none", seed=seed).fit(build(seed))
    p = res.get_map("p").get_fdata()
    g = res.get_map("g").get_fdata()
    p_o.append(float(p[ORIGIN])); p_f.append(float(p[OFFSET]))
    g_o.append(float(g[ORIGIN])); g_f.append(float(g[OFFSET]))

print(f"{'location':>34s} {'median g':>9s} {'median p':>9s} {'detected':>9s}")
print(f"{'origin (images + coordinates)':>34s} {np.median(g_o):9.3f} {np.median(p_o):9.3f} "
      f"{np.mean(np.array(p_o) < 0.05):9.2f}")
print(f"{'offset (images only)':>34s} {np.median(g_f):9.3f} {np.median(p_f):9.3f} "
      f"{np.mean(np.array(p_f) < 0.05):9.2f}")
print("\nBoth locations have the same true effect (g = 0.8) in the images.")
