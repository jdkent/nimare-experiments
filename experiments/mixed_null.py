"""Does the relocation null mean anything when some studies contribute images?

_compute_montecarlo_null relocates the focus table and passes _image_studies_ through
untouched. So in every permutation the image studies sit exactly where they were, contributing
the same map, while only the coordinate studies scatter. The null histogram then pools |z| over
all voxels.

Two consequences to check, because they point opposite ways:

  A  Validity. With no effect anywhere, does holding the images fixed still give a valid test?
     Their contribution is identical in observed and null, so it may simply cancel.
  B  Power. With a real effect that the images see, every null iteration contains that signal
     too, so the null's tail is built partly from the thing being tested. If so the test gets
     *more* conservative the more images you add -- the opposite of what a user expects from
     contributing better data.
"""
import json, sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from pathlib import Path
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import t_to_z

OUT = Path("/tmp/claude-0/mixednull"); OUT.mkdir(parents=True, exist_ok=True)
SHAPE = (14, 14, 14)
AFF = np.diag([4.0, 4.0, 4.0, 1.0]); AFF[:3, 3] = -26.0
MASK = str(OUT / "mask.nii.gz")
nib.save(nib.Nifti1Image(np.ones(SHAPE, np.int32), AFF), MASK)
CENTRE = tuple((np.array(SHAPE) - 1) // 2)
N_SIMS = int(sys.argv[1]) if len(sys.argv) > 1 else 8


def build(seed, true_g, n_image, n_coord):
    """Studyset with n_image studies giving images and n_coord giving peaks."""
    rng = np.random.default_rng(seed)
    grid = np.indices(SHAPE).astype(float)
    blob = np.exp(-sum((grid[i] - CENTRE[i]) ** 2 for i in range(3)) / 8.0)
    truth = true_g * blob
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
        pts = []
        # one focus at the truth when there is one, plus scattered noise foci
        if true_g > 0:
            val = float(truth[CENTRE] + rng.normal(0, 1 / np.sqrt(n)))
            if val > 0:
                z = float(t_to_z(np.array([val * np.sqrt(n)]), n - 1)[0])
                pts.append({"space": "MNI", "coordinates": [0.0, 0.0, 0.0],
                            "values": [{"kind": "Z", "value": max(z, 3.3)}]})
        for _ in range(4):
            xyz = rng.uniform(-24, 24, 3)
            pts.append({"space": "MNI", "coordinates": [float(v) for v in xyz],
                        "values": [{"kind": "Z", "value": float(rng.uniform(3.3, 5.0))}]})
        studies.append({"id": f"c{k}", "name": f"c{k}", "metadata": {"sample_sizes": [n]},
            "analyses": [{"id": f"c{k}-1", "name": "1", "metadata": {"sample_sizes": [n]},
                "points": pts, "images": []}]})
    path = OUT / f"ss_{seed}_{true_g}_{n_image}.json"
    path.write_text(json.dumps({"id": "m", "name": "m", "studies": studies}))
    return Studyset(str(path))


def run(true_g, n_image, n_coord):
    unc, hit = [], []
    for seed in range(N_SIMS):
        ss = build(seed, true_g, n_image, n_coord)
        res = CBES(fwhm=10.0, mask=MASK, use_images=True, null_method="montecarlo",
                   n_iters=150, selection_model="none", seed=seed).fit(ss)
        p = res.get_map("p", return_type="array")
        unc.append(float(np.mean(p < 0.05)))
        pm = res.get_map("p")
        hit.append(float(pm.get_fdata()[CENTRE]))
    return float(np.mean(unc)), float(np.median(hit))


print(f"{N_SIMS} simulations per row; the null relocates coordinates and holds images fixed\n")
print(f"{'scenario':>34s} {'p<.05 rate':>11s} {'p at truth':>11s}", flush=True)
print("--- A. no effect anywhere (nominal p<.05 rate = 0.05) ---", flush=True)
for n_img in (0, 5, 10):
    r, h = run(0.0, n_img, 20)
    print(f"{f'{n_img} image + 20 coordinate studies':>34s} {r:11.4f} {h:11.3f}", flush=True)

print("\n--- B. real effect of g = 0.6 at the origin (want small p at truth) ---", flush=True)
for n_img in (0, 5, 10):
    r, h = run(0.6, n_img, 20)
    print(f"{f'{n_img} image + 20 coordinate studies':>34s} {r:11.4f} {h:11.3f}", flush=True)
