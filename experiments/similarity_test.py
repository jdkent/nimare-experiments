"""Does spatial similarity predict effect magnitude?

The question that decides whether a similarity-weighted mixture over a reference corpus can
supply the effect-size scale coordinates cannot identify. If images that look alike have
similar magnitudes, weighting by similarity beats the corpus average and the mixture is worth
building. If not, it collapses to that average and the simpler borrowed-constant approach is
all that is on offer.

Each map is reduced to a vector on a common 4 mm grid and converted to Hedges' g, so magnitude
is comparable across studies of different sample size. Similarity is the spatial correlation
between two maps; magnitude is the mean |g| over voxels in the top decile of that map.
"""
import json, os, subprocess, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from nilearn.image import resample_to_img
from nimare.transforms import t_to_z, z_to_t
from nimare.utils import get_template

OUT = "/tmp/claude-0/corpus"
manifest = json.load(open(f"{OUT}/manifest.json"))

# Common 4 mm grid, small enough that 260 maps fit comfortably in memory.
template = get_template(space="mni152_2mm", mask="brain")
affine = template.affine.copy()
affine[:3, :3] *= 2.0
shape = tuple(int(np.ceil(s / 2)) for s in template.shape[:3])
target = nib.Nifti1Image(np.ones(shape, dtype=np.int8), affine)
grid = resample_to_img(template, target, interpolation="nearest", force_resample=True,
                       copy_header=True)
mask = grid.get_fdata() > 0
print(f"common grid: {int(mask.sum())} in-mask voxels at 4 mm", flush=True)

def hedges_g(path, map_type, n):
    """Load a t or z map and put it on the Hedges' g scale."""
    img = nib.load(path)
    if img.ndim > 3:
        img = nib.Nifti1Image(img.get_fdata()[..., 0], img.affine)
    img = resample_to_img(img, grid, interpolation="continuous", force_resample=True,
                          copy_header=True)
    data = img.get_fdata()
    values = np.where(np.isfinite(data), data, 0.0)[mask]
    if not np.any(values):
        return None
    t = values if map_type == "t" else z_to_t(values, n - 1)
    bias = 1.0 - 3.0 / (4.0 * (n - 1) - 1)
    return (t / np.sqrt(n)) * bias

vectors, meta = [], []
for i, rec in enumerate(manifest):
    path = f"{OUT}/img_{rec['image']}.nii.gz"
    if not os.path.exists(path):
        subprocess.run(["curl", "-sL", "--max-time", "90", "-o", path, rec["url"]],
                       capture_output=True)
    if not os.path.exists(path) or os.path.getsize(path) < 2000:
        continue
    try:
        g = hedges_g(path, rec["map_type"], rec["n"])
    except Exception:
        g = None
    os.remove(path)                      # keep the disk footprint flat
    if g is None or not np.isfinite(g).all():
        continue
    vectors.append(g.astype(np.float32))
    meta.append(rec)
    if len(vectors) % 25 == 0:
        print(f"  {len(vectors)} maps loaded ({i + 1} tried)", flush=True)

X = np.array(vectors)
np.save(f"{OUT}/g_vectors.npy", X)
json.dump(meta, open(f"{OUT}/g_meta.json", "w"))
print(f"\n{len(X)} maps on the common grid", flush=True)
