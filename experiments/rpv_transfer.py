"""Would a template resels-per-voxel map transfer between studies?

A template is only worth building if local smoothness factorises as

    local FWHM(v) ~ (per-study level) x (shared spatial pattern)

with the pattern set by anatomy and registration rather than by each study's own pipeline.
Testable without building anything: compute each study's local FWHM map, divide by its own
median to strip the level, and see whether what remains agrees across studies. High agreement
means a template captures real structure; low agreement means local smoothness is a property of
the individual analysis and no template can stand in for it.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
import nibabel as nib
from scipy.ndimage import uniform_filter
from load_pain import load_pain
from nimare.transforms import ImageTransformer

ss = ImageTransformer(target="z").transform(load_pain())
mask = ss.masker.mask_img.get_fdata() > 0
zooms = np.asarray(ss.masker.mask_img.header.get_zooms()[:3], dtype=float)
WINDOW = 9

def local_fwhm(path):
    data = nib.load(str(path)).get_fdata()
    if data.shape[:3] != mask.shape:
        return None
    data = np.where(np.isfinite(data), data, 0.0)
    inside = mask & (data != 0)
    if inside.sum() < 1000:
        return None
    values = data[inside]
    data = (data - values.mean()) / values.std()
    lam = np.zeros(data.shape + (3,))
    for axis, spacing in enumerate(zooms):
        diff = np.zeros_like(data)
        slicer = [slice(None)] * 3
        slicer[axis] = slice(0, -1)
        diff[tuple(slicer)] = np.diff(data, axis=axis) / spacing
        lam[..., axis] = uniform_filter(diff**2, size=WINDOW)
    fwhm = np.sqrt(4.0 * np.log(2.0) / np.clip(lam, 1e-12, None))
    return np.exp(np.log(fwhm).mean(axis=-1))[mask]

maps = []
for sid, path in zip(ss.images["id"].astype(str), ss.images["z"]):
    if path is None:
        continue
    field = local_fwhm(path)
    if field is not None and np.isfinite(field).all():
        maps.append(field)
maps = np.array(maps)
# Strip each study's own level, leaving only the spatial pattern.
shape_only = np.log(maps / np.median(maps, axis=1, keepdims=True))

corr = np.corrcoef(shape_only)
off = corr[np.triu_indices(len(maps), k=1)]
print(f"{len(maps)} studies\n")
print(f"pairwise correlation of the normalized pattern:")
print(f"  median {np.median(off):+.3f}   range {off.min():+.3f} to {off.max():+.3f}")

# Leave-one-out: does a template from the others predict the study held out?
held = []
for i in range(len(maps)):
    others = np.delete(shape_only, i, axis=0).mean(axis=0)
    held.append(np.corrcoef(shape_only[i], others)[0, 1])
held = np.array(held)
print(f"\nleave-one-out, template from the other {len(maps) - 1} studies:")
print(f"  median r {np.median(held):+.3f}   range {held.min():+.3f} to {held.max():+.3f}")
print(f"  variance of the pattern a template would explain: {np.median(held) ** 2:.1%}")
