"""How much does smoothness vary *within* a brain, not just between studies?

The regional censoring term counts resels in a sphere around one voxel, so the smoothness that
belongs in it is the smoothness there. A single whole-brain FWHM is only adequate if the field
is stationary, and fMRI is not: smoothness varies with tissue, with distance from the edge of
the brain, and with registration. If it varies by a factor of two, a global value over-predicts
reporting in rough regions and under-predicts in smooth ones -- a spatially structured error
that a single global scale cannot absorb, which is the symptom the calibration shows.

Local FWHM from the derivative variance in a moving window, the resels-per-voxel idea.
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
WINDOW = 9  # voxels; ~18mm at 2mm, comparable to the coverage sphere

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

    # Mean squared derivative in a moving window, per axis -> local lambda -> local FWHM.
    lam = np.zeros(data.shape + (3,), dtype=float)
    for axis, spacing in enumerate(zooms):
        diff = np.zeros_like(data)
        slicer = [slice(None)] * 3
        slicer[axis] = slice(0, -1)
        diff[tuple(slicer)] = np.diff(data, axis=axis) / spacing
        lam[..., axis] = uniform_filter(diff**2, size=WINDOW)
    lam = np.clip(lam, 1e-12, None)
    fwhm = np.sqrt(4.0 * np.log(2.0) / lam)
    return np.exp(np.log(fwhm).mean(axis=-1))   # geometric mean over axes

ratios, spreads, globals_ = [], [], []
for sid, path in zip(ss.images["id"].astype(str), ss.images["z"]):
    if path is None:
        continue
    field = local_fwhm(path)
    if field is None:
        continue
    local = field[mask]
    local = local[np.isfinite(local) & (local > 0)]
    if local.size < 1000:
        continue
    lo, hi = np.percentile(local, [5, 95])
    ratios.append(hi / lo)
    spreads.append((sid, np.median(local), lo, hi))
    globals_.append(np.median(local))

print(f"{'study':>14s} {'median':>8s} {'5th pct':>9s} {'95th pct':>9s} {'ratio':>7s}")
for sid, med, lo, hi in spreads:
    print(f"{sid[:14]:>14s} {med:8.2f} {lo:9.2f} {hi:9.2f} {hi / lo:7.2f}")
print(f"\nwithin-brain 95/5 ratio: median {np.median(ratios):.2f}x "
      f"(range {min(ratios):.2f}-{max(ratios):.2f})")
print(f"between-study median FWHM: {np.median(globals_):.2f} mm, "
      f"range {min(globals_):.2f}-{max(globals_):.2f}")
print("\nresels scale as FWHM^-3, so a 2x spread in local FWHM is an 8x spread in resels,")
print("and the reporting probability the censoring term predicts moves with it.")
