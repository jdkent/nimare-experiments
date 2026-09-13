"""Measure the smoothness of the pain statistic maps, independently of any coordinate.

If the calibrated scale is tuned by choosing a smoothness that makes it look plausible, nothing
has been shown -- rho and smoothness trade off inside the same predicted reporting rate. The
only way to make the test falsifiable is to measure smoothness from something the calibration
never sees, then fix it.

For a smooth Gaussian field, the variance of the spatial derivative fixes the resel size:
lambda_x = Var(dZ/dx) and FWHM_x = sqrt(4 ln 2 / lambda_x). Estimated here per axis from finite
differences of each study's own z map, which uses the images and no coordinates whatsoever.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
import nibabel as nib
from load_pain import load_pain
from nimare.transforms import ImageTransformer

ss = ImageTransformer(target="z").transform(load_pain())
mask = ss.masker.mask_img.get_fdata() > 0
zooms = np.asarray(ss.masker.mask_img.header.get_zooms()[:3], dtype=float)

def fwhm_of(path):
    """Per-axis FWHM from the variance of the normalized spatial derivative."""
    data = nib.load(str(path)).get_fdata()
    if data.shape[:3] != mask.shape:
        return None
    data = np.where(np.isfinite(data), data, 0.0)
    inside = mask & (data != 0)
    if inside.sum() < 1000:
        return None
    # Standardize so the derivative variance is on the unit-field scale.
    values = data[inside]
    data = (data - values.mean()) / values.std()

    out = []
    for axis, spacing in enumerate(zooms):
        diff = np.diff(data, axis=axis)
        keep = np.diff(mask.astype(int), axis=axis) == 0
        keep &= np.take(mask, np.arange(diff.shape[axis]), axis=axis)
        d = diff[keep]
        if d.size < 1000:
            return None
        lam = float(np.var(d)) / spacing**2      # Var(dZ/dx)
        if lam <= 0:
            return None
        out.append(float(np.sqrt(4.0 * np.log(2.0) / lam)))
    return out

rows = []
for sid, path in zip(ss.images["id"].astype(str), ss.images["z"]):
    if path is None:
        continue
    got = fwhm_of(path)
    if got is not None:
        rows.append((sid, got))

per_axis = np.array([r[1] for r in rows])
geometric = np.exp(np.log(per_axis).mean(axis=1))
print(f"{len(rows)} studies with usable z maps\n")
print(f"{'study':>14s} {'FWHM x':>8s} {'y':>7s} {'z':>7s} {'geo mean':>9s}")
for (sid, axes), g in zip(rows, geometric):
    print(f"{sid[:14]:>14s} {axes[0]:8.2f} {axes[1]:7.2f} {axes[2]:7.2f} {g:9.2f}")
print(f"\nacross studies: median {np.median(geometric):.2f} mm, "
      f"range {geometric.min():.2f}-{geometric.max():.2f} mm")
print(f"the estimator assumed 8.00 mm")
