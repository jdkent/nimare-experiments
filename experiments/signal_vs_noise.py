"""Does signal correlation structure predict noise smoothness?

SDM's templates describe how the *signal* covaries between neighbouring voxels, which is what
shapes their anisotropic kernel. The censoring term here needs how the *noise* is smoothed,
which sets how many effectively independent tests a region holds. Both are shaped by the same
physical processes -- applied smoothing, registration blur, anatomy -- so they may share
structure, but they are not the same quantity and one should not be substituted for the other
on the strength of that alone.

Measured directly: a signal-correlation map (across-study correlation of g with its immediate
neighbours) against a noise-smoothness map (local FWHM from derivative variance). If they agree
well, an SDM-style template could stand in for resels; if not, it cannot.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
import nibabel as nib
from scipy.ndimage import uniform_filter
from load_pain import load_pain
from nimare.transforms import ImageTransformer

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
mask = ss.masker.mask_img.get_fdata() > 0
zooms = np.asarray(ss.masker.mask_img.header.get_zooms()[:3], dtype=float)
WINDOW = 9

stack = []
for sid, path in zip(ss.images["id"].astype(str), ss.images["g"]):
    if path is None:
        continue
    data = nib.load(str(path)).get_fdata()
    if data.shape[:3] == mask.shape:
        stack.append(np.where(np.isfinite(data), data, 0.0))
stack = np.array(stack)
print(f"{len(stack)} studies")

# Signal correlation: across studies, correlation of g(v) with each immediate neighbour.
centred = stack - stack.mean(axis=0, keepdims=True)
norm = np.sqrt((centred**2).sum(axis=0))
norm[norm == 0] = np.inf
unit = centred / norm
# At a 2mm lag against ~12mm smoothness the correlation is near 1 everywhere and carries no
# dynamic range, so the lag has to be comparable to the smoothness for the measure to vary.
LAGS = (1, 3, 6, 9)
signals = {}
for lag in LAGS:
    acc = np.zeros(mask.shape)
    count = 0
    for axis in range(3):
        for shift in (lag, -lag):
            acc += (unit * np.roll(unit, shift, axis=axis + 1)).sum(axis=0)
            count += 1
    signals[lag] = acc / count
signal = signals[LAGS[-1]]

# Noise smoothness: local FWHM from the derivative variance of each study's map, averaged.
noise = np.zeros(mask.shape)
for data in stack:
    inside = mask & (data != 0)
    if inside.sum() < 1000:
        continue
    values = data[inside]
    standardized = (data - values.mean()) / values.std()
    lam = np.zeros(data.shape + (3,))
    for axis, spacing in enumerate(zooms):
        diff = np.zeros_like(standardized)
        slicer = [slice(None)] * 3
        slicer[axis] = slice(0, -1)
        diff[tuple(slicer)] = np.diff(standardized, axis=axis) / spacing
        lam[..., axis] = uniform_filter(diff**2, size=WINDOW)
    fwhm = np.sqrt(4.0 * np.log(2.0) / np.clip(lam, 1e-12, None))
    noise += np.log(np.exp(np.log(fwhm).mean(axis=-1)))
noise /= len(stack)

b = noise[mask]
print(f"\n{'lag':>5s} {'lag mm':>7s} {'median r':>9s} {'5-95 pct':>18s} "
      f"{'vs noise':>9s} {'var expl':>9s}")
for lag in LAGS:
    a = signals[lag][mask]
    keep = np.isfinite(a) & np.isfinite(b)
    x, y = a[keep], b[keep]
    r = float(np.corrcoef(x, y)[0, 1])
    lo, hi = np.percentile(x, [5, 95])
    print(f"{lag:5d} {lag * zooms[0]:7.0f} {np.median(x):9.3f} "
          f"{f'{lo:+.3f} to {hi:+.3f}':>18s} {r:9.3f} {r**2:9.1%}")
print("\nlag has to be comparable to the smoothness for the signal measure to carry any"
      "\ndynamic range; a saturated measure cannot correlate with anything.")
