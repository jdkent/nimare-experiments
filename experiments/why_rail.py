"""Which end does it rail to, and what does the model predict at zero effect?

The scale can only raise the fitted effect, never lower it. So if the model already predicts
more reporting than was observed at an effect of zero, no scale can match and the search must
collapse. That is a testable statement about the null reporting rate the RFT term implies.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma.effectsize import _coverage_resels, _ec_peak, _rft_censoring_terms
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.meta.utils import sphere_kernel_offsets
from nimare.utils import mm2vox

U = 3.2905
ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
sizes = dict(zip([str(i) for i in ss.ids], ss.sample_sizes()))
coords = ImagesToCoordinates(
    merge_strategy="demolish", z_threshold=U, two_sided=True, remove_subpeaks=True
).transform(ss).coordinates

zooms = masker.mask_img.header.get_zooms()[:3]
offsets = sphere_kernel_offsets(20.0, zooms)
shape = masker.mask_img.shape[:3]
flat = np.flatnonzero(masker.mask_img.get_fdata().ravel() > 0)
lookup = np.full(int(np.prod(shape)), -1, dtype=np.int64)
lookup[flat] = np.arange(flat.size)

seen = {}
for sid, sub in coords.groupby("id"):
    hit = np.zeros(flat.size, dtype=bool)
    for centre in mm2vox(sub[["x", "y", "z"]].values, masker.mask_img.affine):
        vox = centre + offsets
        inside = np.all((vox >= 0) & (vox < shape), axis=1)
        cols = lookup[np.ravel_multi_index(vox[inside].T, shape)]
        hit[cols[cols >= 0]] = True
    seen[str(sid)] = hit

ids = sorted(seen)
stack = np.array([seen[s] for s in ids])
active = stack.any(axis=0)          # the calibration only looks at covered voxels
sqrt_n = np.array([np.sqrt(sizes[s]) for s in ids])

print(f"brain-wide coverage rate:        {stack.mean():.3f}")
print(f"coverage within *active* voxels: {stack[:, active].mean():.3f}   "
      f"({active.sum():,} of {flat.size:,})")
print("\nthe scale can only raise the effect, so the model's rate at g = 0 is a floor:\n")
print(f"{'FWHM':>6s} {'P(report | g=0)':>16s} {'vs active-voxel observed':>26s}")
for fwhm in (8.0, 10.0, 12.0, 14.0, 16.0):
    resels = _coverage_resels(20.0, fwhm)
    peak = _ec_peak(resels)
    p = _rft_censoring_terms(np.zeros(len(ids)), np.full(len(ids), U), sqrt_n, resels, peak)
    null_rate = float((1.0 - p["prob"]).mean())
    verdict = "over-predicts -> rails" if null_rate > stack[:, active].mean() else "has room"
    print(f"{fwhm:6.1f} {null_rate:16.3f} {verdict:>26s}")
print("\nA study reporting nothing real still 'reports' somewhere in a 20mm sphere this often,")
print("under the model. Real studies control their false positives; the model does not know it.")
