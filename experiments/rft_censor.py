"""Does an RFT extreme-value reporting probability explain the observed rates?

The estimator currently asks whether *this voxel* cleared threshold. The event that actually
occurred is whether the study's field reached threshold anywhere in the neighbourhood it is
scored over -- a maximum, not a point. For a smooth Gaussian field with noncentrality lambda
over a region with R resels,

    E[EC](z) = P(Z > z) + R * A * (z^2 - 1) exp(-z^2 / 2),   A = (4 ln 2)^(3/2) / (2 pi)^2
    P(reach u) = 1 - exp(-2 E[EC](u - lambda))               (two-sided)

with lambda = g sqrt(N). Before plumbing this into the likelihood, check the thing that
matters: can it reproduce how often the pain studies actually reported, where the pointwise
form could not at any scale?
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from scipy.special import ndtr
from load_pain import load_pain
from nimare.meta.utils import sphere_kernel_offsets
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.utils import mm2vox

A_EC = (4.0 * np.log(2.0)) ** 1.5 / (2.0 * np.pi) ** 2

def expected_ec(z, resels):
    """E[EC] of a smooth Gaussian field above z over a region of `resels` resels."""
    z = np.asarray(z, dtype=float)
    return (1.0 - ndtr(z)) + resels * A_EC * (z**2 - 1.0) * np.exp(-0.5 * z**2)

def p_report_rft(g, sqrt_n, u, resels):
    """P(the field reaches u somewhere in the region), two-sided."""
    ec = np.clip(expected_ec(u - g * sqrt_n, resels), 0.0, None)
    return 1.0 - np.exp(-2.0 * ec)

def p_report_pointwise(g, sqrt_n, u, tau=0.1):
    """What the estimator currently uses: one voxel, not a region."""
    sd = np.sqrt(1.0 / sqrt_n**2 + tau**2)
    cut = u / sqrt_n
    return 1.0 - (ndtr((cut - g) / sd) - ndtr((-cut - g) / sd))

U, RADIUS = 3.2905, 20.0
ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
sizes = dict(zip([str(i) for i in ss.ids], ss.sample_sizes()))
coords = ImagesToCoordinates(
    merge_strategy="demolish", z_threshold=U, two_sided=True, remove_subpeaks=True
).transform(ss).coordinates

zooms = masker.mask_img.header.get_zooms()[:3]
offsets = sphere_kernel_offsets(RADIUS, zooms)
voxel_volume = float(np.prod(zooms))
region_volume = len(offsets) * voxel_volume
shape = masker.mask_img.shape[:3]
mask_flat = np.flatnonzero(masker.mask_img.get_fdata().ravel() > 0)
lookup = np.full(int(np.prod(shape)), -1, dtype=np.int64)
lookup[mask_flat] = np.arange(mask_flat.size)

# Observed: fraction of in-mask voxels within RADIUS of any of this study's peaks.
observed = {}
for sid, sub in coords.groupby("id"):
    seen = np.zeros(mask_flat.size, dtype=bool)
    for centre in mm2vox(sub[["x", "y", "z"]].values, masker.mask_img.affine):
        vox = centre + offsets
        inside = np.all((vox >= 0) & (vox < shape), axis=1)
        cols = lookup[np.ravel_multi_index(vox[inside].T, shape)]
        seen[cols[cols >= 0]] = True
    observed[str(sid)] = seen.mean()

rates = np.array([observed[s] for s in sorted(observed)])
sqrt_n = np.array([np.sqrt(sizes[s]) for s in sorted(observed)])
print(f"region: {RADIUS:.0f}mm sphere = {region_volume:.0f} mm^3")
print(f"observed reporting (coverage) rate: mean {rates.mean():.3f}, "
      f"range {rates.min():.3f}-{rates.max():.3f}\n")
print(f"{'g':>5s} {'pointwise':>10s} {'RFT 8mm':>9s} {'RFT 10mm':>9s} {'RFT 12mm':>9s}")
for g in (0.1, 0.2, 0.3, 0.5, 0.8):
    row = [f"{g:5.2f}", f"{p_report_pointwise(g, sqrt_n, U).mean():10.3f}"]
    for fwhm in (8.0, 10.0, 12.0):
        resels = region_volume / fwhm**3
        row.append(f"{p_report_rft(g, sqrt_n, U, resels).mean():9.3f}")
    print(" ".join(row))
print(f"\nobserved mean is {rates.mean():.3f}; the pointwise form cannot reach it at any g,")
print("which is why rate matching ran to the top of its range.")
