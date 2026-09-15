"""Why does the censored likelihood leave the reporting floor in the estimate?

Against held-out HCP subjects, CBES's g came out as truth + offset with the offset close to
u/sqrt(N) -- exactly the censoring floor the zero-inflated model is supposed to remove. Two
mechanisms could explain that, and they are distinguishable.

1. **Covered voxels get no censoring term at all.** ``_coverage_entries`` marks a study as
   "covered" at voxels its foci reach, which suppresses its censoring contribution. So at a
   voxel where studies actually reported, the studies supplying the values contribute no
   truncation correction, and every value they supply is >= u/sqrt(N) by construction. The
   estimate is then the weighted mean of supra-threshold values, which cannot fall below the
   floor. If this is the mechanism, ``selection_model="zero-inflated"`` and ``"none"`` agree at
   well-covered voxels, because the former has nothing to do there.

2. **The prevalence absorbs the silence.** Where studies are silent, the mixture can explain it
   with a low pi rather than a low mu, leaving mu anchored at the reported values. If this is
   the mechanism, prevalence falls where coverage is partial while g stays at the floor.

Measured on the same HCP data and the same held-out truth, stratified by how many studies
report near the voxel.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy.ndimage import maximum_filter
from nilearn.datasets import load_mni152_brain_mask
from nilearn.maskers import NiftiMasker
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import t_to_z

CONTRAST, THRESHOLD, N_PER, N_STUDIES = "MOTOR_LH", 3.2905, 40, 9
mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
S = np.load(f"/tmp/claude-0/hcp/{CONTRAST}_masked.npy")
rng = np.random.default_rng(0)


def hedges(mean, sd, n):
    return (1.0 - 3.0 / (4.0 * (n - 1) - 1.0)) * mean / np.maximum(sd, 1e-9)


order = rng.permutation(S.shape[0])
used, held = order[:N_PER * N_STUDIES], order[N_PER * N_STUDIES:]
truth = np.abs(hedges(S[held].mean(axis=0), S[held].std(axis=0, ddof=1), len(held)))

studies, peak_g = [], []
for k in range(N_STUDIES):
    block = S[used[k * N_PER:(k + 1) * N_PER]]
    mean, sd = block.mean(axis=0), block.std(axis=0, ddof=1)
    z = np.nan_to_num(t_to_z(mean / np.maximum(sd / np.sqrt(N_PER), 1e-9), dof=N_PER - 1))
    volume = np.zeros(shape); volume[mask_bool] = z
    magnitude = np.abs(volume)
    idx = np.argwhere((magnitude == maximum_filter(magnitude, size=3))
                      & (magnitude >= THRESHOLD) & mask_bool)
    if len(idx) < 2:
        continue
    values = volume[tuple(idx.T)]
    meta = {"sample_sizes": [N_PER]}
    studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
        {"id": f"s{k}", "name": "1", "metadata": meta,
         "points": [{"space": "MNI",
                     "coordinates": [float(c) for c in nib.affines.apply_affine(
                         affine, np.asarray(i, dtype=float))],
                     "values": [{"kind": "Z", "value": float(v)}]}
                    for i, v in zip(idx, values)]}]})
    peak_g.append(np.abs(values).mean() / np.sqrt(N_PER))

floor = THRESHOLD / np.sqrt(N_PER)
print(f"{len(studies)} studies, threshold {THRESHOLD}, N {N_PER}")
print(f"censoring floor u/sqrt(N) = {floor:.3f}; mean reported |g| = {np.mean(peak_g):.3f}\n")

fits = {}
for model in ("zero-inflated", "none"):
    est = CBES(fwhm=10.0, mask=masker, peak_bias=None, null_method="none",
               threshold="study-min", selection_model=model)
    r = est.fit(Studyset({"id": "h", "name": "h", "studies": studies}, target=None,
                         mask=mask_img))
    fits[model] = {
        "g": np.abs(r.get_map("g", return_type="array").ravel()),
        "n": r.get_map("n_studies", return_type="array").ravel(),
        "pi": (r.get_map("prevalence", return_type="array").ravel()
               if "prevalence" in r.maps else None),
    }

zi, none = fits["zero-inflated"], fits["none"]
covered = zi["n"] > 0
print(f"{'studies at voxel':>17} {'vox':>7} {'truth':>7} {'zero-infl':>10} {'none':>7} "
      f"{'zi-none':>8} {'zi-truth':>9} {'prev':>6}")
for lo, hi in ((1, 2), (2, 4), (4, 6), (6, 10), (10, 99)):
    use = covered & (zi["n"] >= lo) & (zi["n"] < hi)
    if use.sum() < 50:
        continue
    print(f"{f'{lo}-{hi-1}':>17} {int(use.sum()):>7} {truth[use].mean():7.3f} "
          f"{zi['g'][use].mean():10.3f} {none['g'][use].mean():7.3f} "
          f"{zi['g'][use].mean()-none['g'][use].mean():8.3f} "
          f"{zi['g'][use].mean()-truth[use].mean():9.3f} "
          f"{zi['pi'][use].mean() if zi['pi'] is not None else float('nan'):6.3f}")
use = covered
print(f"{'all covered':>17} {int(use.sum()):>7} {truth[use].mean():7.3f} "
      f"{zi['g'][use].mean():10.3f} {none['g'][use].mean():7.3f} "
      f"{zi['g'][use].mean()-none['g'][use].mean():8.3f} "
      f"{zi['g'][use].mean()-truth[use].mean():9.3f} "
      f"{zi['pi'][use].mean() if zi['pi'] is not None else float('nan'):6.3f}")
print(f"\nif 'zi-none' is ~0 the censoring is doing nothing where studies reported;")
print(f"if 'zi-truth' tracks the floor {floor:.3f}, the estimate is anchored at it.")
