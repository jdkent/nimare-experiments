"""How much does the conditional reference over-state mu? Measured where mu is known.

The real-data reference counts a study as active at a voxel when its own map clears its own cut
there. That selects on the same noise it then measures, so it is biased up, and the size of the
bias decides how to read a CBES/reference ratio near 1: if the reference runs 20% high, a ratio
of 1.0 means CBES runs 20% high too.

The bias cannot be measured on real data, where mu is unknown. It can be measured in the field
simulator, where the true effect at every voxel is recorded, by building exactly the same
reference and comparing it against that truth.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.transforms import t_to_z

TRUE_G = 0.8
N_SIMS = 40
N_STUDIES = 24
U = 3.2905
shape, step = (21, 21, 21), 4.0
affine = np.eye(4)
affine[:3, :3] *= step
affine[:3, 3] = -step * (np.array(shape) - 1) / 2

print(f"true g = {TRUE_G} at the blob centre, {N_SIMS} simulations, {N_STUDIES} studies\n")
print(f"{'active if |z| >':>15} {'true-effect band':>14} {'voxels':>8} "
      f"{'reference':>10} {'true mu there':>14} {'bias':>8}")
for cut in (3.2905, 4.0, 4.5, 5.0):
    refs, truths = [], []
    for seed in range(N_SIMS):
        rng = np.random.default_rng(1000 + seed)
        g_maps, z_maps, true_maps = [], [], []
        for k in range(N_STUDIES):
            n = int(rng.integers(20, 41))
            # Same construction the simulator uses: a smooth field with a blob of known height.
            noise = rng.normal(size=shape)
            from scipy.ndimage import gaussian_filter
            noise = gaussian_filter(noise, sigma=10.0 / 2.355 / step)
            noise /= noise.std()
            grid = np.stack(np.indices(shape), axis=-1) * step + affine[:3, 3]
            blob_sigma = 10.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
            squared = (grid ** 2).sum(axis=-1)
            signal = TRUE_G * np.exp(-squared / (2.0 * blob_sigma ** 2))
            z = t_to_z(signal * np.sqrt(n) + noise, n - 1)
            g_maps.append(signal + noise / np.sqrt(n))   # observed g at each voxel
            z_maps.append(z)
            true_maps.append(signal)
        G, Z, T = np.array(g_maps), np.array(z_maps), np.array(true_maps)
        active = np.abs(Z) >= cut
        n_active = active.sum(axis=(0,)) if active.ndim == 3 else active.sum(axis=0)
        with np.errstate(invalid="ignore"):
            ref = np.where(n_active > 0,
                           (np.abs(G) * active).sum(axis=0) / np.maximum(n_active, 1), np.nan)
        use = (n_active >= 2) & np.isfinite(ref)
        if use.sum() < 20:
            continue
        refs.append((ref[use], T[0][use]))
    if not refs:
        print(f"{cut:15.2f}   too few voxels at any signal level")
        continue
    ref_all = np.concatenate([r for r, _ in refs])
    true_all = np.concatenate([t for _, t in refs])
    # Stratified by the true effect at the voxel. Averaging over every voxel with two active
    # studies is dominated by the ones where there is no effect at all -- conditioning a study
    # on clearing its threshold at a null voxel selects a pure noise excursion, so the
    # reference is large there while the truth is ~0, and the ratio is unbounded. That says
    # the reference cannot be trusted where there is no signal, which is true but useless; the
    # question is how far it runs above a *real* effect.
    bands = [(0.0, 0.05), (0.05, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 10.0)]
    for lo, hi in bands:
        band = (true_all >= lo * TRUE_G) & (true_all < hi * TRUE_G)
        if band.sum() < 50:
            continue
        r_mean, t_mean = ref_all[band].mean(), true_all[band].mean()
        print(f"{cut:15.2f} {f'{lo:.2f}-{hi:.2f} x g':>14} {int(band.sum()):>8} "
              f"{r_mean:10.3f} {t_mean:14.3f} {r_mean / max(t_mean, 1e-9):8.3f}")
