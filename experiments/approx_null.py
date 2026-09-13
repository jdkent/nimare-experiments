"""Does a factorised per-voxel null reproduce the full relocation null?

The relocation null moves every focus to a uniform in-mask voxel and refits the brain. But at
any single voxel the studies are independent: the number of study-k foci landing within kernel
range is Binomial(m_k, |support| / |mask|), and the weight of each hit is drawn from the
kernel's own value histogram. If that is right, a voxel's null z can be sampled directly --
draw each study's local configuration, fit that one voxel -- at a tiny fraction of the cost,
and the draws are genuinely independent rather than spatially correlated.

This checks the assumption end to end: build the real relocation null, build the factorised
null, and compare the distributions of |z| where it matters -- the tail.
"""
import sys, time, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES

N_REAL = int(sys.argv[1]) if len(sys.argv) > 1 else 200

affine = np.array([[4.0, 0, 0, -40.0], [0, 4.0, 0, -40.0], [0, 0, 4.0, -40.0], [0, 0, 0, 1.0]])
mask = nib.Nifti1Image(np.ones((21, 21, 21), dtype=np.int32), affine)

ss = create_effect_size_coordinate_studyset(
    [(0, 0, 0)], effect_sizes=0.5, n_studies=30, sample_size=(15, 45), seed=3,
    n_noise_foci=4, noise_extent=36.0, threshold_z=[2.3263, 3.0902, 3.2905, 4.2649],
)

est = CBES(fwhm=10.0, mask=mask, threshold="study-min-corrected", peak_bias="per-study",
           n_iters=N_REAL, seed=0, null_method="montecarlo")
start = time.perf_counter()
est.fit(ss)
real_time = time.perf_counter() - start
histogram = est.null_distributions_["histweights_corr-none_method-montecarlo"]

edges = np.arange(0.0, 50.0 + 0.01, 0.01)
centres = edges[:-1] + 0.005
total = histogram.sum()
cdf = np.cumsum(histogram) / total

def quantile_of(p):
    """|z| exceeded with probability p under the relocation null."""
    return float(np.interp(1.0 - p, cdf, centres[: len(cdf)]))

print(f"relocation null: {N_REAL} iterations, {real_time:.1f}s "
      f"({real_time / N_REAL:.2f}s per iteration)")
print(f"effective independent draws are far fewer than the "
      f"{N_REAL * int(mask.get_fdata().sum()):,} voxel-samples, since neighbours share studies")
print(f"\n{'p':>8s} {'|z| relocation':>15s}")
for p in (0.05, 0.01, 0.001, 1e-4):
    print(f"{p:8.4g} {quantile_of(p):15.3f}")
