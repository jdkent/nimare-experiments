"""Factorised null: sample each study's local configuration instead of refitting the brain.

At a single voxel the studies are independent under relocation, so a voxel's null statistic
can be simulated directly. For study k with m_k foci, each focus independently lands

  * inside the kernel support of this voxel, with probability |support| / |mask|
  * inside the coverage sphere but outside the kernel, with the remaining sphere probability
  * nowhere near, otherwise

and a focus that lands is a uniformly chosen one of that study's foci, carrying its own g and
variance. The kernel weight of a hit is drawn from the kernel's own value histogram. Nothing
here touches the brain: every quantity comes from the study's foci count and the kernel
geometry, so a million independent draws cost about what one refit costs.

Run against the relocation null built by approx_null.py, comparing where it matters: the tail.
"""
import sys, time, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.meta.cbma.effectsize import _local_dersimonian_laird
from nimare.meta.utils import sphere_kernel_offsets

N_DRAWS = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
N_REAL = int(sys.argv[2]) if len(sys.argv) > 2 else 200

affine = np.array([[4.0, 0, 0, -40.0], [0, 4.0, 0, -40.0], [0, 0, 4.0, -40.0], [0, 0, 0, 1.0]])
mask = nib.Nifti1Image(np.ones((21, 21, 21), dtype=np.int32), affine)
n_mask = int(mask.get_fdata().sum())

ss = create_effect_size_coordinate_studyset(
    [(0, 0, 0)], effect_sizes=0.5, n_studies=30, sample_size=(15, 45), seed=3,
    n_noise_foci=4, noise_extent=36.0, threshold_z=[2.3263, 3.0902, 3.2905, 4.2649],
)

est = CBES(fwhm=10.0, mask=mask, threshold="study-min-corrected", peak_bias="per-study",
           null_method="parametric")
est.fit(ss)                                   # populates the tables the sampler needs
table, roster = est._focus_table_, est._sample_sizes_
thresholds, study_ids = est._thresholds_, list(est._sample_sizes_.index)
n_studies = len(study_ids)

# Kernel support and its weight histogram, exactly as the estimator truncates it.
offsets, kernel_weights = est._kernel_support()
n_kernel = len(kernel_weights)
radius = est.coverage_radius or 2.0 * (est.fwhm or 10.0)
n_sphere = len(sphere_kernel_offsets(radius, mask.header.get_zooms()[:3]))

p_kernel = n_kernel / n_mask
p_sphere = n_sphere / n_mask
print(f"kernel support {n_kernel} voxels (p={p_kernel:.4f}); "
      f"coverage sphere {n_sphere} (p={p_sphere:.4f}); mask {n_mask}")

rng = np.random.default_rng(0)
start = time.perf_counter()

weights = np.zeros((n_studies, N_DRAWS))
g_obs = np.zeros((n_studies, N_DRAWS))
var_obs = np.ones((n_studies, N_DRAWS))
covered = np.zeros((n_studies, N_DRAWS), dtype=bool)

by_study = {str(k): v for k, v in table.groupby("id")}
for position, study_id in enumerate(study_ids):
    sub = by_study.get(str(study_id))
    if sub is None or not len(sub):
        continue                                   # silent everywhere; never covered
    m = len(sub)
    g_k = sub["g"].to_numpy()
    var_k = sub["var_g"].to_numpy()

    # How many of this study's foci land in the kernel support / the coverage sphere.
    n_hit = rng.binomial(m, p_kernel, size=N_DRAWS)
    covered[position] = rng.binomial(m, p_sphere, size=N_DRAWS) > 0

    hit = n_hit > 0
    covered[position] |= hit                       # a kernel hit is inside the sphere too
    if not hit.any():
        continue
    # One hit dominates in practice; take the strongest weight among the hits.
    which = rng.integers(0, m, size=hit.sum())
    drawn = rng.choice(kernel_weights, size=hit.sum())
    weights[position, hit] = drawn
    g_obs[position, hit] = g_k[which]
    var_obs[position, hit] = var_k[which]

# tau^2, from the same moment sums the estimator accumulates.
a = weights / np.where(var_obs > 0, var_obs, 1.0)
sums = dict(
    w=weights.sum(0), w2=(weights**2).sum(0), a=a.sum(0), a2=(a**2).sum(0),
    ag=(a * g_obs).sum(0), ag2=(a * g_obs**2).sum(0),
    w2_over_s2=((weights**2) / np.where(var_obs > 0, var_obs, 1.0)).sum(0),
    n=(weights > 0).sum(0).astype(float),
)
tau2 = _local_dersimonian_laird(sums["w"], sums["a"], sums["a2"], sums["ag"], sums["ag2"],
                                sums["w2_over_s2"], sums["n"])

total_var = var_obs + tau2
pooling = np.where(weights > 0, weights / total_var, 0.0)
denominator = pooling.sum(0)
start_g = np.divide((pooling * g_obs).sum(0), denominator,
                    out=np.zeros(N_DRAWS), where=denominator > 0)

mu, pi, se = est._fit_chunk(
    weights=weights, g_obs=g_obs, var_obs=var_obs, covered=covered, tau2=tau2,
    null_var=np.array([[v] for v in __import__("nimare.meta.cbma.effectsize", fromlist=["x"])
                       .null_effect_variance(roster.values, design=est.design)]),
    cutoffs=np.abs(thresholds.loc[study_ids].to_numpy())[:, None], start=start_g,
)
z = np.divide(mu, se, out=np.zeros_like(mu), where=np.isfinite(se) & (se > 0))
# The relocation null pools every in-mask voxel, and a voxel no relocated focus reached
# contributes |z| = 0. The factorised draws have to be scored the same way or the two are
# not comparable: keeping only covered draws discards exactly the mass at zero.
z_all = np.abs(z)
z_all[denominator <= 0] = 0.0
z_covered = np.abs(z[denominator > 0])
elapsed = time.perf_counter() - start

print(f"\nfactorised null: {N_DRAWS:,} independent draws in {elapsed:.1f}s "
      f"({len(z_covered):,} covered, {100 * len(z_covered) / N_DRAWS:.1f}%)")
print(f"relocation null for comparison cost 84.5s for {N_REAL} iterations\n")
print(f"{'p':>8s} {'|z| relocation':>15s} {'|z| factorised':>15s} {'ratio':>7s}")
reference = {0.05: 2.057, 0.01: 3.096, 0.001: 4.538, 1e-4: 5.111}
for p, ref in reference.items():
    got = float(np.quantile(z_all, 1.0 - p))
    print(f"{p:8.4g} {ref:15.3f} {got:15.3f} {got / ref:7.3f}")
