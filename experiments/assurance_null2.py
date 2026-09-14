"""E, properly posed: is the relocation null wrong, or is the *mask* wrong?

The first attempt varied how tightly null foci clump while holding the mask at the full box.
That conflated two different things. `noise_extent` is the half-width of the cube foci are drawn
from, so "clumped 10mm" meant foci occupying 1.6% of the mask -- which tests whether the null
survives a mask far larger than the region foci can actually occupy, not whether it survives
non-uniformity as such.

Separated here, with no true effect anywhere in any condition:

  matched      foci uniform over the whole mask            the null's assumption holds exactly
  mismatched   foci in a small region, mask is the box     assumption broken by the mask
  remasked     the same small region, mask is that region  assumption restored by fixing it

If `mismatched` fails and `remasked` passes, the null is sound and the guidance is that the mask
must be the space foci can occupy -- which is why ALE ships a gray matter mask. If `remasked`
also fails, the problem is the null itself.
"""
import sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES

AFFINE = np.array([[4.0, 0, 0, -40.0], [0, 4.0, 0, -40.0], [0, 0, 4.0, -40.0], [0, 0, 0, 1.0]])
SHAPE = (21, 21, 21)
BIG = nib.Nifti1Image(np.ones(SHAPE, dtype=np.int32), AFFINE)

# voxels whose world coordinate lies within +/- 20 mm, the region the clumped foci occupy
ijk = np.indices(SHAPE).reshape(3, -1).T
world = nib.affines.apply_affine(AFFINE, ijk)
inside = (np.abs(world) <= 20.0).all(axis=1)
small_data = np.zeros(SHAPE, dtype=np.int32)
small_data[tuple(ijk[inside].T)] = 1
SMALL = nib.Nifti1Image(small_data, AFFINE)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 30

CONDITIONS = [
    ("matched    (foci 40mm, mask 40mm)", 40.0, BIG),
    ("mismatched (foci 20mm, mask 40mm)", 20.0, BIG),
    ("remasked   (foci 20mm, mask 20mm)", 20.0, SMALL),
]

print(f"{N} simulations per row; no true effect in any condition.")
print(f"big mask {int(np.prod(SHAPE))} voxels, small mask {int(small_data.sum())} voxels\n")
print(f"{'condition':>36s} {'unc p<.05':>10s} {'unc p<.001':>11s} {'FWE .05':>9s} "
      f"{'covered':>9s} {'verdict':>9s}", flush=True)
print("(rates are over *covered* voxels only)", flush=True)

for label, extent, mask in CONDITIONS:
    u5, u1, fwe, frac_cov = [], [], [], []
    for seed in range(N):
        ss = create_effect_size_coordinate_studyset(
            [(0, 0, 0)], effect_sizes=0.0, n_studies=30, sample_size=(20, 40),
            prevalence=0.0, n_noise_foci=8, noise_extent=extent, seed=seed,
            threshold_z=[2.3263, 3.0902, 3.2905, 4.2649])
        est = CBES(fwhm=12.0, mask=mask, null_method="montecarlo", n_iters=150,
                   threshold="study-min", peak_bias="per-study", seed=seed)
        res = est.fit(ss)
        p = res.get_map("p", return_type="array")
        # Only voxels some study actually speaks to are testable; the rest carry z = 0 and
        # p = 1 by construction. Scoring over the whole mask divides by untestable voxels and
        # makes a clumped configuration look conservative when it is merely sparse -- which is
        # what the first version of this measured.
        covered = res.get_map("n_studies", return_type="array") > 0
        if covered.sum() < 50:
            continue
        u5.append(float(np.mean(p[covered] < 0.05)))
        u1.append(float(np.mean(p[covered] < 0.001)))
        frac_cov.append(float(covered.mean()))
        maps, _, _ = est.correct_fwe_montecarlo(res, vfwe_only=True)
        fwe.append(bool(np.any(maps["logp_level-voxel"] > -np.log10(0.05))))
    rate = float(np.mean(fwe))
    # exact binomial upper bound at 0.05, not a normal approximation on few events
    from scipy import stats as st
    upper = st.binomtest(int(np.sum(fwe)), N, 0.05).pvalue
    ok = rate <= 0.05 or upper > 0.05
    print(f"{label:>36s} {np.mean(u5):10.4f} {np.mean(u1):11.5f} {rate:9.2f} "
          f"{np.mean(frac_cov):9.3f} {'PASS' if ok else 'FAIL':>9s}", flush=True)

print("\n(nominal 0.050 / 0.00100 / 0.05; verdict uses an exact binomial test at 0.05)")
