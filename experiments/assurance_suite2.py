"""Round two: the assumptions a skeptic attacks after the basic validity checks pass.

Round one asked whether the estimator is valid when its assumptions hold. These ask what
happens when they do not, which is the condition every real meta-analysis is in.

  E  The spatial null. The Monte Carlo null relocates foci *uniformly* in the mask. Real foci
     are not uniform -- they concentrate in gray matter and in a handful of much-studied
     regions. If the data cluster for reasons that have nothing to do with an effect, a
     uniform null is too diffuse and everything looks significant. This is the classic
     critique of convergence nulls and the single most dangerous assumption CBES inherits.
  F  Independence (A7). Overlapping samples and re-used datasets are endemic; duplicated
     studies are the extreme case.
  G  Influence. One study reporting far more foci than the rest should not dictate the map.
  H  Degeneracy. Inputs that are legal but pathological should not produce a confident answer.
"""
import sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from nimare.correct import FWECorrector
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import mm2vox

AFFINE = np.array([[4.0, 0, 0, -40.0], [0, 4.0, 0, -40.0], [0, 0, 4.0, -40.0], [0, 0, 0, 1.0]])
MASK = nib.Nifti1Image(np.ones((21, 21, 21), dtype=np.int32), AFFINE)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 12
results = []


def report(name, ok, detail):
    results.append((name, ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def null_set(seed, extent, n_foci=8, n_studies=30):
    """No true effect anywhere; `extent` controls how tightly the noise foci cluster."""
    return create_effect_size_coordinate_studyset(
        [(0, 0, 0)], effect_sizes=0.0, n_studies=n_studies, sample_size=(20, 40),
        prevalence=0.0, n_noise_foci=n_foci, noise_extent=extent, seed=seed,
        threshold_z=[2.3263, 3.0902, 3.2905, 4.2649],
    )


print("\nE. THE SPATIAL NULL UNDER NON-UNIFORM FOCI", flush=True)
print("   (no true effect anywhere; only how tightly the null foci clump changes)", flush=True)
for extent, label in ((36.0, "diffuse 36mm"), (18.0, "clumped 18mm"), (10.0, "tight 10mm")):
    unc, fwe = [], []
    for seed in range(N):
        est = CBES(fwhm=12.0, mask=MASK, null_method="montecarlo", n_iters=150,
                   threshold="study-min", peak_bias="per-study", seed=seed)
        res = est.fit(null_set(seed, extent))
        unc.append(float(np.mean(res.get_map("p", return_type="array") < 0.05)))
        maps, _, _ = est.correct_fwe_montecarlo(res, vfwe_only=True)
        fwe.append(bool(np.any(maps["logp_level-voxel"] > -np.log10(0.05))))
    rate = float(np.mean(fwe))
    ok = rate <= 0.05 + 3 * np.sqrt(0.05 * 0.95 / N)
    report(f"E {label}", ok,
           f"uncorrected {np.mean(unc):.4f} (nominal .05), FWE {rate:.2f} over {N} sims")

print("\nF. NON-INDEPENDENCE", flush=True)
base = null_set(0, 36.0, n_studies=15)
dup = base.copy()
rows = dup.coordinates.copy()
clone = rows.copy()
clone["id"] = clone["id"].astype(str) + "_dup"
dup.coordinates = __import__("pandas").concat([rows, clone], ignore_index=True)
try:
    r_base = CBES(fwhm=12.0, mask=MASK, null_method="none", threshold="study-min",
                  peak_bias="per-study").fit(base)
    r_dup = CBES(fwhm=12.0, mask=MASK, null_method="none", threshold="study-min",
                 peak_bias="per-study").fit(dup)
    n_base = float(np.nanmax(r_base.get_map("n_studies", return_type="array")))
    n_dup = float(np.nanmax(r_dup.get_map("n_studies", return_type="array")))
    z_base = float(np.nanmax(np.abs(r_base.get_map("z", return_type="array"))))
    z_dup = float(np.nanmax(np.abs(r_dup.get_map("z", return_type="array"))))
    # Duplicating every study must not be treated as twice the evidence. It will inflate
    # something -- the question is whether z grows roughly as sqrt(2) (as independent data
    # would) or not at all. Anything at or above sqrt(2) means duplication buys false certainty.
    ratio = z_dup / max(z_base, 1e-9)
    report("F1 duplicated studies do not double confidence", ratio < 1.414,
           f"max|z| {z_base:.2f} -> {z_dup:.2f} (ratio {ratio:.2f}); "
           f"n_studies {n_base:.0f} -> {n_dup:.0f}")
except Exception as exc:
    report("F1 duplicated studies do not double confidence", False, f"{type(exc).__name__}: {exc}")

print("\nG. INFLUENCE OF ONE HEAVY STUDY", flush=True)
balanced = create_effect_size_coordinate_studyset(
    [(0, 0, 0)], effect_sizes=0.6, n_studies=30, sample_size=(20, 40),
    n_noise_foci=2, noise_extent=36.0, seed=3)
heavy = create_effect_size_coordinate_studyset(
    [(0, 0, 0)], effect_sizes=0.6, n_studies=30, sample_size=(20, 40),
    n_noise_foci=2, noise_extent=36.0, seed=3)
h = heavy.coordinates
one = str(h["id"].iloc[0])
extra = h[h["id"] == one]
pd = __import__("pandas")
heavy.coordinates = pd.concat([h] + [extra] * 24, ignore_index=True)
centre = tuple(mm2vox(np.array([[0.0, 0.0, 0.0]]), AFFINE)[0])
g_bal = CBES(fwhm=12.0, mask=MASK, null_method="none",
             peak_bias="per-study").fit(balanced).get_map("g").get_fdata()[centre]
g_hvy = CBES(fwhm=12.0, mask=MASK, null_method="none",
             peak_bias="per-study").fit(heavy).get_map("g").get_fdata()[centre]
shift = abs(g_hvy - g_bal) / max(abs(g_bal), 1e-9)
report("G1 one study with 25x the foci does not take over", shift < 0.25,
       f"g at truth {g_bal:.3f} -> {g_hvy:.3f} ({100*shift:.0f}% shift)")

print("\nH. DEGENERATE INPUTS", flush=True)
for label, kwargs in (("one study", dict(n_studies=1)), ("two studies", dict(n_studies=2))):
    try:
        ss = create_effect_size_coordinate_studyset(
            [(0, 0, 0)], effect_sizes=0.6, sample_size=(20, 40), n_noise_foci=1,
            noise_extent=36.0, seed=1, **kwargs)
        res = CBES(fwhm=12.0, mask=MASK, null_method="none", peak_bias="per-study").fit(ss)
        g = res.get_map("g", return_type="array")
        tau2 = res.get_map("tau2", return_type="array")
        finite = bool(np.all(np.isfinite(g)) and np.all(np.isfinite(tau2)))
        # With fewer than two studies heterogeneity is not estimable and must read zero.
        expected_tau = float(np.nanmax(tau2)) == 0.0 if kwargs["n_studies"] < 2 else True
        report(f"H {label} runs and stays finite", finite and expected_tau,
               f"max|g| {np.nanmax(np.abs(g)):.3f}, max tau2 {np.nanmax(tau2):.4f}")
    except Exception as exc:
        report(f"H {label} runs and stays finite", False, f"{type(exc).__name__}: {exc}")

print("\n" + "=" * 72)
passed = sum(1 for _, ok in results if ok)
print(f"{passed} of {len(results)} checks passed")
for name, ok in results:
    if not ok:
        print(f"  FAILED: {name}")
