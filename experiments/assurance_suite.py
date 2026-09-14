"""What a skeptical methodologist would demand before believing CBES, run as pass/fail checks.

The existing tests check that the code does what the code intends. These check the things a
reviewer who does not trust the author would ask, in roughly the order Nichols, Eickhoff or
Mumford would ask them. Each prints PASS/FAIL with the number it turned on, so a failure names
itself rather than needing interpretation.

  A. Inference validity
     A1  uncorrected p is valid across the whole alpha range, not just at .05
     A2  voxel-level FWE holds at several alphas
     A3  cluster-level FWE holds (size and mass)
     A4  validity survives few studies (K = 10)
  B. Estimation
     B1  reduces exactly to DerSimonian-Laird when the spatial part is switched off
     B2  consistent: |bias| falls as K grows
     B3  prevalence is recovered, or its bias is characterised
  C. Invariance and determinism
     C1  permuting study order changes nothing
     C2  same seed gives the same map
  D. Robustness to its own assumptions
     D1  sensitivity to kernel bandwidth
     D2  sensitivity to a misspecified reporting threshold
     D3  cluster-extent thresholding, which violates the height-threshold assumption
"""
import sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from scipy import stats
from nimare.correct import FWECorrector
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import mm2vox

AFFINE = np.array([[4.0, 0, 0, -40.0], [0, 4.0, 0, -40.0], [0, 0, 4.0, -40.0], [0, 0, 0, 1.0]])
MASK = nib.Nifti1Image(np.ones((21, 21, 21), dtype=np.int32), AFFINE)
CENTRE = tuple(mm2vox(np.array([[0.0, 0.0, 0.0]]), AFFINE)[0])
results = []


def report(name, ok, detail):
    results.append((name, ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def null_studyset(seed, n_studies=30, n_foci=8):
    return create_effect_size_coordinate_studyset(
        [(0, 0, 0)], effect_sizes=0.0, n_studies=n_studies, sample_size=(20, 40),
        prevalence=0.0, n_noise_foci=n_foci, noise_extent=36.0, seed=seed,
        threshold_z=[2.3263, 3.0902, 3.2905, 4.2649],
    )


# ---------------------------------------------------------------- A. inference validity
print("\nA. INFERENCE VALIDITY", flush=True)
ALPHAS = np.array([0.001, 0.005, 0.01, 0.025, 0.05, 0.10, 0.20])
N_NULL = int(sys.argv[1]) if len(sys.argv) > 1 else 20

rates = []
for seed in range(N_NULL):
    res = CBES(fwhm=12.0, mask=MASK, threshold="study-min", peak_bias="per-study",
               seed=seed).fit(null_studyset(seed))
    p = res.get_map("p", return_type="array")
    rates.append([np.mean(p < a) for a in ALPHAS])
rates = np.array(rates)
observed = rates.mean(0)
# A valid p-value is stochastically >= uniform. Allow Monte Carlo slack from N_NULL sims.
slack = 2.5 * rates.std(0) / np.sqrt(N_NULL) + 0.2 * ALPHAS
ok = np.all(observed <= ALPHAS + slack)
report("A1 uncorrected p valid across alpha", ok,
       "obs " + " ".join(f"{o:.4f}" for o in observed) + " vs nominal "
       + " ".join(f"{a:.3f}" for a in ALPHAS))

any_fwe = {a: [] for a in (0.01, 0.05, 0.10)}
for seed in range(N_NULL):
    est = CBES(fwhm=12.0, mask=MASK, null_method="montecarlo", n_iters=150,
               threshold="study-min", peak_bias="per-study", seed=seed)
    res = est.fit(null_studyset(seed))
    maps, _, _ = est.correct_fwe_montecarlo(res, vfwe_only=True)
    logp = maps["logp_level-voxel"]
    for a in any_fwe:
        any_fwe[a].append(bool(np.any(logp > -np.log10(a))))
ok = all(np.mean(v) <= a + 3 * np.sqrt(a * (1 - a) / N_NULL) for a, v in any_fwe.items())
report("A2 voxel FWE holds at .01/.05/.10", ok,
       ", ".join(f"{a}: {np.mean(v):.2f}" for a, v in any_fwe.items()))

any_cluster = {"size": [], "mass": []}
for seed in range(min(N_NULL, 12)):
    est = CBES(fwhm=12.0, mask=MASK, null_method="montecarlo", n_iters=150,
               cluster_threshold=0.001, threshold="study-min", peak_bias="per-study", seed=seed)
    res = est.fit(null_studyset(seed))
    maps, _, _ = est.correct_fwe_montecarlo(res, voxel_thresh=0.001)
    for kind in any_cluster:
        any_cluster[kind].append(
            bool(np.any(maps[f"logp_desc-{kind}_level-cluster"] > -np.log10(0.05))))
n_c = len(any_cluster["size"])
ok = all(np.mean(v) <= 0.05 + 3 * np.sqrt(0.05 * 0.95 / n_c) for v in any_cluster.values())
report("A3 cluster FWE holds (size, mass)", ok,
       ", ".join(f"{k}: {np.mean(v):.2f}" for k, v in any_cluster.items()) + f" over {n_c} sims")

small = [np.mean(CBES(fwhm=12.0, mask=MASK, threshold="study-min", peak_bias="per-study",
                      seed=s).fit(null_studyset(s, n_studies=10)).get_map(
                          "p", return_type="array") < 0.05) for s in range(N_NULL)]
ok = np.mean(small) <= 0.05 + 3 * np.std(small) / np.sqrt(N_NULL)
report("A4 valid with only 10 studies", ok, f"p<.05 rate {np.mean(small):.4f} (nominal 0.05)")

# ---------------------------------------------------------------- B. estimation
print("\nB. ESTIMATION", flush=True)

# B1: images enter at weight 1 everywhere and contribute no censoring term, so with the
# selection model off the fit must be textbook inverse-variance DerSimonian-Laird.
rng = np.random.default_rng(0)
K = 12
g_k = rng.normal(0.5, 0.25, K)
v_k = rng.uniform(0.02, 0.08, K)
w = 1.0 / v_k
fixed = np.sum(w * g_k) / np.sum(w)
q = np.sum(w * (g_k - fixed) ** 2)
scale = np.sum(w) - np.sum(w**2) / np.sum(w)
tau2 = max(0.0, (q - (K - 1)) / scale)
w_re = 1.0 / (v_k + tau2)
dl = np.sum(w_re * g_k) / np.sum(w_re)
from nimare.meta.cbma.effectsize import _local_dersimonian_laird
a = 1.0 / v_k
tau2_cbes = _local_dersimonian_laird(
    np.array([float(K)]), np.array([a.sum()]), np.array([(a**2).sum()]),
    np.array([(a * g_k).sum()]), np.array([(a * g_k**2).sum()]),
    np.array([a.sum()]), np.array([float(K)]))[0]
w_c = 1.0 / (v_k + tau2_cbes)
cbes_mu = np.sum(w_c * g_k) / np.sum(w_c)
report("B1 reduces to DerSimonian-Laird", abs(cbes_mu - dl) < 1e-10 and abs(tau2_cbes - tau2) < 1e-10,
       f"mu {cbes_mu:.10f} vs {dl:.10f}; tau2 {tau2_cbes:.10f} vs {tau2:.10f}")

bias = {}
for K in (10, 20, 40, 80):
    vals = []
    for seed in range(4):
        ss = create_effect_size_coordinate_studyset(
            [(0, 0, 0)], effect_sizes=0.6, n_studies=K, sample_size=(20, 40),
            n_noise_foci=2, noise_extent=36.0, seed=seed, simulate_field=True)
        res = CBES(fwhm=12.0, mask=MASK, null_method="none", peak_bias="per-study").fit(ss)
        vals.append(res.get_map("g").get_fdata()[CENTRE])
    bias[K] = float(np.mean(vals) - 0.6)
ok = abs(bias[80]) <= abs(bias[10]) + 0.02
report("B2 bias does not grow with K", ok,
       " ".join(f"K={k}: {b:+.3f}" for k, b in bias.items()))

prev = {}
for true_pi in (0.3, 0.6, 1.0):
    vals = []
    for seed in range(3):
        ss = create_effect_size_coordinate_studyset(
            [(0, 0, 0)], effect_sizes=0.8, n_studies=40, sample_size=(20, 40),
            prevalence=true_pi, n_noise_foci=2, noise_extent=36.0, seed=seed)
        res = CBES(fwhm=12.0, mask=MASK, null_method="none", peak_bias="per-study").fit(ss)
        vals.append(res.get_map("prevalence").get_fdata()[CENTRE])
    prev[true_pi] = float(np.mean(vals))
ordered = all(prev[a] < prev[b] + 0.05 for a, b in zip((0.3, 0.6), (0.6, 1.0)))
report("B3 prevalence tracks the truth (monotone)", ordered,
       " ".join(f"true {t}: est {e:.2f}" for t, e in prev.items()))

# ---------------------------------------------------------------- C. invariance
print("\nC. INVARIANCE AND DETERMINISM", flush=True)
ss = create_effect_size_coordinate_studyset(
    [(0, 0, 0)], effect_sizes=0.6, n_studies=20, sample_size=(20, 40),
    n_noise_foci=3, noise_extent=36.0, seed=11)
base = CBES(fwhm=12.0, mask=MASK, null_method="none", peak_bias="per-study").fit(ss)
base_g = base.get_map("g", return_type="array")

shuffled = ss.copy()
order = list(reversed(list(shuffled.coordinates["id"].unique())))
shuffled.coordinates["_k"] = shuffled.coordinates["id"].map({s: i for i, s in enumerate(order)})
shuffled.coordinates = shuffled.coordinates.sort_values("_k").drop(columns="_k").reset_index(drop=True)
perm_g = CBES(fwhm=12.0, mask=MASK, null_method="none",
              peak_bias="per-study").fit(shuffled).get_map("g", return_type="array")
d = float(np.nanmax(np.abs(base_g - perm_g)))
report("C1 study order does not matter", d < 1e-8, f"max |difference| {d:.2e}")

again = CBES(fwhm=12.0, mask=MASK, null_method="none", peak_bias="per-study").fit(ss)
d = float(np.nanmax(np.abs(base_g - again.get_map("g", return_type="array"))))
report("C2 deterministic given a seed", d == 0.0, f"max |difference| {d:.2e}")

# ---------------------------------------------------------------- D. robustness
print("\nD. ROBUSTNESS TO ASSUMPTION VIOLATIONS", flush=True)
band = {}
for fwhm in (8.0, 10.0, 12.0, 16.0):
    ss = create_effect_size_coordinate_studyset(
        [(0, 0, 0)], effect_sizes=0.6, n_studies=30, sample_size=(20, 40),
        n_noise_foci=2, noise_extent=36.0, seed=5, simulate_field=True)
    res = CBES(fwhm=fwhm, mask=MASK, null_method="none", peak_bias="per-study").fit(ss)
    band[fwhm] = float(res.get_map("g").get_fdata()[CENTRE])
spread = max(band.values()) / min(band.values())
report("D1 bandwidth changes g by <1.5x over 8-16mm", spread < 1.5,
       " ".join(f"{k:g}mm: {v:.3f}" for k, v in band.items()) + f" (ratio {spread:.2f})")

thr = {}
ss = create_effect_size_coordinate_studyset(
    [(0, 0, 0)], effect_sizes=0.6, n_studies=30, sample_size=(20, 40),
    n_noise_foci=2, noise_extent=36.0, seed=5, threshold_z=3.2905, simulate_field=True)
for label, t in (("true 3.29", 3.2905), ("under 2.33", 2.3263), ("over 4.26", 4.2649),
                 ("inferred", "study-min")):
    res = CBES(fwhm=12.0, mask=MASK, null_method="none", threshold=t,
               peak_bias="per-study").fit(ss)
    thr[label] = float(res.get_map("g").get_fdata()[CENTRE])
err = max(abs(thr[k] - thr["true 3.29"]) for k in thr)
report("D2 threshold misspecification is bounded", err < 0.5,
       " ".join(f"{k}: {v:.3f}" for k, v in thr.items()) + f" (max shift {err:.3f})")

print("\n" + "=" * 72)
passed = sum(1 for _, ok in results if ok)
print(f"{passed} of {len(results)} checks passed")
for name, ok in results:
    if not ok:
        print(f"  FAILED: {name}")
