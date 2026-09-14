"""The assurance checks that need no Monte Carlo null, so they can run beside the slow ones.

Sections B, C, D from round one and F, G, H from round two. Everything here uses
null_method="none" or no fitting at all, so it finishes in minutes rather than the hour the
FWE checks need. Split out so the answers arrive while the relocation-based checks grind.
"""
import warnings; warnings.simplefilter("ignore")
import numpy as np
import pandas as pd
import nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.meta.cbma.effectsize import _local_dersimonian_laird
from nimare.utils import mm2vox

AFFINE = np.array([[4.0, 0, 0, -40.0], [0, 4.0, 0, -40.0], [0, 0, 4.0, -40.0], [0, 0, 0, 1.0]])
MASK = nib.Nifti1Image(np.ones((21, 21, 21), dtype=np.int32), AFFINE)
CENTRE = tuple(mm2vox(np.array([[0.0, 0.0, 0.0]]), AFFINE)[0])
results = []


def report(name, ok, detail):
    results.append((name, ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def fit(ss, **kw):
    kw.setdefault("peak_bias", "per-study")
    return CBES(fwhm=12.0, mask=MASK, null_method="none", **kw).fit(ss)


print("\nB. ESTIMATION", flush=True)
rng = np.random.default_rng(0)
K = 12
g_k, v_k = rng.normal(0.5, 0.25, K), rng.uniform(0.02, 0.08, K)
w = 1.0 / v_k
fixed = np.sum(w * g_k) / np.sum(w)
q = np.sum(w * (g_k - fixed) ** 2)
tau2 = max(0.0, (q - (K - 1)) / (np.sum(w) - np.sum(w**2) / np.sum(w)))
dl = np.sum(g_k / (v_k + tau2)) / np.sum(1.0 / (v_k + tau2))
a = 1.0 / v_k
tau2_c = _local_dersimonian_laird(
    np.array([float(K)]), np.array([a.sum()]), np.array([(a**2).sum()]),
    np.array([(a * g_k).sum()]), np.array([(a * g_k**2).sum()]),
    np.array([a.sum()]), np.array([float(K)]))[0]
mu_c = np.sum(g_k / (v_k + tau2_c)) / np.sum(1.0 / (v_k + tau2_c))
report("B1 reduces exactly to DerSimonian-Laird",
       abs(mu_c - dl) < 1e-12 and abs(tau2_c - tau2) < 1e-12,
       f"mu diff {abs(mu_c - dl):.2e}, tau2 diff {abs(tau2_c - tau2):.2e}")

bias, ratio = {}, {}
for K in (10, 20, 40, 80):
    vals = []
    for seed in range(4):
        ss = create_effect_size_coordinate_studyset(
            [(0, 0, 0)], effect_sizes=0.6, n_studies=K, sample_size=(20, 40),
            n_noise_foci=2, noise_extent=36.0, seed=seed, simulate_field=True)
        vals.append(fit(ss).get_map("g").get_fdata()[CENTRE])
    bias[K] = float(np.mean(vals) - 0.6)
    ratio[K] = float(np.mean(vals) / 0.6)
ok = abs(bias[80]) <= abs(bias[10]) + 0.02
report("B2 bias stabilises rather than growing with K", ok,
       " ".join(f"K={k}: {b:+.3f} (x{ratio[k]:.2f})" for k, b in bias.items()))

prev = {}
for true_pi in (0.3, 0.6, 1.0):
    vals = []
    for seed in range(3):
        ss = create_effect_size_coordinate_studyset(
            [(0, 0, 0)], effect_sizes=0.8, n_studies=40, sample_size=(20, 40),
            prevalence=true_pi, n_noise_foci=2, noise_extent=36.0, seed=seed)
        vals.append(fit(ss).get_map("prevalence").get_fdata()[CENTRE])
    prev[true_pi] = float(np.mean(vals))
ok = prev[0.3] < prev[0.6] < prev[1.0] + 1e-9
report("B3 prevalence is monotone in the truth", ok,
       " ".join(f"true {t}: est {e:.2f}" for t, e in prev.items()))

print("\nC. INVARIANCE AND DETERMINISM", flush=True)
ss = create_effect_size_coordinate_studyset(
    [(0, 0, 0)], effect_sizes=0.6, n_studies=20, sample_size=(20, 40),
    n_noise_foci=3, noise_extent=36.0, seed=11)
base_g = fit(ss).get_map("g", return_type="array")

shuffled = ss.copy()
rows = shuffled.coordinates
order = {s: i for i, s in enumerate(reversed(list(rows["id"].unique())))}
shuffled.coordinates = (rows.assign(_k=rows["id"].map(order))
                        .sort_values("_k").drop(columns="_k").reset_index(drop=True))
d = float(np.nanmax(np.abs(base_g - fit(shuffled).get_map("g", return_type="array"))))
report("C1 study order does not change the map", d < 1e-8, f"max |difference| {d:.2e}")

d = float(np.nanmax(np.abs(base_g - fit(ss).get_map("g", return_type="array"))))
report("C2 repeated fits are identical", d == 0.0, f"max |difference| {d:.2e}")

print("\nD. ROBUSTNESS TO ITS OWN ASSUMPTIONS", flush=True)
band = {}
for fwhm in (8.0, 10.0, 12.0, 16.0):
    s = create_effect_size_coordinate_studyset(
        [(0, 0, 0)], effect_sizes=0.6, n_studies=30, sample_size=(20, 40),
        n_noise_foci=2, noise_extent=36.0, seed=5, simulate_field=True)
    band[fwhm] = float(CBES(fwhm=fwhm, mask=MASK, null_method="none",
                            peak_bias="per-study").fit(s).get_map("g").get_fdata()[CENTRE])
spread = max(band.values()) / min(band.values())
report("D1 kernel bandwidth 8-16mm shifts g by <1.5x", spread < 1.5,
       " ".join(f"{k:g}mm {v:.3f}" for k, v in band.items()) + f" (ratio {spread:.2f})")

s = create_effect_size_coordinate_studyset(
    [(0, 0, 0)], effect_sizes=0.6, n_studies=30, sample_size=(20, 40),
    n_noise_foci=2, noise_extent=36.0, seed=5, threshold_z=3.2905, simulate_field=True)
thr = {}
for label, t in (("true 3.29", 3.2905), ("under 2.33", 2.3263), ("over 4.26", 4.2649),
                 ("inferred", "study-min")):
    thr[label] = float(fit(s, threshold=t).get_map("g").get_fdata()[CENTRE])
err = max(abs(v - thr["true 3.29"]) for v in thr.values())
report("D2 threshold misspecification is bounded", err < 0.5,
       " ".join(f"{k} {v:.3f}" for k, v in thr.items()) + f" (max shift {err:.3f})")

print("\nF. NON-INDEPENDENCE", flush=True)
base = create_effect_size_coordinate_studyset(
    [(0, 0, 0)], effect_sizes=0.0, n_studies=15, sample_size=(20, 40), prevalence=0.0,
    n_noise_foci=8, noise_extent=36.0, seed=0, threshold_z=[2.3263, 3.0902, 3.2905, 4.2649])
dup = base.copy()
rows = dup.coordinates
clone = rows.copy(); clone["id"] = clone["id"].astype(str) + "_dup"
dup.coordinates = pd.concat([rows, clone], ignore_index=True)
z_b = float(np.nanmax(np.abs(fit(base, threshold="study-min").get_map("z", return_type="array"))))
z_d = float(np.nanmax(np.abs(fit(dup, threshold="study-min").get_map("z", return_type="array"))))
r = z_d / max(z_b, 1e-9)
report("F1 duplicating every study does not buy sqrt(2) confidence", r < 1.414,
       f"max|z| {z_b:.2f} -> {z_d:.2f} (ratio {r:.2f}; independent data would give 1.41)")

print("\nG. INFLUENCE OF ONE HEAVY STUDY", flush=True)
bal = create_effect_size_coordinate_studyset(
    [(0, 0, 0)], effect_sizes=0.6, n_studies=30, sample_size=(20, 40),
    n_noise_foci=2, noise_extent=36.0, seed=3)
hvy = bal.copy()
h = hvy.coordinates
extra = h[h["id"] == str(h["id"].iloc[0])]
hvy.coordinates = pd.concat([h] + [extra] * 24, ignore_index=True)
g_b = float(fit(bal).get_map("g").get_fdata()[CENTRE])
g_h = float(fit(hvy).get_map("g").get_fdata()[CENTRE])
shift = abs(g_h - g_b) / max(abs(g_b), 1e-9)
report("G1 one study with 25x the foci does not take over", shift < 0.25,
       f"g at truth {g_b:.3f} -> {g_h:.3f} ({100 * shift:.0f}% shift)")

print("\nH. DEGENERATE INPUTS", flush=True)
for label, n in (("one study", 1), ("two studies", 2)):
    try:
        s = create_effect_size_coordinate_studyset(
            [(0, 0, 0)], effect_sizes=0.6, n_studies=n, sample_size=(20, 40),
            n_noise_foci=1, noise_extent=36.0, seed=1)
        res = fit(s)
        g = res.get_map("g", return_type="array")
        t2 = res.get_map("tau2", return_type="array")
        finite = bool(np.all(np.isfinite(g)) and np.all(np.isfinite(t2)))
        tau_ok = (float(np.nanmax(t2)) == 0.0) if n < 2 else True
        report(f"H {label} runs and stays finite", finite and tau_ok,
               f"max|g| {np.nanmax(np.abs(g)):.3f}, max tau2 {np.nanmax(t2):.4f}")
    except Exception as exc:
        report(f"H {label} runs and stays finite", False, f"{type(exc).__name__}: {exc}")

print("\n" + "=" * 72)
print(f"{sum(1 for _, ok in results if ok)} of {len(results)} checks passed")
for name, ok in results:
    if not ok:
        print(f"  FAILED: {name}")
