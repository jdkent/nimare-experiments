"""Does the HCP floor reproduce on other real datasets?

The held-out-HCP design found CBES's ``g`` sitting near ``(threshold + overshoot) / sqrt(N)``
almost regardless of the truth underneath it. That was one contrast of one collection, and the
synthetic studies were all drawn from a single population, so before drawing any conclusion
about the model the same measurement has to be made where the studies are real and different
from each other.

Two corpora, one protocol.

  * **NIDM pain**, 21 studies with full images.
  * **NeuroVault paradigm groups**, separate collections annotated with the same cognitive
    paradigm, one map per collection.

In both, the studies are split in half at random: one half is reduced to thresholded peak
coordinates and handed to CBES, the other half is pooled by inverse variance to give the truth.
Neither half sees the other, so the comparison does not condition the reference on the same
noise that produced the peaks -- the flaw that invalidated every earlier same-map check.

What decides the question is not the overall ratio but its *shape*: if CBES's ``g`` tracks the
truth it should rise across truth strata roughly one-for-one, and if it is reading the
threshold instead it will be flat near ``u / sqrt(N)`` with the strata ratios falling as the
truth rises. The predicted floor is printed alongside so the two can be told apart by eye.
"""
import json, logging, os, subprocess, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from scipy.ndimage import maximum_filter
from nilearn.datasets import load_mni152_brain_mask
from nilearn.image import resample_to_img
from nilearn.maskers import NiftiMasker
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import d_to_g, t_to_d, t_to_z, ImageTransformer
from load_pain import load_pain

THRESHOLDS = [3.2905, 4.5]
MAX_PEAKS = 10
N_SPLITS = 6
MIN_HALF = 5
CACHE = "/tmp/claude-0/paradigms/nii"

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine


def to_g(z, n):
    """Whole map on the subject-level Hedges' g scale."""
    t = np.sign(z) * np.abs(stats.t.isf(stats.norm.sf(np.abs(z)), n - 1))
    t = np.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0)
    return d_to_g(t_to_d(t, n), n)


def peaks_of(z, cut, cap):
    volume = np.zeros(shape, dtype=float)
    volume[mask_bool] = z
    magnitude = np.abs(volume)
    is_peak = (magnitude == maximum_filter(magnitude, size=3)) & (magnitude >= cut) & mask_bool
    idx = np.argwhere(is_peak)
    if not len(idx):
        return []
    values = volume[tuple(idx.T)]
    order = np.argsort(-np.abs(values))[:cap]
    return [(idx[i], float(values[i])) for i in order]


def pooled_truth(maps, sizes):
    """Inverse-variance pooled Hedges' g over the held-out half."""
    g_stack = np.array([to_g(z, n) for z, n in zip(maps, sizes)])
    weights = 1.0 / np.maximum(1.0 / np.asarray(sizes)[:, None]
                               + g_stack**2 / (2.0 * np.asarray(sizes)[:, None]), 1e-9)
    return np.sum(g_stack * weights, axis=0) / np.maximum(weights.sum(axis=0), 1e-12)


def run_cbes(maps, sizes, cut, selection):
    studies, kept = [], 0
    for k, (z, n) in enumerate(zip(maps, sizes)):
        found = peaks_of(z, cut, MAX_PEAKS)
        if len(found) < 2:
            continue
        kept += 1
        meta = {"sample_sizes": [int(n)]}
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}", "name": "1", "metadata": meta, "points": [
                {"space": "MNI",
                 "coordinates": [float(c) for c in nib.affines.apply_affine(
                     affine, np.asarray(ijk, dtype=float))],
                 "values": [{"kind": "Z", "value": value}]} for ijk, value in found]}]})
    if kept < MIN_HALF:
        return None, None, kept
    est = CBES(fwhm=10.0, mask=masker, peak_bias=None, null_method="none",
               threshold="study-min", selection_model=selection)
    result = est.fit(Studyset({"id": "x", "name": "x", "studies": studies},
                              target=None, mask=mask_img))
    g = np.abs(result.get_map("g", return_type="array").ravel())
    covered = result.get_map("n_studies", return_type="array").ravel() > 0
    return g, covered, kept


def report(label, maps, sizes, rng):
    sizes = np.asarray(sizes, dtype=float)
    n = len(maps)
    if n < 2 * MIN_HALF:
        print(f"{label}: {n} studies, too few to split\n", flush=True)
        return
    for cut in THRESHOLDS:
        rows = {name: [] for name in ("truth", "cbes", "none")}
        strata, kept_all = [], []
        for _ in range(N_SPLITS):
            order = rng.permutation(n)
            lo, hi = order[: n // 2], order[n // 2:]
            truth = np.abs(pooled_truth([maps[i] for i in hi], sizes[hi]))
            g, covered, kept = run_cbes([maps[i] for i in lo], sizes[lo], cut, "zero-inflated")
            if g is None:
                continue
            g0, _, _ = run_cbes([maps[i] for i in lo], sizes[lo], cut, "none")
            kept_all.append(kept)
            use = covered & np.isfinite(g) & (g > 0)
            if use.sum() < 100:
                continue
            rows["truth"].append(truth[use].mean())
            rows["cbes"].append(g[use].mean())
            rows["none"].append(g0[use].mean() if g0 is not None else np.nan)
            cells = []
            for a, b in ((0, 50), (50, 75), (75, 90), (90, 99), (99, 100)):
                band = (truth >= np.percentile(truth, a)) & (
                    truth < np.percentile(truth, b) if b < 100 else np.ones_like(truth, bool))
                pick = use & band
                cells.append((pick.sum(), truth[pick].mean() if pick.sum() >= 30 else np.nan,
                              g[pick].mean() if pick.sum() >= 30 else np.nan))
            strata.append(cells)
        if not rows["truth"]:
            print(f"{label} @ {cut:.2f}: no usable split\n", flush=True)
            continue
        mean_n = sizes.mean()
        floor = cut / np.sqrt(mean_n)
        print(f"--- {label} @ threshold {cut:.2f} "
              f"({n} studies, {np.mean(kept_all):.0f} reported per half, "
              f"mean N {mean_n:.0f}, floor u/sqrt(N) = {floor:.3f}) ---")
        print(f"  {'truth stratum':>16} {'truth g':>9} {'CBES g':>9} {'ratio':>7}")
        block = np.array([[c[1:] for c in cells] for cells in strata], dtype=float)
        for j, name in enumerate(("0-50%", "50-75%", "75-90%", "90-99%", "99-100%")):
            t_bar, g_bar = np.nanmean(block[:, j, 0]), np.nanmean(block[:, j, 1])
            print(f"  {name:>16} {t_bar:9.3f} {g_bar:9.3f} {g_bar / max(t_bar, 1e-9):7.2f}")
        t_bar, g_bar = np.mean(rows["truth"]), np.mean(rows["cbes"])
        print(f"  {'all covered':>16} {t_bar:9.3f} {g_bar:9.3f} {g_bar / max(t_bar, 1e-9):7.2f}"
              f"   (selection_model='none' gives {np.nanmean(rows['none']):.3f};"
              f" g - floor = {g_bar - floor:+.3f} against truth {t_bar:.3f})\n", flush=True)


rng = np.random.default_rng(0)

ss = ImageTransformer(target="z").transform(load_pain())
pain_maps, pain_sizes = [], []
for row, n in zip(ss.images.itertuples(), ss.sample_sizes()):
    path = getattr(row, "z", None)
    if path is None or not os.path.isfile(str(path)) or not np.isfinite(float(n)):
        continue
    img = resample_to_img(nib.load(str(path)), mask_img, interpolation="continuous",
                          force_resample=True, copy_header=True)
    pain_maps.append(np.nan_to_num(masker.transform(img).ravel().astype(float)))
    pain_sizes.append(float(n))
report("NIDM pain", pain_maps, pain_sizes, rng)

EXCLUDE = ("none", "other", "null", "rest eyes open", "rest eyes closed", "none / other",
           "resting state")
entries = json.load(open("/tmp/claude-0/paradigms/maps.json"))
by = {}
for m in entries:
    if (m["paradigm"] or "").strip().lower() in EXCLUDE:
        continue
    by.setdefault(m["paradigm"], {}).setdefault(m["collection"], m)
groups = sorted(((p, list(c.values())) for p, c in by.items()), key=lambda kv: -len(kv[1]))


def fetch(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 2000:
        return True
    subprocess.run(["curl", "-sL", "--max-time", "120", "-o", path, url], capture_output=True)
    return os.path.exists(path) and os.path.getsize(path) > 2000


os.makedirs(CACHE, exist_ok=True)
for paradigm, members in groups[:5]:
    maps, sizes = [], []
    for entry in members:
        path = os.path.join(CACHE, f"{entry['image']}.nii.gz")
        if not entry.get("url") or not fetch(entry["url"], path):
            continue
        try:
            img = nib.load(path)
            if img.ndim > 3:
                img = nib.Nifti1Image(np.asarray(img.dataobj)[..., 0], img.affine, img.header)
            img = resample_to_img(img, mask_img, interpolation="continuous",
                                  force_resample=True, copy_header=True)
            data = np.nan_to_num(masker.transform(img).ravel().astype(float))
        except Exception:
            continue
        if not np.isfinite(data).any() or np.allclose(data, 0):
            continue
        n = float(entry["n"])
        maps.append(t_to_z(data, dof=n - 1) if entry["map_type"] == "t" else data)
        sizes.append(n)
    report(paradigm[:40], maps, sizes, rng)
