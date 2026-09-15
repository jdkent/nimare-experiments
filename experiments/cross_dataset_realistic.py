"""The cross-dataset check again, with coordinates produced the way papers produce them.

``cross_dataset_floor`` reached its conclusion from peaks taken at a fixed uncorrected cut and
capped at the strongest ten per study. Both choices are wrong in ways that could have made the
result: nobody reports an uncorrected 3.29, and a cap fixes the count while letting the
effective cut float to wherever the tenth peak landed, which makes the cut a function of the
signal. This redoes the measurement through :mod:`reporting`, which corrects for multiplicity,
keeps every local maximum 8 mm apart that survives, and caps nothing.

Three reporting schemes, because real papers do not agree on one: voxelwise FDR at q = 0.05,
voxelwise family-wise error by Bonferroni, and a p < 0.001 cluster-forming cut kept only where
the cluster reaches ten voxels. A study whose map has nothing surviving reports nothing and
drops out of the meta-analysis, which is what happens in practice and which the capped
extraction could never produce.

The real height threshold is handed to the estimator as metadata rather than inferred from the
smallest reported value, so ``prevalence`` is read under the conditions its documentation asks
for.

Truth is unchanged and still independent: the studies are split in half, one half reduced to
coordinates, the other half pooled by inverse variance.
"""
import json, logging, os, subprocess, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn.datasets import load_mni152_brain_mask
from nilearn.image import resample_to_img
from nilearn.maskers import NiftiMasker
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import d_to_g, t_to_d, t_to_z, ImageTransformer
from load_pain import load_pain
from reporting import report_peaks

SCHEMES = ("fdr", "fwe", "cluster")
N_SPLITS = 6
MIN_HALF = 5
CACHE = "/tmp/claude-0/paradigms/nii"

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
zooms = mask_img.header.get_zooms()[:3]


def to_g(z, n):
    t = np.sign(z) * np.abs(stats.t.isf(stats.norm.sf(np.abs(z)), n - 1))
    return d_to_g(t_to_d(np.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0), n), n)


def pooled_truth(maps, sizes):
    sizes = np.asarray(sizes, dtype=float)[:, None]
    g_stack = np.array([to_g(z, float(n)) for z, n in zip(maps, sizes.ravel())])
    weights = 1.0 / np.maximum(1.0 / sizes + g_stack**2 / (2.0 * sizes), 1e-9)
    return np.sum(g_stack * weights, axis=0) / np.maximum(weights.sum(axis=0), 1e-12)


def fit(maps, sizes, scheme, selection="zero-inflated"):
    """Coordinates the way a paper would print them, then CBES on nothing else."""
    studies, counts, heights = [], [], []
    for k, (z, n) in enumerate(zip(maps, sizes)):
        found, height = report_peaks(z, mask_bool, shape, zooms, scheme=scheme)
        counts.append(len(found))
        if len(found) < 2:
            continue           # nothing survived: this paper reports no table
        heights.append(height)
        meta = {"sample_sizes": [int(n)], "reporting_threshold": float(height)}
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta, "analyses": [
            {"id": f"s{k}", "name": "1", "metadata": meta, "points": [
                {"space": "MNI",
                 "coordinates": [float(c) for c in nib.affines.apply_affine(
                     affine, np.asarray(ijk, dtype=float))],
                 "values": [{"kind": "Z", "value": v}]} for ijk, v in found]}]})
    if len(studies) < MIN_HALF:
        return None
    est = CBES(fwhm=10.0, mask=masker, peak_bias=None, null_method="none",
               threshold="reporting_threshold", selection_model=selection)
    res = est.fit(Studyset({"id": "x", "name": "x", "studies": studies},
                           target=None, mask=mask_img))
    g = np.abs(res.get_map("g", return_type="array").ravel())
    pi = (res.get_map("prevalence", return_type="array").ravel()
          if "prevalence" in res.maps else np.ones_like(g))
    covered = res.get_map("n_studies", return_type="array").ravel() > 0
    return {"g": g, "pi": pi, "covered": covered, "n_studies": len(studies),
            "peaks": float(np.mean([c for c in counts if c >= 2]) if studies else 0),
            "silent": int(sum(c < 2 for c in counts)),
            "height": float(np.median(heights))}


def report(label, maps, sizes, rng):
    sizes = np.asarray(sizes, dtype=float)
    n = len(maps)
    if n < 2 * MIN_HALF:
        print(f"{label}: {n} studies, too few to split\n", flush=True)
        return
    for scheme in SCHEMES:
        acc = {k: [] for k in ("truth", "g", "marg", "r_g", "r_marg", "strata",
                               "kept", "peaks", "silent", "height")}
        for _ in range(N_SPLITS):
            order = rng.permutation(n)
            lo, hi = order[: n // 2], order[n // 2:]
            truth = np.abs(pooled_truth([maps[i] for i in hi], sizes[hi]))
            got = fit([maps[i] for i in lo], sizes[lo], scheme)
            if got is None:
                continue
            g, pi, use = got["g"], got["pi"], got["covered"]
            use = use & np.isfinite(g) & (g > 0)
            if use.sum() < 100:
                continue
            marg = pi * g
            for key, value in (("kept", got["n_studies"]), ("peaks", got["peaks"]),
                               ("silent", got["silent"]), ("height", got["height"])):
                acc[key].append(value)
            acc["truth"].append(truth[use].mean())
            acc["g"].append(g[use].mean())
            acc["marg"].append(marg[use].mean())
            acc["r_g"].append(stats.pearsonr(g[use], truth[use])[0])
            acc["r_marg"].append(stats.pearsonr(marg[use], truth[use])[0])
            cells = []
            for a, b in ((0, 50), (50, 75), (75, 90), (90, 99), (99, 100)):
                band = (truth >= np.percentile(truth, a)) & (
                    truth < np.percentile(truth, b) if b < 100 else np.ones_like(truth, bool))
                pick = use & band
                ok = pick.sum() >= 30
                cells.append((truth[pick].mean() if ok else np.nan,
                              g[pick].mean() if ok else np.nan,
                              marg[pick].mean() if ok else np.nan))
            acc["strata"].append(cells)
        if not acc["truth"]:
            print(f"{label} / {scheme}: no usable split\n", flush=True)
            continue
        print(f"--- {label} / {scheme} ({n} studies, {np.mean(acc['kept']):.0f} with a table, "
              f"{np.mean(acc['silent']):.1f} reporting nothing, "
              f"{np.mean(acc['peaks']):.0f} peaks each, "
              f"median height z = {np.mean(acc['height']):.2f}) ---")
        print(f"  {'truth stratum':>16} {'truth g':>9} {'CBES g':>9} {'ratio':>7} "
              f"{'pi*g':>8} {'ratio':>7}")
        block = np.array(acc["strata"], dtype=float)
        for j, name in enumerate(("0-50%", "50-75%", "75-90%", "90-99%", "99-100%")):
            t, g_, m = (np.nanmean(block[:, j, c]) for c in range(3))
            print(f"  {name:>16} {t:9.3f} {g_:9.3f} {g_ / max(t, 1e-9):7.2f} "
                  f"{m:8.3f} {m / max(t, 1e-9):7.2f}")
        t, g_, m = np.mean(acc["truth"]), np.mean(acc["g"]), np.mean(acc["marg"])
        print(f"  {'all covered':>16} {t:9.3f} {g_:9.3f} {g_ / max(t, 1e-9):7.2f} "
              f"{m:8.3f} {m / max(t, 1e-9):7.2f}")
        print(f"  correlation with the truth: g {np.mean(acc['r_g']):+.3f}, "
              f"pi*g {np.mean(acc['r_marg']):+.3f}\n", flush=True)


if __name__ == "__main__":
    rng = np.random.default_rng(0)

    ss = ImageTransformer(target="z").transform(load_pain())
    maps, sizes = [], []
    for row, n in zip(ss.images.itertuples(), ss.sample_sizes()):
        path = getattr(row, "z", None)
        if path is None or not os.path.isfile(str(path)) or not np.isfinite(float(n)):
            continue
        img = resample_to_img(nib.load(str(path)), mask_img, interpolation="continuous",
                              force_resample=True, copy_header=True)
        maps.append(np.nan_to_num(masker.transform(img).ravel().astype(float)))
        sizes.append(float(n))
    report("NIDM pain", maps, sizes, rng)

    EXCLUDE = ("none", "other", "null", "rest eyes open", "rest eyes closed", "none / other",
               "resting state")
    by = {}
    for m in json.load(open("/tmp/claude-0/paradigms/maps.json")):
        if (m["paradigm"] or "").strip().lower() in EXCLUDE:
            continue
        by.setdefault(m["paradigm"], {}).setdefault(m["collection"], m)
    groups = sorted(((p, list(c.values())) for p, c in by.items()), key=lambda kv: -len(kv[1]))

    def fetch(url, path):
        if os.path.exists(path) and os.path.getsize(path) > 2000:
            return True
        subprocess.run(["curl", "-sL", "--max-time", "120", "-o", path, url],
                       capture_output=True)
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
                    img = nib.Nifti1Image(np.asarray(img.dataobj)[..., 0], img.affine,
                                          img.header)
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
