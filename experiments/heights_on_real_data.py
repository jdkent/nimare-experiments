"""Does deleting the reported magnitudes improve the map on a real collection too?

On simulated fields, replacing every reported height with a single constant improved the fitted
map's correlation with a known truth from 0.429 to 0.713 -- the largest accuracy gain measured
anywhere in this program, from deleting an input. That is a strong enough claim to need a real
check: a bed can favour deletion for reasons peculiar to how it generates fields.

So the same contrast on the 21-study NIDM pain collection, with the protocol's design: one half
supplies coordinates extracted the way a paper would print them (multiplicity-corrected, whole
surviving clusters, one focus each, nothing capped, the real height threshold handed to the
estimator through metadata), the other half supplies the inverse-variance-pooled reference, so no
study is on both sides.

Three inputs, identical in every respect but the reported statistic:

  as reported      the table as extracted
  study-flattened  every focus in a study carries that study's own mean reported value
  all-flattened    every focus in the collection carries the collection's mean

Scored on magnitude accuracy and on localisation separately, because a biased but monotone map
can be good at the second while failing the first, and only the second is what a reader uses a
meta-analytic map for.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
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
from nimare.transforms import d_to_g, t_to_d, ImageTransformer
from load_pain import load_pain
from reporting import report_peaks

SCHEME, FOCUS = "cluster", "max"
N_SPLITS = int(os.environ.get("NSPLITS", 10))

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)


def to_g(z, n):
    t = np.sign(z) * np.abs(stats.t.isf(stats.norm.sf(np.abs(z)), n - 1))
    return d_to_g(t_to_d(np.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0), n), n)


def pooled(members, maps, sizes):
    stack, weights = [], []
    for i in members:
        g = to_g(maps[i], sizes[i])
        stack.append(g)
        weights.append(1.0 / np.maximum(1.0 / sizes[i] + g**2 / (2.0 * sizes[i]), 1e-9))
    stack, weights = np.array(stack), np.array(weights)
    return np.sum(stack * weights, axis=0) / np.maximum(weights.sum(axis=0), 1e-12)


def extract(members, maps, sizes):
    """One table per study: the foci a paper would print, and the threshold behind them."""
    out = []
    for i in members:
        foci, height = report_peaks(maps[i], mask_bool, shape, zooms,
                                    scheme=SCHEME, focus=FOCUS)
        if foci:
            out.append((i, float(sizes[i]), float(height), foci))
    return out


def fit(tables, mode):
    """Fit CBES on the extracted tables, with the reported statistics altered as `mode` says."""
    if mode == "all-flattened":
        grand = float(np.mean([abs(v) for _, _, _, foci in tables for _, v in foci]))
    studies = []
    for i, size, height, foci in tables:
        if mode == "study-flattened":
            level = float(np.mean([abs(v) for _, v in foci]))
            used = [(ijk, np.sign(v) * level) for ijk, v in foci]
        elif mode == "all-flattened":
            used = [(ijk, np.sign(v) * grand) for ijk, v in foci]
        else:
            used = foci
        meta = {"sample_sizes": [int(size)], "reporting_threshold": height}
        studies.append({"id": f"c{i}", "name": f"c{i}", "metadata": meta, "analyses": [
            {"id": f"c{i}", "name": "1", "metadata": meta, "points": [
                {"space": "MNI",
                 "coordinates": [float(c) for c in nib.affines.apply_affine(
                     affine, np.asarray(ijk, dtype=float))],
                 "values": [{"kind": "Z", "value": float(v)}]} for ijk, v in used]}]})
    if len(studies) < 2:
        return None
    est = CBES(fwhm=10.0, mask=masker, peak_bias=None, null_method="none",
               threshold="reporting_threshold")
    res = est.fit(Studyset({"id": "h", "name": "h", "studies": studies},
                           target=None, mask=mask_img))
    g = np.abs(res.get_map("g", return_type="array").ravel())
    pi = (res.get_map("prevalence", return_type="array").ravel()
          if "prevalence" in res.maps else np.full(g.size, np.nan))
    covered = res.get_map("n_studies", return_type="array").ravel() > 0
    return g, pi, covered


def score(est, truth, use):
    a, t = est[use], truth[use]
    top = t >= np.percentile(t, 90)
    order = stats.rankdata(a)
    n_pos, n_neg = int(top.sum()), int((~top).sum())
    auc = ((order[top].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
           if n_pos and n_neg else np.nan)
    return (stats.pearsonr(a, t)[0], float(stats.spearmanr(a, t)[0]), float(auc),
            a.mean() / max(t.mean(), 1e-9))


if __name__ == "__main__":
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
    sizes = np.asarray(sizes, dtype=float)
    total = len(maps)
    print(f"NIDM pain: {total} studies, {SCHEME}/{FOCUS} reporting, "
          f"{total // 2} per half over {N_SPLITS} splits\n")

    modes = ("as reported", "study-flattened", "all-flattened")
    rows = {m: [] for m in modes}
    rng = np.random.default_rng(0)
    for _ in range(N_SPLITS):
        order = rng.permutation(total)
        work, hold = order[: total // 2], order[total // 2:]
        truth = np.abs(pooled(hold, maps, sizes))
        tables = extract(list(work), maps, sizes)
        fits = {m: fit(tables, m) for m in modes}
        if any(v is None for v in fits.values()):
            continue
        use = np.logical_and.reduce([v[2] for v in fits.values()]) & np.isfinite(truth)
        if use.sum() < 100:
            continue
        for m in modes:
            g, pi, _ = fits[m]
            r, rho, auc, ratio = score(g, truth, use)
            pi_r = stats.pearsonr(pi[use], truth[use])[0] if np.isfinite(pi[use]).all() else np.nan
            rows[m].append((r, rho, auc, ratio, pi_r))

    print(f"{'height input':18s} {'r(g,truth)':>11s} {'rank r':>8s} {'AUC top':>8s} "
          f"{'g/truth':>8s} {'r(pi,truth)':>12s}")
    for m in modes:
        a = np.array(rows[m])
        if not a.size:
            print(f"{m:18s}   no usable splits")
            continue
        print(f"{m:18s} {a[:,0].mean():11.3f} {a[:,1].mean():8.3f} {a[:,2].mean():8.3f} "
              f"{a[:,3].mean():8.2f} {np.nanmean(a[:,4]):12.3f}")
    base = np.array(rows["as reported"])
    for m in modes[1:]:
        a = np.array(rows[m])
        if a.size and base.size and a.shape == base.shape:
            d = a[:, 0] - base[:, 0]
            p = stats.ttest_rel(a[:, 0], base[:, 0]).pvalue
            print(f"\n  {m:18s} r changes {d.mean():+.4f} "
                  f"(sd {d.std(ddof=1):.4f}, paired p {p:.3f})")
    print("\nOn simulated fields, all-flattened improved r by +0.284 at paired p < 0.001.")
    print("If it does not here, the simulated result was about the bed and not the method.")
