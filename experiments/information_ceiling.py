"""Is the compression the estimator's fault, or is it not in the coordinates to begin with?

Seven interventions have now failed to move it: a truncated-normal selection correction, the
per-study peak-bias rescaling, subtracting the censoring floor, three reporting schemes, two
focus conventions, oracle cluster-extent coverage, and an assumed cluster size driving the
kernel. When that many independent changes move nothing, the suspect is the input.

This measures the input directly, bypassing the estimator entirely. For every focus a paper
would print, take two numbers:

  * the effect size the table reports, converted to Hedges' g -- the only magnitude information
    the estimator ever receives about that location;
  * the held-out truth at that same location, from studies whose data played no part in
    producing the focus.

If the reported magnitudes track the held-out truth across foci, the information is present and
the estimator is failing to use it. If they do not, no procedure reading those numbers can
recover the truth, and the compression is a property of coordinate tables rather than of CBES.

The comparison is made against two references for honesty: the held-out pooled truth, and the
reporting study's *own* map at the focus, which the value trivially agrees with and which shows
what a perfect measurement of the wrong thing looks like.
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
from nimare.transforms import d_to_g, t_to_d, ImageTransformer
from load_pain import load_pain
from reporting import report_peaks

N_SPLITS = 12
mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape = mask_img.shape
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)
flat = np.flatnonzero(mask_bool.ravel())
position = np.full(mask_bool.size, -1, dtype=np.int64)
position[flat] = np.arange(flat.size)


def to_g(z, n):
    t = np.sign(z) * np.abs(stats.t.isf(stats.norm.sf(np.abs(z)), n - 1))
    return d_to_g(t_to_d(np.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0), n), n)


def pooled_truth(maps, sizes):
    sizes = np.asarray(sizes, dtype=float)[:, None]
    stack = np.array([to_g(z, float(n)) for z, n in zip(maps, sizes.ravel())])
    weights = 1.0 / np.maximum(1.0 / sizes + stack**2 / (2.0 * sizes), 1e-9)
    return np.sum(stack * weights, axis=0) / np.maximum(weights.sum(axis=0), 1e-12)


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
n = len(maps)
print(f"NIDM pain: {n} studies\n", flush=True)

rng = np.random.default_rng(0)
for scheme, focus in (("cluster", "max"), ("cluster", "com"), ("fdr", "max")):
    reported, held, own = [], [], []
    for _ in range(N_SPLITS):
        order = rng.permutation(n)
        lo, hi = order[: n // 2], order[n // 2:]
        truth = np.abs(pooled_truth([maps[i] for i in hi], sizes[hi]))
        for i in lo:
            foci, _ = report_peaks(maps[i], mask_bool, shape, zooms,
                                   scheme=scheme, focus=focus)
            if not foci:
                continue
            own_g = np.abs(to_g(maps[i], sizes[i]))
            for ijk, value in foci:
                v = position[np.ravel_multi_index(tuple(np.asarray(ijk)), shape)]
                if v < 0:
                    continue
                reported.append(abs(float(to_g(np.array([value]), sizes[i])[0])))
                held.append(truth[v])
                own.append(own_g[v])
    reported, held, own = map(np.asarray, (reported, held, own))
    if reported.size < 50:
        print(f"{scheme}/{focus}: too few foci\n")
        continue
    print(f"--- {scheme} reporting, focus at the {focus}: {reported.size} foci ---")
    print(f"  reported g at the focus:  mean {reported.mean():.3f}, "
          f"sd {reported.std():.3f}, 5th-95th {np.percentile(reported, 5):.3f}"
          f"-{np.percentile(reported, 95):.3f}")
    print(f"  held-out truth there:     mean {held.mean():.3f}, "
          f"sd {held.std():.3f}, 5th-95th {np.percentile(held, 5):.3f}"
          f"-{np.percentile(held, 95):.3f}")
    print(f"  correlation of the reported value with the held-out truth: "
          f"{stats.pearsonr(reported, held)[0]:+.3f}")
    print(f"  correlation with the reporting study's own map there:      "
          f"{stats.pearsonr(reported, own)[0]:+.3f}")
    # Regression of truth on the reported value: the slope is how much of a unit of reported
    # magnitude is real, and the intercept is what a table gives when the truth is zero.
    slope, intercept, r, _, _ = stats.linregress(reported, held)
    print(f"  truth = {slope:.3f} * reported {intercept:+.3f}   "
          f"(a slope near 0 means the value carries no magnitude)\n", flush=True)
