"""Does sample-size spread explain why CBES wins on pain and loses on HCP? A real prediction.

Four dials were tested for that disagreement and all four failed: prevalence, the reported
statistic, `tau2`, and reference construction (which moved 2 of 22 points). The identifiability
result supplies a fifth, and this time it is a derivation rather than a guess.

A silence constrains `(pi, mu)` only through `P = pi S(mu) + (1 - pi) S(0)`, and because the
cutoff in sampling-sd units is just the reported statistic again (`c / sigma ~ z`),

    S(0)  = 2 Phi(z) - 1                                 threshold only
    S(mu) = Phi(z - mu sqrt(n)) - Phi(-z - mu sqrt(n))    threshold and mu sqrt(n)

so only *sample-size* spread separates the parameters. On the synthetic bed that moved the
bounded fraction of the profile interval from 0.42 to 0.69 and the fitted prevalence from 0.82
to 0.91, while threshold spread did nothing at all.

**And the HCP bed is built from equal-sized studies** -- `n_per_study` is one number -- while the
NIDM pain collection has real, widely varying sample sizes. So the prediction is that HCP's
verdict is a property of its uniformity, and that spreading the sizes should move CBES toward
images-only without changing the data volume.

Held fixed so the comparison is only about spread:

  * the same total subject budget, so neither arm has more data;
  * the same number of table studies;
  * **the image studies pinned at the same size in both arms**, so the image channel -- and
    therefore the images-only baseline -- is identical and cannot explain a difference;
  * the same held-out subjects for the truth, per split.

Truth is Hedges' g on subjects that made no coordinate, which is the only reference here that is
independent of the selection being measured. No cap on peaks: each study reports whatever
survives its own correction.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import logging; logging.getLogger("nimare").setLevel(logging.ERROR)
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn.datasets import load_mni152_brain_mask
from nilearn.maskers import NiftiMasker
from nimare.meta.cbma.effectsize import CBES
from nimare.studyset import Studyset
from nimare.transforms import t_to_z
from reporting import report_peaks

CONTRAST = os.environ.get("CONTRAST", "MOTOR_LH")
CACHE = f"/tmp/claude-0/hcp/{CONTRAST}_masked.npy"
SCHEME, FOCUS = "cluster", "max"
N_IMAGES, N_TABLES, IMAGE_N = 2, 14, 30
N_SPLITS = int(os.environ.get("NSPLITS", 4))
PROFILE = os.environ.get("PROFILE", "1") == "1"

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
ZOOMS = mask_img.header.get_zooms()[:3]


def hedges(mean, sd, n):
    """Hedges' g from a one-sample mean and sd over subjects."""
    return (1.0 - 3.0 / (4.0 * (n - 1) - 1.0)) * mean / np.maximum(sd, 1e-9)


def table_sizes(kind, budget, count):
    """`count` study sizes summing to `budget`: all equal, or spread over a wide range."""
    if kind == "uniform":
        sizes = np.full(count, budget // count)
    else:
        raw = np.geomspace(10.0, 80.0, count)
        sizes = np.maximum(np.round(raw * budget / raw.sum()), 8).astype(int)
    sizes = sizes.astype(int)
    # Absorb the rounding remainder into the largest study, keeping the budget exact.
    sizes[int(np.argmax(sizes))] += budget - int(sizes.sum())
    return sizes


def build(S, used, sizes):
    """One collection: `N_IMAGES` image studies of IMAGE_N, then tables of the given sizes."""
    studies, donor_g, donor_v, peaks = [], [], [], 0
    at = 0
    for k in range(N_IMAGES + len(sizes)):
        n = IMAGE_N if k < N_IMAGES else int(sizes[k - N_IMAGES])
        block = S[used[at:at + n]]
        at += n
        mean, sd = block.mean(0), block.std(0, ddof=1)
        t = mean / np.maximum(sd / np.sqrt(n), 1e-9)
        z = np.nan_to_num(t_to_z(t, dof=n - 1))
        found, height = report_peaks(z, mask_bool, shape, ZOOMS, scheme=SCHEME, focus=FOCUS)
        meta = {"sample_sizes": [n], "reporting_threshold": float(height)}
        analysis = {"id": f"s{k}", "name": "1", "metadata": meta, "points": [], "images": []}
        if k < N_IMAGES:
            g = hedges(mean, sd, n)
            var = 1.0 / n + g**2 / (2.0 * n)
            donor_g.append(g)
            donor_v.append(var)
            work = f"/tmp/claude-0/spread_{os.getpid()}"
            os.makedirs(work, exist_ok=True)
            gp, vp = f"{work}/s{k}_g.nii.gz", f"{work}/s{k}_v.nii.gz"
            nib.save(masker.inverse_transform(g), gp)
            nib.save(masker.inverse_transform(var), vp)
            analysis["images"] = [
                {"url": gp, "filename": "g", "space": "MNI", "value_type": "g"},
                {"url": vp, "filename": "v", "space": "MNI", "value_type": "g_var"}]
        elif found:
            peaks += len(found)
            analysis["points"] = [
                {"space": "MNI",
                 "coordinates": [float(c) for c in nib.affines.apply_affine(
                     affine, np.asarray(ijk, dtype=float))],
                 "values": [{"kind": "Z", "value": float(value)}]}
                for ijk, value in found]
        else:
            continue
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta,
                        "analyses": [analysis]})
    collection = Studyset({"id": "hcp", "name": "hcp", "studies": studies},
                          target=None, mask=mask_img)
    images_only = np.abs(sum(g / np.maximum(v, 1e-9) for g, v in zip(donor_g, donor_v))
                         / sum(1.0 / np.maximum(v, 1e-9) for v in donor_v))
    return collection, images_only, peaks, len(studies)


def score(est, truth, use):
    """Pearson r, AUC for the top decile, and the magnitude ratio in the top quartile."""
    a, t = est[use], truth[use]
    top = t >= np.percentile(t, 90)
    order = stats.rankdata(a)
    npos, nneg = int(top.sum()), int((~top).sum())
    auc = ((order[top].sum() - npos * (npos + 1) / 2) / (npos * nneg)
           if npos and nneg else np.nan)
    quartile = t >= np.percentile(t, 75)
    return (float(stats.pearsonr(a, t)[0]), float(auc),
            float(np.median(a[quartile] / np.maximum(t[quartile], 1e-6))))


if __name__ == "__main__":
    S = np.load(CACHE)
    budget = N_TABLES * 30
    print(f"{CONTRAST}: {S.shape[0]} subjects, {S.shape[1]} voxels at 4 mm")
    print(f"{N_IMAGES} image studies of {IMAGE_N} + {N_TABLES} tables sharing {budget} "
          f"subjects, {N_SPLITS} splits, profile={PROFILE}")
    for kind in ("uniform", "spread"):
        print(f"  {kind:>8}: sizes {list(table_sizes(kind, budget, N_TABLES))}")
    print()
    rng = np.random.default_rng(0)
    out = {}
    for split in range(N_SPLITS):
        need = N_IMAGES * IMAGE_N + budget
        order = rng.permutation(S.shape[0])
        used, held = order[:need], order[need:]
        truth = np.abs(hedges(S[held].mean(0), S[held].std(0, ddof=1), len(held)))
        for kind in ("uniform", "spread"):
            sizes = table_sizes(kind, budget, N_TABLES)
            collection, images_only, peaks, n_used = build(S, used, sizes)
            est = CBES(mask=masker, null_method="none", threshold="reporting_threshold",
                       interval="profile" if PROFILE else "wald")
            result = est.fit(collection)
            g = np.abs(result.get_map("g", return_type="array").ravel())
            covered = result.get_map("n_studies", return_type="array").ravel() > 0
            use = covered & np.isfinite(g) & (g > 0) & np.isfinite(truth)
            pi = result.get_map("prevalence", return_type="array").ravel()
            bounded = np.nan
            if PROFILE:
                lo = result.get_map("g_lower", return_type="array").ravel()
                hi = result.get_map("g_upper", return_type="array").ravel()
                strong = use & (truth >= np.percentile(truth[use], 90))
                bounded = float((np.isfinite(lo) & np.isfinite(hi))[strong].mean())
            row = dict(zip(("r", "auc", "ratio"), score(g, truth, use)))
            row.update(zip(("r_img", "auc_img", "ratio_img"), score(images_only, truth, use)))
            row["pi"] = float(pi[use].mean())
            row["bounded"] = bounded
            out.setdefault(kind, []).append(row)
            print(f"  split {split + 1} {kind:>8}: {peaks} peaks, {n_used} studies, "
                  f"pi {row['pi']:.3f}, bounded {bounded:.3f}, "
                  f"CBES ratio {row['ratio']:.3f} vs images {row['ratio_img']:.3f}", flush=True)

    print(f"\n{'arm':>10} {'bounded':>8} {'pi':>6} | {'CBES r':>7} {'CBES auc':>9} "
          f"{'CBES ratio':>11} | {'img r':>7} {'img auc':>8} {'img ratio':>10}")
    for kind, rows in out.items():
        def m(key):
            return float(np.mean([r[key] for r in rows]))
        print(f"{kind:>10} {m('bounded'):8.3f} {m('pi'):6.3f} | {m('r'):7.3f} {m('auc'):9.3f} "
              f"{m('ratio'):11.3f} | {m('r_img'):7.3f} {m('auc_img'):8.3f} "
              f"{m('ratio_img'):10.3f}")
