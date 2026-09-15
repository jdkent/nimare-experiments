"""Is the coordinates-only bias in `g` one scale error, or does it distort the pattern?

The coverage run put the coordinates-only bias at +0.51 on a truth of 0.80 -- a factor of about
1.63 at one voxel. The whole defence of `g` without images rests on that factor being *one
number*: the PR says the magnitude is compressed but the map is still readable relatively, and a
reader comparing two regions is on safe ground. That is only true if `g/truth` is the same at
every site. If the factor varies with the true strength, the relative map is distorted too and
there is nothing left to defend.

The bed carries five sites at true weights 1.0, 1.0, 0.75, 0.75, 0.5, so the question is
answered by reading `g` at all five and asking whether `g/truth` is constant across them --
and, separately, whether the *ordering* of the five is preserved, which is the weaker claim the
docstring actually makes.

A constant ratio would mean one multiplicative constant repairs everything, which is exactly
what an image donor supplies. A ratio that shrinks as the truth grows would mean the map is
compressed non-linearly -- strong and weak regions pulled together -- and differences between
regions are understated by an amount that depends on the regions.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import ndimage, stats
from joblib import Parallel, delayed
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
import reporting

SHAPE, VOXEL_MM, SMOOTH_VOX, RADIUS_VOX = (30, 30, 30), 4.0, 0.8, 2.5
PEAK_G, N_SIMS, N_STUDIES = 0.8, 40, 12
AFF = np.eye(4); AFF[:3, :3] *= VOXEL_MM
AFF[:3, 3] = -VOXEL_MM * (np.array(SHAPE) - 1) / 2.0
MASK = nib.Nifti1Image(np.ones(SHAPE, dtype=np.int32), AFF)
MASK_BOOL = np.ones(SHAPE, dtype=bool); ZOOMS = np.full(3, VOXEL_MM)

SITE_IJK = [(10, 15, 15), (20, 15, 15), (15, 9, 18), (15, 21, 12), (15, 15, 22)]
SITE_WEIGHT = [1.0, 0.875, 0.75, 0.625, 0.5]   # spread out, so the ordering test has resolution


def truth_field():
    out = np.zeros(SHAPE)
    grid = np.indices(SHAPE).astype(float)
    for (i, j, k), w in zip(SITE_IJK, SITE_WEIGHT):
        d2 = (grid[0] - i) ** 2 + (grid[1] - j) ** 2 + (grid[2] - k) ** 2
        out = np.maximum(out, PEAK_G * w * np.exp(-d2 / (2 * RADIUS_VOX**2)))
    return out


TRUTH = truth_field()
TRUE_AT = np.array([TRUTH[s] for s in SITE_IJK])


def one(seed, n_image):
    rng = np.random.default_rng(seed)
    studies, reported = [], 0
    for k in range(N_STUDIES):
        n = int(rng.integers(20, 41))
        noise = ndimage.gaussian_filter(rng.standard_normal(SHAPE), SMOOTH_VOX)
        noise *= 1.0 / (noise.std() + 1e-12)
        gmap = TRUTH + noise / np.sqrt(n)
        foci, _ = reporting.report_peaks((gmap * np.sqrt(n))[MASK_BOOL], MASK_BOOL, SHAPE,
                                         ZOOMS, "cluster", "max")
        reported += bool(foci)
        meta = {"sample_sizes": [n]}
        analysis = {"id": f"s{k}-1", "name": "1", "metadata": meta, "points": [
            {"space": "MNI",
             "coordinates": [float(v) for v in nib.affines.apply_affine(AFF, ijk)],
             "values": [{"kind": "Z", "value": float(zv)}]} for ijk, zv in foci]}
        if k < n_image:
            tag = f"g{seed}_{n_image}_{k}"
            gp = f"/tmp/claude-0/cov_imgs/{tag}_g.nii.gz"
            vp = f"/tmp/claude-0/cov_imgs/{tag}_v.nii.gz"
            nib.save(nib.Nifti1Image(gmap.astype(np.float32), AFF), gp)
            nib.save(nib.Nifti1Image(np.full(SHAPE, 1.0 / n, np.float32), AFF), vp)
            analysis["images"] = [
                {"url": gp, "filename": "g.nii.gz", "space": "MNI", "value_type": "g"},
                {"url": vp, "filename": "v.nii.gz", "space": "MNI", "value_type": "g_var"}]
        studies.append({"id": f"s{k}", "name": f"s{k}", "metadata": meta,
                        "analyses": [analysis]})
    if reported < 2:
        return None
    est = CBES(fwhm=10.0, mask=MASK, null_method="none", use_images=n_image > 0, peak_bias=None)
    res = est.fit(Studyset({"id": "w", "name": "w", "studies": studies}, target=None, mask=MASK))
    g = res.get_map("g", return_type="array").ravel()
    return np.array([g[int(np.ravel_multi_index(s, SHAPE))] for s in SITE_IJK])


if __name__ == "__main__":
    print(f"{N_STUDIES} studies, {N_SIMS} replications, five sites")
    print("true g at the sites:", " ".join(f"{v:.3f}" for v in TRUE_AT), "\n")
    for n_image in (0, 2, 12):
        rows = [r for r in Parallel(n_jobs=8)(delayed(one)(s, n_image) for s in range(N_SIMS))
                if r is not None]
        g = np.array(rows)
        ok = np.isfinite(g).all(axis=1)
        g = g[ok]
        ratio = g / TRUE_AT
        label = "coordinates only" if n_image == 0 else f"{n_image} of {N_STUDIES} with images"
        print(f"--- {label} ({g.shape[0]} usable) ---")
        print("      true g   " + " ".join(f"{v:7.3f}" for v in TRUE_AT))
        print("      mean g   " + " ".join(f"{v:7.3f}" for v in g.mean(axis=0)))
        print("  g / true g   " + " ".join(f"{v:7.3f}" for v in ratio.mean(axis=0)))
        rho = np.array([stats.spearmanr(row, TRUE_AT)[0] for row in g])
        exact = np.mean([np.array_equal(np.argsort(row), np.argsort(TRUE_AT)) for row in g])
        r = ratio.mean(axis=0)
        print(f"  ratio spread across sites: {r.max()/r.min():.2f}x "
              f"(1.00 would mean one constant repairs the map)")
        print(f"  within-map Spearman vs truth: {np.nanmean(rho):+.3f} "
              f"(sd {np.nanstd(rho):.3f}); exact ordering in {exact:.0%} of maps\n", flush=True)
    print("A ratio spread far from 1 means the relative map is distorted, not merely rescaled,")
    print("and no single peak_bias_scale -- however well calibrated -- can put it right.")
