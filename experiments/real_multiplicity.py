"""How bad is the relocation null's density mismatch on a real collection?

The toy that exposed it used a cubic mask with foci confined to the middle 58% of it, which is
an unrealistically severe mismatch. What decides whether this needs a code change is the size
of the effect on data people will actually run: a real brain mask, and foci reported by real
studies, which concentrate in gray matter and in whatever the literature was looking at.

The driver is one number -- mean studies contributing per covered voxel -- measured on the
observed foci and on relocations of them. If observed and relocated agree on a real collection,
this is a documentation matter and a fit-time warning; if they do not, the relocation domain
itself has to change.

Three domains are compared, as in relocation_domain.py: the uniform relocation the estimator
does now, the dilated bounding box of the observed foci, and a smoothed foci density.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import ndimage
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.transforms import ImageTransformer, ImagesToCoordinates

ITERS = int(sys.argv[1]) if len(sys.argv) > 1 else 15
U = 3.2905  # one-tailed p < .001, the modal reporting threshold


def domain(est, kind):
    ijk_all = est._in_mask_ijk()
    foci = est._focus_table_[["i", "j", "k"]].values.astype(int)
    if kind == "uniform":
        return ijk_all, None
    if kind == "hull":
        # the dilated *union* of the foci, not their bounding box: on a real brain the box is
        # the whole mask and the restriction does nothing
        shape = np.asarray(est.masker.mask_img.shape[:3])
        occupied = np.zeros(shape, dtype=bool)
        occupied[tuple(foci.T)] = True
        radius = (est.coverage_radius or 2 * est.fwhm)
        zooms = np.array(est.masker.mask_img.header.get_zooms()[:3], dtype=float)
        grid = np.stack(np.meshgrid(*[np.arange(-int(np.ceil(radius / z)),
                                                int(np.ceil(radius / z)) + 1) for z in zooms],
                                    indexing="ij"), axis=-1)
        ball = (((grid * zooms) ** 2).sum(-1) <= radius ** 2)
        occupied = ndimage.binary_dilation(occupied, structure=ball)
        keep = occupied[tuple(ijk_all.T)]
        return ijk_all[keep], None
    if kind == "density":
        shape = np.asarray(est.masker.mask_img.shape[:3])
        vol = np.zeros(shape, dtype=float)
        np.add.at(vol, tuple(foci.T), 1.0)
        sigma = 3.0 * est.fwhm / np.array(est.masker.mask_img.header.get_zooms()[:3])
        vol = ndimage.gaussian_filter(vol, sigma=sigma)
        p = vol[tuple(ijk_all.T)]
        p = p + p.max() * 1e-3
        return ijk_all, p / p.sum()
    raise ValueError(kind)


def multiplicity(est, table):
    fit, z = est._statistic(table, est._sample_sizes_, est._thresholds_, est._image_studies_)
    c = fit["covered"]
    return float(fit["n_studies"][c].mean()), float(c.mean()), float(np.abs(z[c]).mean())


ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
ss = ImagesToCoordinates(merge_strategy="demolish", z_threshold=U, two_sided=True,
                         remove_subpeaks=True).transform(ss)
print(f"NIDM pain peaks (z > {U:.4f} out of the 21 full images), real MNI152 mask, "
      f"{ITERS} relocations.\n", flush=True)
print(f"{'source':>10s} {'n_studies':>10s} {'covered':>8s} {'mean |z|':>9s} {'ratio':>7s}",
      flush=True)

est = CBES(fwhm=10.0, null_method="none", threshold="study-min", peak_bias="per-study",
           use_images=False, seed=0)
est.fit(ss)
n_obs, c_obs, z_obs = multiplicity(est, est._focus_table_)
print(f"{'observed':>10s} {n_obs:10.2f} {c_obs:8.3f} {z_obs:9.3f} {'':>7s}", flush=True)

for kind in ("uniform", "hull", "density", "values"):
    rng = np.random.default_rng(7)
    ijk, p = (None, None) if kind == "values" else domain(est, kind)
    rows = []
    for _ in range(ITERS):
        t = est._focus_table_.copy()
        if kind == "values":
            # positions and study membership stay put; only the magnitudes move, so the
            # spatial design -- coverage and studies per voxel -- is preserved exactly
            order = rng.permutation(len(t))
            t[["g", "var_g"]] = t[["g", "var_g"]].values[order]
        else:
            idx = rng.choice(len(ijk), size=len(t), p=p) if p is not None \
                else rng.integers(0, len(ijk), size=len(t))
            t[["i", "j", "k"]] = ijk[idx]
        rows.append(multiplicity(est, t))
    n, c, z = np.mean(np.array(rows), axis=0)
    print(f"{kind:>10s} {n:10.2f} {c:8.3f} {z:9.3f} {n_obs / n:7.2f}", flush=True)

print(f"\nmask voxels {int(est._mask_bool().sum())}, foci {len(est._focus_table_)}", flush=True)
