"""Padded dilation must give exactly the old coverage set, and be faster.

The old implementation is inlined here so the comparison does not depend on a stash. Checked
on the real pain table, on a 3x replication with jittered coordinates, and on foci deliberately
pushed outside the image -- the case the padding argument turns on.
"""
import sys, time, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.meta.utils import _get_mask_flat_to_masked, sphere_kernel_offsets
from nimare.transforms import ImagesToCoordinates, ImageTransformer

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
ss = ImagesToCoordinates(merge_strategy="demolish", z_threshold=3.2905, two_sided=True,
                         remove_subpeaks=True).transform(ss)


def old_coverage(self, table, study_ids, active, n_voxels, image_ids=()):
    mask_img = self.masker.mask_img
    shape = np.asarray(mask_img.shape[:3], dtype=np.int64)
    mask_flat_to_masked = _get_mask_flat_to_masked(mask_img)
    radius = self.coverage_radius
    if radius is None:
        radius = 2.0 * (self.fwhm if self.fwhm is not None else 10.0)
    offsets = sphere_kernel_offsets(radius, mask_img.header.get_zooms()[:3])
    active_lookup = np.full(n_voxels, -1, dtype=np.int64)
    active_lookup[active] = np.arange(active.size)
    seen = np.zeros(active.size, dtype=bool)
    image_ids = set(image_ids)
    cols, positions = [], []
    for position, study_id in enumerate(study_ids):
        if study_id in image_ids:
            cols.append(np.arange(active.size, dtype=np.int64))
            positions.append(np.full(active.size, position, dtype=np.int64))
            continue
        ijk = table.loc[table["id"] == study_id, ["i", "j", "k"]].values.astype(np.int64)
        if not ijk.size:
            continue
        candidates = ijk[:, None, :] + offsets[None, :, :].astype(np.int64)
        in_bounds = np.all((candidates >= 0) & (candidates < shape), axis=-1)
        flat = (candidates[..., 0] * shape[1] * shape[2]
                + candidates[..., 1] * shape[2] + candidates[..., 2])
        reached = mask_flat_to_masked[np.where(in_bounds, flat, 0)]
        reached = reached[in_bounds & (reached >= 0)].astype(np.int64)
        if not reached.size:
            continue
        local = active_lookup[reached]
        local = local[local >= 0]
        if not local.size:
            continue
        seen[local] = True
        local = np.flatnonzero(seen)
        seen[local] = False
        cols.append(local)
        positions.append(np.full(local.size, position, dtype=np.int64))
    if not cols:
        empty = np.array([], dtype=np.int64)
        return empty, empty
    return np.concatenate(cols), np.concatenate(positions)


# Coordinates only: with an image for every study the focus table is empty and there is no
# coverage to compare.
from nimare.meta.cbma import effectsize as es
es.CBES._load_image_studies = lambda self, dataset: {}

estimator = CBES(fwhm=10.0, null_method="none", peak_bias="per-study")
estimator.fit(ss)
table = estimator._focus_table_
study_ids = sorted(table["id"].unique())
n_voxels = int(estimator.masker.mask_img.get_fdata().astype(bool).sum())
active = np.arange(n_voxels)

cases = {"real": table}
shape = np.asarray(estimator.masker.mask_img.shape[:3])
outside = table.copy()
# Three foci pushed out: just inside the sphere's reach, just outside it, and far away.
outside.loc[outside.index[0], ["i", "j", "k"]] = [-3, 40, 40]
outside.loc[outside.index[1], ["i", "j", "k"]] = [-40, 40, 40]
outside.loc[outside.index[2], ["i", "j", "k"]] = [shape[0] + 5, 40, 40]
cases["out of image"] = outside

for label, tbl in cases.items():
    t = time.perf_counter()
    c_old, p_old = old_coverage(estimator, tbl, study_ids, active, n_voxels)
    t_old = time.perf_counter() - t
    t = time.perf_counter()
    c_new, p_new = estimator._coverage_entries(tbl, study_ids, active, n_voxels)
    t_new = time.perf_counter() - t
    same = np.array_equal(c_old, c_new) and np.array_equal(p_old, p_new)
    print(f"{label:14s} {'identical' if same else 'DIFFERENT'}  "
          f"{c_old.size} entries   old {t_old:5.2f}s -> new {t_new:5.2f}s ({t_old/t_new:.2f}x)")
    if not same:
        print(f"    sizes {c_old.size} vs {c_new.size}")
