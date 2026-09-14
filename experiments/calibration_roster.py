"""Does a calibration fit that receives the whole roster get a different scale?

``_calibrate_peak_bias_scale`` fits the coordinates alone and then each image donor alone,
but hands every one of those fits the full ``sample_sizes``/``thresholds`` roster. A study
with no foci in the table it was given falls through ``_coverage_entries`` as "reported
nothing anywhere: silent at every voxel", so the donor-only fit carries one censored-silent
observation per non-donor study -- and the coordinate-only fit carries one per image donor.

Measures the per-donor ratios and the pooled scale with the roster as-is against the roster
restricted to the studies each fit actually contains.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
import pandas as pd
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import ImageTransformer, ImagesToCoordinates

U = 3.2905
ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
masker = ss.masker
sizes = dict(zip([str(i) for i in ss.ids], ss.sample_sizes()))
coords = ImagesToCoordinates(
    merge_strategy="demolish", z_threshold=U, two_sided=True, remove_subpeaks=True
).transform(ss).coordinates
paths = {str(r.id): (r.g, r.g_var) for r in ss.images.itertuples() if r.g is not None}
ids = sorted(set(str(i) for i in coords["id"]))
N_DONORS = 3

studies = []
for position, sid in enumerate(ids):
    sub = coords[coords["id"].astype(str) == sid]
    meta = {"sample_sizes": [int(sizes[sid])], "reporting_threshold": U}
    analysis = {
        "id": sid, "name": "1", "metadata": meta,
        "points": [{"space": "MNI", "coordinates": [float(r.x), float(r.y), float(r.z)],
                    "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                   for r in sub.itertuples()],
    }
    if position < N_DONORS and sid in paths:
        g_path, var_path = paths[sid]
        analysis["images"] = [
            {"url": str(g_path), "filename": "g", "space": "MNI", "value_type": "g"},
            {"url": str(var_path), "filename": "gv", "space": "MNI", "value_type": "g_var"},
        ]
    studies.append({"id": sid, "name": sid, "metadata": meta, "analyses": [analysis]})
studyset = Studyset({"id": "pain", "name": "pain", "studies": studies},
                    target=None, mask=masker.mask_img)

est = CBES(fwhm=10.0, null_method="none", threshold="reporting_threshold",
           peak_bias="per-study", peak_bias_scale="images", use_images=True, mask=masker)
est.fit(studyset)
print(f"roster: {len(est._sample_sizes_)} studies, {len(est._image_studies_)} image donors, "
      f"{len(est._focus_table_)} foci")

table, sample_sizes, thresholds = est._focus_table_, est._sample_sizes_, est._thresholds_
image_studies = est._image_studies_
with_foci = set(str(v) for v in table["id"].unique())


def one_ratio(coordinate_only, single):
    both = (coordinate_only["covered"] & single["covered"]
            & np.isfinite(coordinate_only["g"]) & np.isfinite(single["g"]))
    if not both.any():
        return np.nan
    from_coordinates, from_image = np.abs(coordinate_only["g"][both]), np.abs(single["g"][both])
    strong = from_image >= np.percentile(from_image, 75)
    if strong.sum() < 50:
        strong = np.ones_like(from_image, dtype=bool)
    ratios = from_image[strong] / np.clip(from_coordinates[strong], 1e-12, None)
    ratios = ratios[np.isfinite(ratios) & (ratios > 0)]
    return float(np.median(ratios)) if ratios.size else np.nan


def run(restrict):
    """Per-donor ratios, with the roster restricted to each fit's own studies or not."""
    est._coverage_ = None
    coord_ids = sorted(with_foci - set(image_studies)) if restrict else list(sample_sizes.index)
    coordinate_only = est._statistic(
        table[table["id"].astype(str).isin(coord_ids)] if restrict else table,
        sample_sizes.loc[coord_ids], thresholds.loc[coord_ids], image_studies=None)[0]
    est._coverage_ = None
    out = {}
    for study_id, payload in image_studies.items():
        rows = [study_id] if restrict else list(sample_sizes.index)
        single = est._statistic(table.iloc[:0], sample_sizes.loc[rows], thresholds.loc[rows],
                                image_studies={study_id: payload})[0]
        est._coverage_ = None
        out[study_id] = (one_ratio(coordinate_only, single), float(np.nanmedian(
            np.abs(single["g"][single["covered"]]))))
    return coordinate_only, out


whole_coord, whole = run(restrict=False)
own_coord, own = run(restrict=True)

print(f"\ncoordinate-only fit, median |g| over covered voxels:")
print(f"  whole roster (image donors read as silent everywhere): "
      f"{np.nanmedian(np.abs(whole_coord['g'][whole_coord['covered']])):.4f}")
print(f"  only the studies with foci:                            "
      f"{np.nanmedian(np.abs(own_coord['g'][own_coord['covered']])):.4f}")
print(f"\n{'donor':>8} {'|g| whole':>11} {'|g| own':>9} {'ratio whole':>13} {'ratio own':>11}")
for study_id in image_studies:
    rw, gw = whole[study_id]
    ro, go = own[study_id]
    print(f"{study_id:>8} {gw:>11.4f} {go:>9.4f} {rw:>13.4f} {ro:>11.4f}")
rw = [v[0] for v in whole.values() if np.isfinite(v[0])]
ro = [v[0] for v in own.values() if np.isfinite(v[0])]
print(f"\npooled scale, whole roster:      {np.median(rw):.4f}  "
      f"(interval {min(rw):.4f}-{max(rw):.4f})")
print(f"pooled scale, own studies only:  {np.median(ro):.4f}  "
      f"(interval {min(ro):.4f}-{max(ro):.4f})")
print(f"packaged estimator reported:     {est._peak_bias_scale_:.4f}")
