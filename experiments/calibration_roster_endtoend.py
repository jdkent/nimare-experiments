"""The calibrated scale the packaged estimator reports, with and without the roster fix.

The replication in calibration_roster.py measures the mechanism but not the shipped number:
it reimplements the ratio and came out at 0.72 where the estimator said 0.61. This asks the
estimator itself, monkeypatching the roster restriction off to get the old behaviour, so the
two numbers are produced by the same code path and differ only in what each fit receives.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.meta.cbma import effectsize as module
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


def build(n_donors):
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
        if position < n_donors and sid in paths:
            g_path, var_path = paths[sid]
            analysis["images"] = [
                {"url": str(g_path), "filename": "g", "space": "MNI", "value_type": "g"},
                {"url": str(var_path), "filename": "gv", "space": "MNI", "value_type": "g_var"},
            ]
        studies.append({"id": sid, "name": sid, "metadata": meta, "analyses": [analysis]})
    return Studyset({"id": "pain", "name": "pain", "studies": studies},
                    target=None, mask=masker.mask_img)


fixed = module.CBES._calibrate_peak_bias_scale


def whole_roster(self, table, sample_sizes, thresholds, image_studies):
    """The old behaviour: every fit receives the full roster, as before the fix."""
    real_statistic = module.CBES._statistic

    def widened(inner, sub_table, sub_sizes, sub_thresholds, image_studies=None):
        return real_statistic(inner, sub_table, sample_sizes, thresholds, image_studies)

    module.CBES._statistic = widened
    try:
        # The coordinate fit also got the donors' foci, which is what the unrestricted table was.
        return fixed(self, table, sample_sizes, thresholds, image_studies)
    finally:
        module.CBES._statistic = real_statistic


print(f"{'donors':>7s} {'whole roster':>14s} {'own studies':>13s} {'factor':>8s}")
for n_donors in (2, 3, 5):
    studyset = build(n_donors)
    common = dict(fwhm=10.0, null_method="none", threshold="reporting_threshold",
                  peak_bias="per-study", peak_bias_scale="images", use_images=True, mask=masker)
    module.CBES._calibrate_peak_bias_scale = whole_roster
    try:
        old = CBES(**common)
        old.fit(studyset)
    finally:
        module.CBES._calibrate_peak_bias_scale = fixed
    new = CBES(**common)
    new.fit(studyset)
    a, b = old._peak_bias_scale_, new._peak_bias_scale_
    print(f"{n_donors:7d} {a:14.4f} {b:13.4f} {b / a:8.3f}", flush=True)
