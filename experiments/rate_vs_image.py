"""Does rate calibration agree with image calibration, on data where the bias is real?

The simulator has no peak-height inflation, so it cannot test a peak-bias correction. The pain
collection can: the peaks are real local maxima of real images, and the images give an
independent value of rho to check against. Two routes to the same constant, one using images
and one using only who reported where.
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

def build(with_images):
    """Coordinate studyset; the first half also carries its images when asked."""
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
        if with_images and position < len(ids) // 2 and sid in paths:
            g_path, var_path = paths[sid]
            analysis["images"] = [
                {"url": str(g_path), "filename": "g", "space": "MNI", "value_type": "g"},
                {"url": str(var_path), "filename": "gv", "space": "MNI", "value_type": "g_var"},
            ]
        studies.append({"id": sid, "name": sid, "metadata": meta, "analyses": [analysis]})
    return Studyset({"id": "pain", "name": "pain", "studies": studies},
                    target=None, mask=masker.mask_img)

common = dict(fwhm=10.0, null_method="parametric", threshold="reporting_threshold",
              peak_bias="per-study", mask=masker, censoring="rft", smoothness_fwhm=12.16)
from_images = CBES(peak_bias_scale="images", use_images=True, **common)
from_images.fit(build(True))
from_rates = CBES(peak_bias_scale="rate-match", use_images=False, **common)
from_rates.fit(build(False))

print(f"rho from images (half the collection supplies them): "
      f"{from_images._peak_bias_scale_:.3f}")
print(f"rho, rate match + quadrature, smoothness 12.16mm:   "
      f"{from_rates._peak_bias_scale_:.3f}")
print(f"ratio rate-match/images: {from_rates._peak_bias_scale_ / from_images._peak_bias_scale_:.2f}")
print("\nfor reference, earlier independent estimates of rho on this collection:")
print("  leave-one-out pooled effect at reported peaks:  ~0.20")
print("  split-half image calibration:                   0.479 +- 0.056")
print("  means-matched ratio against the LOO truth:      0.203")
