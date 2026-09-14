"""How common are unestimable voxels at realistic peak counts?

The full pain tables carry ~130 peaks per study, which no paper reports. A typical coordinate
table has 3-20, and a thinner table means fewer studies reach any given voxel -- which is
exactly the regime where a per-voxel mask on g would bite. Peaks kept strongest-first, as a
paper reporting only its top peaks gives.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.meta.cbma import effectsize as es
from nimare.studyset import Studyset
from nimare.transforms import ImagesToCoordinates, ImageTransformer

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
ss = ImagesToCoordinates(merge_strategy="demolish", z_threshold=3.2905, two_sided=True,
                         remove_subpeaks=True).transform(ss)
coords = ss.coordinates.copy()
coords["_abs"] = coords["z_stat"].abs()
sizes = {r.study_id: float(np.mean(r.sample_sizes))
         for r in ss.metadata.itertuples() if r.sample_sizes}
mask = ss.masker.mask_img
es.CBES._load_image_studies = lambda self, dataset: {}


def build(cap):
    keep = coords.sort_values("_abs", ascending=False).groupby("study_id", sort=False).head(cap)
    studies = []
    for sid, rows in keep.groupby("study_id"):
        n = sizes.get(sid, 20.0)
        points = [{"space": "MNI", "coordinates": [float(r.x), float(r.y), float(r.z)],
                   "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                  for r in rows.itertuples()]
        studies.append({"id": sid, "name": sid, "metadata": {"sample_sizes": [n]},
                        "analyses": [{"id": f"{sid}-1", "name": "1",
                                      "metadata": {"sample_sizes": [n]},
                                      "points": points, "images": []}]})
    return Studyset({"id": "pain", "name": "pain", "studies": studies})


n_brain = int(np.asarray(mask.dataobj).astype(bool).sum())
print(f"{'peaks/study':>11s} {'covered':>8s} {'k=1':>7s} {'k<=2':>7s} {'k>=5':>7s} "
      f"{'median k':>9s} {'|z|>2 at k=1':>13s}")
for cap in (3, 5, 10, 20, 1000):
    result = CBES(fwhm=10.0, mask=mask, null_method="none", peak_bias="per-study").fit(build(cap))
    z = result.get_map("z", return_type="array").ravel()
    k = result.get_map("n_studies", return_type="array").ravel()
    covered = k > 0
    c = covered.sum()
    label = "all" if cap == 1000 else str(cap)
    single = (k == 1) & covered
    print(f"{label:>11s} {100 * c / n_brain:7.1f}% {100 * single.sum() / c:6.1f}% "
          f"{100 * ((k <= 2) & covered).sum() / c:6.1f}% {100 * (k >= 5).sum() / c:6.1f}% "
          f"{np.median(k[covered]):9.0f} "
          f"{100 * np.mean(np.abs(z[single]) > 2) if single.any() else float('nan'):12.1f}%")
