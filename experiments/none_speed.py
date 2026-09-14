"""How much faster is selection_model='none' really? The docstring claims about 3x."""
import sys, time, warnings; warnings.simplefilter("ignore")
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
coords = ss.coordinates
sizes = {r.study_id: float(np.mean(r.sample_sizes))
         for r in ss.metadata.itertuples() if r.sample_sizes}
mask = ss.masker.mask_img
es.CBES._load_image_studies = lambda self, dataset: {}
rng = np.random.default_rng(0)
studies = []
for c in range(3):
    jitter = rng.normal(0.0, 4.0, (len(coords), 3)) if c else np.zeros((len(coords), 3))
    offset = dict(zip(coords.index, jitter))
    for sid, rows in coords.groupby("study_id"):
        name = f"{sid}_r{c}"
        points = [{"space": "MNI",
                   "coordinates": [float(r.x + offset[i][0]), float(r.y + offset[i][1]),
                                   float(r.z + offset[i][2])],
                   "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                  for i, r in zip(rows.index, rows.itertuples())]
        studies.append({"id": name, "name": name, "metadata": {"sample_sizes": [sizes.get(sid, 20.0)]},
                        "analyses": [{"id": f"{name}-1", "name": "1",
                                      "metadata": {"sample_sizes": [sizes.get(sid, 20.0)]},
                                      "points": points, "images": []}]})
dataset = Studyset({"id": "big", "name": "big", "studies": studies})

times = {"zero-inflated": [], "none": []}
for _ in range(3):
    for model in ("zero-inflated", "none"):
        t = time.perf_counter()
        CBES(fwhm=10.0, mask=mask, null_method="none", peak_bias="per-study",
             selection_model=model).fit(dataset)
        times[model].append(time.perf_counter() - t)
best = {k: min(v) for k, v in times.items()}
print(f"zero-inflated {best['zero-inflated']:6.2f}s   none {best['none']:6.2f}s   "
      f"{best['zero-inflated'] / best['none']:.2f}x faster")
print("  all runs:", {k: [round(x, 2) for x in v] for k, v in times.items()})
