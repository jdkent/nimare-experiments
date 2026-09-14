"""Where does a whole-brain coordinate-only fit actually spend its time?

The censoring kernel measures 21% of the fit, not the 63% the docstring claims, so the
remaining 79% is unaccounted for. cProfile, cumulative time, on the 63-study replication --
big enough that per-call overhead is not what is being measured.
"""
import cProfile, pstats, sys, warnings; warnings.simplefilter("ignore")
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
        n = sizes.get(sid, 20.0)
        points = [{"space": "MNI",
                   "coordinates": [float(r.x + offset[i][0]), float(r.y + offset[i][1]),
                                   float(r.z + offset[i][2])],
                   "values": [{"kind": "Z", "value": float(r.z_stat)}]}
                  for i, r in zip(rows.index, rows.itertuples())]
        studies.append({"id": name, "name": name, "metadata": {"sample_sizes": [n]},
                        "analyses": [{"id": f"{name}-1", "name": "1",
                                      "metadata": {"sample_sizes": [n]},
                                      "points": points, "images": []}]})
dataset = Studyset({"id": "big", "name": "big", "studies": studies})

estimator = CBES(fwhm=10.0, mask=mask, null_method="none", peak_bias="per-study")
estimator.fit(dataset)  # warm the geometry cache and any import-time work
profiler = cProfile.Profile()
profiler.enable()
CBES(fwhm=10.0, mask=mask, null_method="none", peak_bias="per-study").fit(dataset)
profiler.disable()
stats = pstats.Stats(profiler)
stats.sort_stats("tottime").print_stats(25)
