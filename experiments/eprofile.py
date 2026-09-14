"""Line-level split of the two hot EM kernels, so the next optimization is aimed, not guessed."""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from line_profiler import LineProfiler
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
        n = sizes.get(sid, 20.0)
        studies.append({"id": name, "name": name, "metadata": {"sample_sizes": [n]},
                        "analyses": [{"id": f"{name}-1", "name": "1",
                                      "metadata": {"sample_sizes": [n]},
                                      "points": points, "images": []}]})
dataset = Studyset({"id": "big", "name": "big", "studies": studies})

lp = LineProfiler()
lp.add_function(es.CBES._update_prevalence)
lp.add_function(es._censoring_terms)
lp.add_function(es.CBES._coverage_entries)
lp.enable_by_count()
CBES(fwhm=10.0, mask=mask, null_method="none", peak_bias="per-study").fit(dataset)
lp.disable_by_count()
lp.print_stats()
