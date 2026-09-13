"""With 10 images in hand, do extra coordinate studies add anything?

Truth is the 21-study image pooling. Compare an images-only meta-analysis against the same
images plus coordinate studies, and against simply having had more images -- which gives an
exchange rate between the two kinds of study.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from scipy import stats
from load_pain import load_pain
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.meta.cbma import CBES

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
cs = ImagesToCoordinates(merge_strategy="demolish", z_threshold=3.2905, two_sided=True,
                         remove_subpeaks=True).transform(ss)
ref = np.load("/tmp/claude-0/cmp/reference_g.npy")
ids = np.array(list(cs.ids))

def restrict(s, keep): return s.slice(analyses=[a for a in s.ids if a in set(keep)])

def run(image_ids, coord_ids, peak_bias=None):
    keep = set(image_ids) | set(coord_ids)
    sub = restrict(cs, keep)
    est = CBES(fwhm=10.0, null_method="parametric", use_images=bool(len(image_ids)),
               peak_bias=peak_bias)
    original = est._load_image_studies
    est._load_image_studies = lambda d: {a: b for a, b in original(d).items()
                                         if a in set(image_ids)}
    g = est.fit(sub).get_map("g", return_type="array").ravel()
    cov = g != 0
    return (stats.pearsonr(g[cov], ref[cov])[0], g[cov & (ref > 0.2)].mean())

rng = np.random.default_rng(0)
REPS = 5
rows = {}
for rep in range(REPS):
    order = rng.permutation(ids)
    imgs10, rest11 = order[:10], order[10:]
    trials = {
        "10 images only": (imgs10, []),
        "10 images + 11 coords (raw)": (imgs10, rest11),
        "10 images + 11 coords (peak_bias)": (imgs10, rest11),
        "15 images only": (order[:15], []),
        "21 images only": (order, []),
        "0 images, 21 coords (peak_bias)": ([], order),
    }
    for label, (im, co) in trials.items():
        pb = 0.479 if "peak_bias" in label else None
        r, m = run(im, co, pb)
        rows.setdefault(label, []).append((r, m))
    print(f"  rep {rep} done", flush=True)

print(f"\n{'configuration':36s} {'r with truth':>13s} {'mean g':>8s}   (truth {ref[ref>0.2].mean():.3f})")
for label in ["10 images only", "10 images + 11 coords (raw)",
              "10 images + 11 coords (peak_bias)", "15 images only", "21 images only",
              "0 images, 21 coords (peak_bias)"]:
    v = np.array(rows[label])
    print(f"{label:36s} {v[:,0].mean():13.3f} {v[:,1].mean():8.3f}")
