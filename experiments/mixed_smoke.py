import sys, warnings, json; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
from scipy import stats
from load_pain import load_pain
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.meta.cbma import CBES

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
# Coordinates for every study, images also present for every study.
cs = ImagesToCoordinates(merge_strategy="demolish", z_threshold=3.2905, two_sided=True,
                         remove_subpeaks=True).transform(ss)
ref = np.load("/tmp/claude-0/cmp/reference_g.npy")
ids = list(cs.ids)

for label, keep in (("coordinates only", []), ("5 images", ids[:5]), ("all 21 images", ids)):
    est = CBES(fwhm=10.0, null_method="parametric", use_images=bool(keep))
    original = est._load_image_studies
    keep_set = set(keep)
    est._load_image_studies = lambda ds: {
        k: v for k, v in original(ds).items() if k in keep_set}
    try:
        res = est.fit(cs)
        g = res.get_map("g", return_type="array").ravel()
        cov = g != 0
        print(f"{label:18s} r={stats.pearsonr(g[cov], ref[cov])[0]:.3f} "
              f"mean|ref>.2={g[(ref>0.2)&cov].mean():.3f} (ref 0.432) coverage={cov.mean()*100:.1f}%")
    except Exception as exc:
        print(f"{label:18s} FAILED: {type(exc).__name__}: {exc}")
