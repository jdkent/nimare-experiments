"""Is the NeuroVault failure too few reporting studies, or something deeper?"""
import warnings; warnings.simplefilter("ignore")
import numpy as np
from nimare.generate import create_neurovault_studyset
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.meta.cbma import CBES

base = ImageTransformer(target=["g", "g_var"]).transform(create_neurovault_studyset())
base = base.slice(analyses=[a for a, g in zip(base.ids, base.images["g"]) if g is not None])

def fit(s, use_images, peak_bias=None):
    return CBES(fwhm=10.0, null_method="parametric", use_images=use_images,
                peak_bias=peak_bias).fit(s).get_map("g", return_type="array").ravel()

print(f"{'threshold z':>11s} {'studies w/ peaks':>17s} {'peaks':>7s} {'raw/image':>10s} "
      f"{'rho':>6s} {'corrected/image':>16s}")
for zt in (3.2905, 2.5758, 2.3263, 1.96):
    cs = ImagesToCoordinates(merge_strategy="demolish", z_threshold=zt, two_sided=True,
                             remove_subpeaks=True).transform(base)
    reporting = sorted(set(cs.coordinates["id"]))
    if len(reporting) < 3:
        print(f"{zt:11.3f} {len(reporting):17d} {len(cs.coordinates):7d}   (too few)")
        continue
    sub = cs.slice(analyses=reporting)
    im, co = fit(sub, True), fit(sub, False)
    ok = (im != 0) & (co != 0)
    rho = float(im[ok].mean() / co[ok].mean())
    fixed = fit(sub, False, peak_bias=rho)
    ok2 = (im != 0) & (fixed != 0)
    print(f"{zt:11.3f} {len(reporting):17d} {len(cs.coordinates):7d} "
          f"{co[ok].mean()/im[ok].mean():10.2f}x {rho:6.3f} "
          f"{fixed[ok2].mean()/im[ok2].mean():15.2f}x", flush=True)
