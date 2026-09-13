"""The same protocol on an independent real collection: 11 NeuroVault studies."""
import warnings; warnings.simplefilter("ignore")
import numpy as np
from scipy import stats
from nimare.generate import create_neurovault_studyset
from nimare.transforms import ImageTransformer, ImagesToCoordinates
from nimare.meta.cbma import CBES

ss = ImageTransformer(target=["g", "g_var"]).transform(create_neurovault_studyset())
have = [a for a, g in zip(ss.ids, ss.images["g"]) if g is not None]
ss = ss.slice(analyses=have)
cs = ImagesToCoordinates(merge_strategy="demolish", z_threshold=3.2905, two_sided=True,
                         remove_subpeaks=True).transform(ss)
ids = np.array([a for a in cs.ids if a in set(cs.coordinates["id"])])
print(f"{len(ss.ids)} studies with g images; {len(ids)} also report peaks; "
      f"{len(cs.coordinates)} peaks; N range "
      f"{np.nanmin(ss.sample_sizes()):.0f}-{np.nanmax(ss.sample_sizes()):.0f}")

def restrict(s, keep): return s.slice(analyses=[a for a in s.ids if a in set(keep)])
def fit(s, use_images, peak_bias=None):
    return CBES(fwhm=10.0, null_method="parametric", use_images=use_images,
                peak_bias=peak_bias).fit(s).get_map("g", return_type="array").ravel()

rng = np.random.default_rng(0)
rhos, before, after = [], [], []
for rep in range(4):
    perm = rng.permutation(len(ids))
    calib, test = ids[perm[: len(ids) // 2]], ids[perm[len(ids) // 2:]]
    cal = restrict(cs, calib)
    im, co = fit(cal, True), fit(cal, False)
    ok = (im != 0) & (co != 0)
    rho = float(im[ok].mean() / co[ok].mean())
    rhos.append(rho)

    ts = restrict(cs, test)
    truth, raw, fixed = fit(ts, True), fit(ts, False), fit(ts, False, peak_bias=rho)
    m = (truth != 0) & (raw != 0)
    before.append(raw[m].mean() / truth[m].mean())
    after.append(fixed[m].mean() / truth[m].mean())
    print(f"  rep {rep}: rho={rho:.3f}  raw/image={before[-1]:.2f}x  "
          f"corrected/image={after[-1]:.2f}x", flush=True)

print(f"\nNeuroVault rho = {np.mean(rhos):.3f} +- {np.std(rhos):.3f}   (pain collection: 0.479)")
print(f"coordinate/image ratio  before {np.mean(before):.2f}x   after {np.mean(after):.2f}x")
