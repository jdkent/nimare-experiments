"""Does a flat single-study likelihood protect the z map on its own?

g at a one-reporting-study voxel sits somewhere on a plateau, so it is arbitrary. But se under
the selection model comes from the curvature of that same likelihood, and a plateau has almost
none -- so se should blow up and z = g/se should collapse toward zero without anyone masking
anything. If that holds, the z map is not merely permutation-valid at those voxels, it is
actively down-weighting them.
"""
import sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/home/user/nimare-experiments/experiments")
import numpy as np
from load_pain import load_pain
from nimare.meta.cbma import CBES
from nimare.meta.cbma import effectsize as es
from nimare.transforms import ImagesToCoordinates, ImageTransformer

ss = ImageTransformer(target=["g", "g_var"]).transform(load_pain())
ss = ImagesToCoordinates(merge_strategy="demolish", z_threshold=3.2905, two_sided=True,
                         remove_subpeaks=True).transform(ss)
es.CBES._load_image_studies = lambda self, dataset: {}  # coordinates only

result = CBES(fwhm=10.0, null_method="none", peak_bias="per-study").fit(ss)
g = result.get_map("g", return_type="array").ravel()
se = result.get_map("se", return_type="array").ravel()
z = result.get_map("z", return_type="array").ravel()
k = result.get_map("n_studies", return_type="array").ravel()

covered = k > 0
print(f"{'n_studies':>10s} {'voxels':>8s} {'% brain':>8s} {'median |g|':>11s} "
      f"{'median se':>10s} {'median |z|':>11s} {'% |z|>2':>8s}")
for label, sel in [("1", k == 1), ("2", k == 2), ("3-4", (k >= 3) & (k <= 4)),
                   ("5+", k >= 5)]:
    sel = sel & covered
    if not sel.any():
        continue
    print(f"{label:>10s} {sel.sum():8d} {100 * sel.sum() / covered.sum():7.1f}% "
          f"{np.median(np.abs(g[sel])):11.3f} {np.median(se[sel]):10.3f} "
          f"{np.median(np.abs(z[sel])):11.3f} {100 * np.mean(np.abs(z[sel]) > 2):7.1f}%")

single = (k == 1) & covered
print(f"\nsingle-study voxels are {100 * single.sum() / covered.sum():.1f}% of the covered brain")
print(f"they hold {100 * np.abs(z[single]).sum() / np.abs(z[covered]).sum():.1f}% of the total |z|")
print(f"and {100 * np.mean(np.abs(z[single]) > 2):.2f}% of them exceed |z| = 2, "
      f"against {100 * np.mean(np.abs(z[k >= 5]) > 2):.1f}% at five or more studies")
