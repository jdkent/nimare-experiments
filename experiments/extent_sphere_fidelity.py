"""How much of a real cluster does a sphere of its reported volume actually cover?

If cluster extent is going to drive the coverage radius, the question is not whether extent
helps in principle -- the causal test answers that with the true significant set -- but how much
of that is recoverable from the one number a table prints. A paper gives ``k`` voxels, or a
volume in mm3. The estimator would turn it into a radius, ``r = (3V / 4pi)^(1/3)``, and place a
sphere of that radius at the reported focus.

Clusters are not spheres. They are elongated, they bend around sulci, and the reported focus
sits at the maximum rather than at the centre, so a sphere of the right *volume* centred at the
wrong *place* can miss a great deal. Two numbers decide whether this is worth building:

  * **recall** -- what fraction of the cluster the sphere covers. Everything missed stays
    misread as silence.
  * **precision** -- what fraction of the sphere is really in the cluster. Everything extra is
    territory where the study is wrongly excused from silence, which costs the censoring term.

Both are reported against the fixed 20 mm sphere the estimator uses now, and for a focus at the
maximum and at the centre of mass, since the two place the sphere differently.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import ndimage
from nilearn.datasets import load_mni152_brain_mask
from nilearn.image import resample_to_img
from nilearn.maskers import NiftiMasker
from nimare.transforms import ImageTransformer
from load_pain import load_pain
from reporting import cluster_extent_threshold, CLUSTER_FORMING_Z

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape = mask_img.shape
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)
voxel_volume = float(np.prod(zooms))
grid = np.stack(np.meshgrid(*[np.arange(s) for s in shape], indexing="ij"), axis=-1)

ss = ImageTransformer(target="z").transform(load_pain())
rows = []
for row, n in zip(ss.images.itertuples(), ss.sample_sizes()):
    path = getattr(row, "z", None)
    if path is None or not os.path.isfile(str(path)) or not np.isfinite(float(n)):
        continue
    img = resample_to_img(nib.load(str(path)), mask_img, interpolation="continuous",
                          force_resample=True, copy_header=True)
    rows.append(np.nan_to_num(masker.transform(img).ravel().astype(float)))

records = []
for z in rows:
    volume = np.zeros(shape); volume[mask_bool] = z
    supra = (np.abs(volume) >= CLUSTER_FORMING_Z) & mask_bool
    labels, n_labels = ndimage.label(supra)
    if not n_labels:
        continue
    sizes = np.bincount(labels.ravel())
    critical = cluster_extent_threshold(volume, mask_bool, zooms, CLUSTER_FORMING_Z)
    for cluster in range(1, n_labels + 1):
        if sizes[cluster] < critical:
            continue
        member = labels == cluster
        size = int(member.sum())
        radius = (3.0 * size * voxel_volume / (4.0 * np.pi)) ** (1.0 / 3.0)
        peak = np.array(np.unravel_index(
            np.argmax(np.where(member, np.abs(volume), -np.inf)), shape), dtype=float)
        com = np.array(ndimage.center_of_mass(member), dtype=float)
        for focus_name, centre in (("max", peak), ("com", com)):
            distance = np.linalg.norm((grid - centre) * zooms, axis=-1)
            for label, r in ((f"extent sphere ({focus_name})", radius),
                             (f"fixed 20 mm ({focus_name})", 20.0)):
                ball = distance <= r
                hit = int((ball & member).sum())
                records.append((label, size, hit / size, hit / max(int(ball.sum()), 1)))

print(f"{len(rows)} studies, {len(records) // 4} surviving clusters "
      f"at p < 0.001 with a family-wise extent test\n")
sizes_of = np.array([r[1] for r in records[::4]])
print(f"cluster size: median {np.median(sizes_of):.0f} voxels, "
      f"largest {sizes_of.max():.0f}, {(sizes_of > 200).sum()} above 200\n")
print(f"  {'placement':>26} {'recall':>8} {'precision':>10}   "
      f"{'recall, clusters > 200 vox':>27}")
for label in dict.fromkeys(r[0] for r in records):
    subset = [r for r in records if r[0] == label]
    big = [r for r in subset if r[1] > 200]
    print(f"  {label:>26} {np.mean([r[2] for r in subset]):8.1%} "
          f"{np.mean([r[3] for r in subset]):10.1%}   "
          f"{(np.mean([r[2] for r in big]) if big else float('nan')):27.1%}")
print("\nRecall is the share of the cluster that stops being misread as silence; precision is"
      "\nhow much genuinely silent territory the sphere wrongly excuses.")
