"""With two images in hand, *where* do the coordinate studies help and where do they hurt?

On real data the whole-brain answer is already in: coordinates lose. Holding the image count
fixed and adding the rest of the half as thresholded tables, scored against a held-out half of
the NIDM pain collection, images alone win on magnitude *and* on localisation at every image
count -- at two images, r +0.723 against +0.545 and AUC 0.905 against 0.868.

But a whole-brain average cannot answer the design question, which is not "are coordinates
better" but "is there a circumstance in which they earn their place". A pooled number hides a
trade: if the tables help where two images know almost nothing and hurt where the images are
already strong, that is a real finding and a reason to keep the pathway. If they hurt everywhere
the pathway is dead weight whenever two images exist.

So the same comparison, scored per voxel and stratified two ways.

**By the truth's own magnitude.** Deciles of the held-out reference. Coordinates are known to be
compressed, so the prior expectation is that they hurt most where the truth is large -- but the
interesting cell is the bottom, where two images are mostly noise.

**By how many coordinate studies actually reported near the voxel.** This is the circumstance that
matters, because a coordinate study contributes nothing at a voxel no focus of its reaches. The
mixed fit's ``n_studies`` less the two images gives the count directly. Where it is zero the two
arms should agree exactly, which doubles as a check that the stratification is real rather than
an artefact of the fit.

Design, unchanged from `do_coordinates_help.py` because it is the one that survived review:
truth is the inverse-variance pooling of the *other* half, so no study is on both sides; the
image arm and the coordinate arm are *different* studies, so this is not degrading an available
image to a table; coordinates come from `reporting` -- multiplicity-corrected, whole surviving
clusters, one focus each, **nothing capped**; the real per-study threshold is passed through
metadata rather than inferred; and the mixed fit uses the documented configuration,
``peak_bias="per-study"`` with ``peak_bias_scale="images"``.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn.datasets import load_mni152_brain_mask
from nilearn.maskers import NiftiMasker
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset
from nimare.transforms import d_to_g, t_to_d
from load_pain import load_pain
from reporting import report_peaks

SCHEME, FOCUS = "cluster", "max"
N_IMAGES = int(os.environ.get("NIMAGES", 2))
N_SPLITS = int(os.environ.get("NSPLITS", 8))
WORKDIR = f"/tmp/claude-0/wdch_{os.getpid()}"
os.makedirs(WORKDIR, exist_ok=True)

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)


def to_g(z, n):
    t = np.sign(z) * np.abs(stats.t.isf(stats.norm.sf(np.abs(z)), n - 1))
    return d_to_g(t_to_d(np.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0), n), n)


def g_and_var(z, n):
    g = to_g(z, n)
    return g, 1.0 / n + g**2 / (2.0 * n)


def pooled(members, maps, sizes):
    stack, weights = [], []
    for i in members:
        g, var = g_and_var(maps[i], sizes[i])
        stack.append(g)
        weights.append(1.0 / np.maximum(var, 1e-9))
    stack, weights = np.array(stack), np.array(weights)
    return np.sum(stack * weights, axis=0) / np.maximum(weights.sum(axis=0), 1e-12)


def write_image(values, path):
    volume = np.zeros(shape, dtype=np.float32)
    volume[mask_bool] = values
    nib.save(nib.Nifti1Image(volume, affine), path)
    return path


def mixed_fit(image_members, coord_members, maps, sizes):
    studies = []
    for i in image_members:
        g, var = g_and_var(maps[i], sizes[i])
        meta = {"sample_sizes": [int(sizes[i])]}
        studies.append({"id": f"i{i}", "name": f"i{i}", "metadata": meta, "analyses": [
            {"id": f"i{i}", "name": "1", "metadata": meta, "points": [], "images": [
                {"url": write_image(g, f"{WORKDIR}/g_{i}.nii.gz"),
                 "filename": f"g_{i}.nii.gz", "value_type": "g", "space": "MNI"},
                {"url": write_image(var, f"{WORKDIR}/gvar_{i}.nii.gz"),
                 "filename": f"gvar_{i}.nii.gz", "value_type": "g_var", "space": "MNI"}]}]})
    for i in coord_members:
        foci, height = report_peaks(maps[i], mask_bool, shape, zooms,
                                    scheme=SCHEME, focus=FOCUS)
        if not foci:
            continue
        meta = {"sample_sizes": [int(sizes[i])], "reporting_threshold": float(height)}
        studies.append({"id": f"c{i}", "name": f"c{i}", "metadata": meta, "analyses": [
            {"id": f"c{i}", "name": "1", "metadata": meta, "points": [
                {"space": "MNI",
                 "coordinates": [float(c) for c in nib.affines.apply_affine(
                     affine, np.asarray(ijk, dtype=float))],
                 "values": [{"kind": "Z", "value": v}]} for ijk, v in foci]}]})
    if len(studies) < 2:
        return None
    est = CBES(fwhm=10.0, mask=masker, peak_bias="per-study", peak_bias_scale="images",
               null_method="none", threshold="reporting_threshold")
    res = est.fit(Studyset({"id": "x", "name": "x", "studies": studies},
                           target=None, mask=mask_img))
    return (np.abs(res.get_map("g", return_type="array").ravel()),
            res.get_map("n_studies", return_type="array").ravel())


ss = load_pain()
maps, sizes = [], []
for row, n in zip(ss.images.itertuples(), ss.sample_sizes()):
    z = masker.transform(row.z).ravel()
    maps.append(np.nan_to_num(z))
    sizes.append(int(n if np.ndim(n) == 0 else np.asarray(n).ravel()[0]))
maps = np.array(maps)
n_studies_total = len(maps)
print(f"NIDM pain: {n_studies_total} studies, {SCHEME}/{FOCUS} reporting, "
      f"{N_IMAGES} images held fixed, {N_SPLITS} splits\n")

rng = np.random.default_rng(0)
by_truth = {d: {"img": [], "mix": [], "img_signed": [], "mix_signed": [], "n": 0}
            for d in range(10)}
by_count = {c: {"img": [], "mix": [], "img_signed": [], "mix_signed": [], "n": 0}
            for c in (0, 1, 2, 3)}

for split in range(N_SPLITS):
    order = rng.permutation(n_studies_total)
    half = n_studies_total // 2
    work, held = order[:half], order[half:]
    truth = np.abs(pooled(held, maps, sizes))
    image_members, coord_members = work[:N_IMAGES], work[N_IMAGES:]
    images_only = np.abs(pooled(image_members, maps, sizes))
    out = mixed_fit(image_members, coord_members, maps, sizes)
    if out is None:
        continue
    mixed, n_reaching = out
    # Coordinate studies reaching each voxel: the images cover everywhere, so subtract them.
    n_coord = np.maximum(n_reaching - len(image_members), 0)
    use = np.isfinite(truth) & np.isfinite(images_only) & np.isfinite(mixed)
    deciles = np.digitize(truth[use], np.quantile(truth[use], np.linspace(0.1, 0.9, 9)))
    for d in range(10):
        sel = deciles == d
        if not sel.any():
            continue
        by_truth[d]["img"].append(np.abs(images_only[use][sel] - truth[use][sel]).mean())
        by_truth[d]["mix"].append(np.abs(mixed[use][sel] - truth[use][sel]).mean())
        # Signed error too, because the metric above can be gamed by shrinkage. Both arms and
        # the truth are absolute values of a pooled effect, and |x| is upward-biased when x is
        # near zero -- far more so from two studies than from the ten in the held-out half. So
        # anything that pulls the estimate toward zero looks better at low truth for a reason
        # that has nothing to do with information. The signed columns expose that.
        by_truth[d]["img_signed"].append((images_only[use][sel] - truth[use][sel]).mean())
        by_truth[d]["mix_signed"].append((mixed[use][sel] - truth[use][sel]).mean())
        by_truth[d]["n"] += int(sel.sum())
    counts = np.clip(n_coord[use], 0, 3)
    for c in (0, 1, 2, 3):
        sel = counts == c
        if not sel.any():
            continue
        by_count[c]["img"].append(np.abs(images_only[use][sel] - truth[use][sel]).mean())
        by_count[c]["mix"].append(np.abs(mixed[use][sel] - truth[use][sel]).mean())
        by_count[c]["img_signed"].append((images_only[use][sel] - truth[use][sel]).mean())
        by_count[c]["mix_signed"].append((mixed[use][sel] - truth[use][sel]).mean())
        by_count[c]["n"] += int(sel.sum())
    print(f"  split {split + 1}/{N_SPLITS} done", flush=True)


def report(title, table, key_name, signed=False):
    print(f"\n--- {title} ---")
    head = (f"{key_name:>22s} {'voxels':>9s} {'|err| images':>13s} {'|err| mixed':>12s} "
            f"{'mixed - images':>15s}")
    if signed:
        head += f" {'signed images':>14s} {'signed mixed':>13s}"
    print(head)
    for key, cell in table.items():
        if not cell["img"]:
            continue
        a, b = float(np.mean(cell["img"])), float(np.mean(cell["mix"]))
        line = f"{str(key):>22s} {cell['n']:9d} {a:13.3f} {b:12.3f} {b - a:+15.3f}"
        if signed:
            line += (f" {float(np.mean(cell['img_signed'])):+14.3f}"
                     f" {float(np.mean(cell['mix_signed'])):+13.3f}")
        line += "  <-- coordinates help" if b < a else ""
        print(line)


report("by decile of the held-out truth (0 = weakest)", by_truth, "truth decile",
       signed=True)
report("by number of coordinate studies reporting near the voxel", by_count,
       "coord studies", signed=True)
print("\nRead the signed columns before the absolute ones. If images-only is systematically")
print("ABOVE the truth at low deciles, its advantage there is |x| upward bias from two studies")
print("and any shrinkage flatters the mixed arm for free. The interpretable cells are then the")
print("top deciles, where that bias is negligible.")
print("\nThe count-0 row is not a null check after all: a coordinate study that reported nothing")
print("near a voxel still enters the censoring term, so its SILENCE contributes there even with")
print("no focus reaching. That makes count-0 the largest and most interesting stratum rather")
print("than a control -- but it is also where the |x| confound bites hardest.")
