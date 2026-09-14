"""Do images with different fields of view get pooled per voxel, or amputated to an intersection?

Real collections never share a mask: field of view differs, signal dropout differs, and
normalisation trims different edges. Two designs are possible and only one is defensible --
intersect every image and analyse what survives, or estimate each voxel from whatever studies
cover it. Reading the code says CBES does the second. This checks it, and checks that the
per-voxel answer is *exactly* the inverse-variance mean of the covering studies rather than
approximately so.

Three studies, deliberately disjoint coverage, known constant effects:

    study A covers x < 7          g = 0.90, var = 1/20
    study B covers x > 2          g = 0.30, var = 1/40
    study C covers 4 <= x < 5     g = 0.60, var = 1/30

so voxels fall into bands covered by {A}, {A,B}, {A,B,C}, {B}, and by nothing at all.
"""
import json, sys, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from pathlib import Path
from nimare.meta.cbma import CBES
from nimare.studyset import Studyset

OUT = Path("/tmp/claude-0/ragged")
OUT.mkdir(parents=True, exist_ok=True)
SHAPE = (10, 6, 6)
AFF = np.diag([4.0, 4.0, 4.0, 1.0]); AFF[:3, 3] = -10.0
nib.save(nib.Nifti1Image(np.ones(SHAPE, np.int32), AFF), OUT / "mask.nii.gz")

SPEC = [("A", 0.90, 20, lambda x: x < 7),
        ("B", 0.30, 40, lambda x: x > 2),
        ("C", 0.60, 30, lambda x: (x >= 4) & (x < 5))]

xs = np.indices(SHAPE)[0]
studies, cover = [], {}
for name, g_val, n, rule in SPEC:
    mask = rule(xs)
    cover[name] = mask
    g = np.where(mask, g_val, np.nan).astype(np.float32)
    v = np.where(mask, 1.0 / n, np.nan).astype(np.float32)
    nib.save(nib.Nifti1Image(g, AFF), OUT / f"{name}_g.nii.gz")
    nib.save(nib.Nifti1Image(v, AFF), OUT / f"{name}_var.nii.gz")
    studies.append({
        "id": name, "name": name, "metadata": {"sample_sizes": [n]},
        "analyses": [{
            "id": f"{name}-1", "name": "1", "metadata": {"sample_sizes": [n]},
            # CBES requires a coordinates block even when every study supplies images;
            # an image study has its coordinates dropped before pooling, so this is inert.
            "points": [{"space": "MNI", "coordinates": [0.0, 0.0, 0.0],
                        "values": [{"kind": "Z", "value": 4.0}]}],
            "images": [
                {"url": str(OUT / f"{name}_g.nii.gz"), "filename": f"{name}_g.nii.gz",
                 "space": "MNI", "value_type": "g"},
                {"url": str(OUT / f"{name}_var.nii.gz"), "filename": f"{name}_var.nii.gz",
                 "space": "MNI", "value_type": "g_var"},
            ]}]})

source = {"id": "ragged", "name": "ragged", "studies": studies}
(OUT / "studyset.json").write_text(json.dumps(source))
ss = Studyset(str(OUT / "studyset.json"))

res = CBES(fwhm=10.0, mask=str(OUT / "mask.nii.gz"), use_images=True, null_method="none",
           selection_model="none", tau2_method="none").fit(ss)
g_map = res.get_map("g", return_type="array")
n_map = res.get_map("n_studies", return_type="array")
masker = res.estimator.masker
flat_x = masker.transform(nib.Nifti1Image(xs.astype(np.float32), AFF)).ravel()

ok = True
print(f"{'x band':>8s} {'covering':>12s} {'expected g':>11s} {'estimated':>11s} "
      f"{'n_studies':>10s}")
for x in range(SHAPE[0]):
    sel = flat_x == x
    if not sel.any():
        continue
    names = [nm for nm, _, _, rule in SPEC if bool(rule(np.array([x]))[0])]
    if names:
        # inverse-variance mean; with var = 1/n the weights are just n
        ns = np.array([n for nm, _, n, _ in SPEC if nm in names], dtype=float)
        vals = np.array([g for nm, g, _, _ in SPEC if nm in names])
        expected = float(np.sum(vals * ns) / np.sum(ns))
    else:
        expected = 0.0
    got = float(np.nanmean(g_map[sel]))
    n_got = float(np.nanmean(n_map[sel]))
    agree = abs(got - expected) < 1e-6 and abs(n_got - len(names)) < 1e-9
    ok &= agree
    print(f"{x:8d} {','.join(names) or '-':>12s} {expected:11.4f} {got:11.4f} {n_got:10.1f}"
          f"{'' if agree else '   <-- MISMATCH'}")

print()
print(f"finite everywhere: {bool(np.all(np.isfinite(g_map)))}")
print(f"voxels analysed:   {g_map.size} (an intersection mask would give "
      f"{int(np.sum(cover['A'] & cover['B'] & cover['C']))})")
print(f"\n{'PASS' if ok and np.all(np.isfinite(g_map)) else 'FAIL'}: "
      f"per-voxel pooling from available studies, exact to 1e-6")
