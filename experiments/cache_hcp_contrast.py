"""Mask and cache an HCP contrast's per-subject maps, so the comparison scripts can np.load it.

`silence_only_hcp.py` reads `/tmp/claude-0/hcp/{CONTRAST}_masked.npy`. MOTOR_LH was already
cached; EMOTION_FACES was downloaded but never masked. Same procedure as
`hcp_heldout_validation.load_all`: resample each subject map to the 4 mm MNI mask and keep the
masked vector, skipping any file that fails to load rather than aborting the run.
"""
import os, sys, glob, warnings; warnings.simplefilter("ignore")
import numpy as np
import nibabel as nib
from nilearn.datasets import load_mni152_brain_mask
from nilearn.image import resample_to_img
from nilearn.maskers import NiftiMasker

CONTRAST = sys.argv[1] if len(sys.argv) > 1 else "EMOTION_FACES"
DATA = f"/tmp/claude-0/hcp/{CONTRAST}"
CACHE = f"/tmp/claude-0/hcp/{CONTRAST}_masked.npy"

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()

paths = sorted(glob.glob(f"{DATA}/*.nii.gz"))
print(f"{CONTRAST}: {len(paths)} files to mask", flush=True)
rows, skipped = [], 0
for i, path in enumerate(paths):
    try:
        img = nib.load(path)
        if img.ndim > 3:
            img = nib.Nifti1Image(np.asarray(img.dataobj)[..., 0], img.affine, img.header)
        rows.append(masker.transform(resample_to_img(
            img, mask_img, interpolation="continuous", force_resample=True,
            copy_header=True)).ravel().astype(np.float32))
    except Exception:
        skipped += 1
        continue
    if (i + 1) % 100 == 0:
        print(f"  {i + 1}/{len(paths)}", flush=True)
out = np.array(rows, dtype=np.float32)
np.save(CACHE, out)
print(f"cached {out.shape} to {CACHE}; skipped {skipped}")
