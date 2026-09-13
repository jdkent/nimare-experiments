"""The same mixture CBES was given: real images for 5 studies, coordinates for the other 16.

ES-SDM was built for exactly this -- its title is "combines reported peak coordinates and
statistical parametric maps" -- and sdm_parse looks for ``<study>.nii.gz`` alongside the peak
files. A study supplied as an image needs no recreation from peaks, so its lower and upper
bounds should coincide.

The same five studies CBES used (pain_01 to pain_05, the first five by sorted id) are supplied
as t maps on SDM's own template grid; the remaining sixteen keep the peak files already derived
at z = 3.2905. Everything else matches the coordinate-only run so the two are comparable.
"""
import os, shutil, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, "/tmp/claude-0/-home-user-NiMARE/82bada38-540b-5f42-8ab4-86d2423ff73c/scratchpad")
import numpy as np
import nibabel as nib
from nilearn.image import resample_to_img
from load_pain import load_pain

SRC = "/tmp/claude-0/sdm_input"
OUT = "/tmp/claude-0/sdm_mixed"
TEMPLATE = "/tmp/claude-0/sdm/SdmPsiGui-linux64-v6.23/share/sdm_template.nii.gz"
N_IMAGE = 5

shutil.rmtree(OUT, ignore_errors=True)
os.makedirs(OUT)

ss = load_pain()
template = nib.load(TEMPLATE)
t_paths = {str(i): p for i, p in zip(ss.images["id"], ss.images["t"])}
sizes = {str(i): float(n) for i, n in zip(ss.ids, ss.sample_sizes())}
image_ids = sorted(t_paths)[:N_IMAGE]

for sid in image_ids:
    path = t_paths[sid]
    if path is None:
        raise SystemExit(f"no t map for {sid}")
    img = resample_to_img(nib.load(str(path)), template, interpolation="continuous")
    data = np.nan_to_num(np.asarray(img.dataobj, dtype=np.float32))
    out = nib.Nifti1Image(data, template.affine)
    # sdm_parse reads the NIfTI intent code to check the map really is a t statistic, and
    # refuses the study outright if it is missing. Resampling builds a fresh header, so the
    # intent has to be restored -- with the degrees of freedom, one-sample, that SDM will use.
    out.header.set_intent("t test", (sizes[sid] - 1.0,), name="")
    # sdm_parse also checks the NIfTI xform codes and refuses anything not flagged MNI or
    # Talairach. A fresh header defaults to 'aligned' (2), so set both codes to match SDM's own
    # template exactly: sform 4 (MNI152), qform 2.
    out.header.set_sform(template.affine, code=4)
    out.header.set_qform(template.affine, code=2)
    nib.save(out, f"{OUT}/{sid}.nii.gz")

kept = 0
for name in sorted(os.listdir(SRC)):
    if not name.endswith(".spm_mni.txt"):
        continue
    sid = name[: -len(".spm_mni.txt")]
    if sid in image_ids:
        continue          # this study comes in as an image instead
    shutil.copy(f"{SRC}/{name}", f"{OUT}/{name}")
    kept += 1

shutil.copy(f"{SRC}/sdm_table.txt", f"{OUT}/sdm_table.txt")
print(f"{len(image_ids)} studies as t images: {', '.join(image_ids)}")
print(f"{kept} studies as peak files")
print(f"table rows: {len(open(f'{OUT}/sdm_table.txt').read().splitlines()) - 1}")
