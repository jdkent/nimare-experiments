import numpy as np, nibabel as nib
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import mm2vox

affine = np.array([[4.,0,0,-40.],[0,4.,0,-40.],[0,0,4.,-40.],[0,0,0,1.]])
mask = nib.Nifti1Image(np.ones((21,21,21), dtype=np.int32), affine)

def at(res, name):
    img = res.get_map(name); ijk = mm2vox(np.array([(0,0,0)]), img.affine)[0]
    return float(img.get_fdata()[tuple(ijk)])

for true_g in (0.3, 0.5, 0.8):
    for prev in (1.0, 0.5):
        rows = {m: [] for m in ("none","tobit","zero-inflated")}
        pis, ns = [], []
        for seed in range(8):
            ss = create_effect_size_coordinate_studyset(
                [(0,0,0)], effect_sizes=true_g, n_studies=30, sample_size=(20,40),
                tau=0.1, prevalence=prev, seed=seed, n_noise_foci=1, noise_extent=30.,
                spatial_sd=5.0)
            for m in rows:
                r = CBES(fwhm=12.0, mask=mask, selection_model=m).fit(ss)
                rows[m].append(at(r, "g"))
                if m == "zero-inflated":
                    pis.append(at(r,"prevalence")); ns.append(at(r,"n_studies"))
        print(f"g={true_g} prev={prev} | " +
              "  ".join(f"{m}={np.mean(rows[m]):.3f}({np.mean(rows[m])-true_g:+.3f})" for m in rows) +
              f"  pi={np.mean(pis):.2f} reporting={np.mean(ns):.1f}/30", flush=True)
