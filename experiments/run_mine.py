import glob, json, warnings
import numpy as np, nibabel as nib
warnings.simplefilter("ignore")
from nimare.studyset import Studyset
from nimare.meta.cbma import CBES
from nimare.correct import FDRCorrector
affine = np.array([[4.,0,0,-40.],[0,4.,0,-40.],[0,0,4.,-40.],[0,0,0,1.]])
mask = nib.Nifti1Image(np.ones((21,21,21), dtype=np.int32), affine)
unc, fdr, vfwe, csize, cmass, maxg = [], [], [], [], [], []
for path in sorted(glob.glob("/tmp/claude-0/cmp/m_*.json")):
    ss = Studyset(json.load(open(path)))
    est = CBES(fwhm=12.0, mask=mask, null_method="montecarlo", n_iters=100,
               cluster_threshold=0.01, seed=7)
    res = est.fit(ss)
    p = res.get_map("p", return_type="array")
    unc.append(float((p < 0.05).mean()))
    maxg.append(float(np.abs(res.get_map("g", return_type="array")).max()))
    fdr.append(bool((FDRCorrector(method="indep").transform(res)
                     .maps["p_corr-FDR_method-indep"] < 0.05).any()))
    maps, _, _ = est.correct_fwe_montecarlo(res, voxel_thresh=0.01, n_iters=100)
    vfwe.append(bool((10.0**-maps["logp_level-voxel"] < 0.05).any()))
    csize.append(bool((10.0**-maps["logp_desc-size_level-cluster"] < 0.05).any()))
    cmass.append(bool((10.0**-maps["logp_desc-mass_level-cluster"] < 0.05).any()))
print(f"MINE    uncorrected p<.05 {np.mean(unc):.3f}  FDR {np.mean(fdr):.2f}  "
      f"vFWE {np.mean(vfwe):.2f}  cFWE-size {np.mean(csize):.2f}  "
      f"cFWE-mass {np.mean(cmass):.2f}  max|g| {np.mean(maxg):.3f}", flush=True)
