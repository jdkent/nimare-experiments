"""Recovery of a known ground-truth effect size under the three selection models."""
import numpy as np
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import mm2vox

TRUTH = (0, 0, 0)

def value_at(res, name, xyz=TRUTH):
    img = res.get_map(name)
    ijk = mm2vox(np.array([xyz]), img.affine)[0]
    return float(img.get_fdata()[tuple(ijk)])

print(f"{'true g':>7} {'prev':>5} | {'none':>14} {'tobit':>14} {'zero-inflated':>16} {'pi':>5} {'rep':>6}")
for true_g in (0.3, 0.5, 0.8):
    for prevalence in (1.0, 0.5):
        out = {m: [] for m in ("none", "tobit", "zero-inflated")}
        pis, reps = [], []
        for seed in range(6):
            ss = create_effect_size_coordinate_studyset(
                [TRUTH], effect_sizes=true_g, n_studies=30, sample_size=(20, 40),
                tau=0.1, seed=seed, n_noise_foci=2, spatial_sd=6.0,
                prevalence=prevalence,
            )
            for model in out:
                r = CBES(fwhm=12.0, selection_model=model).fit(ss)
                out[model].append(value_at(r, "g"))
                if model == "zero-inflated":
                    pis.append(value_at(r, "prevalence"))
                    reps.append(value_at(r, "n_studies"))
        print(
            f"{true_g:>7} {prevalence:>5} | "
            + " ".join(f"{np.mean(out[m]):>7.3f}({np.mean(out[m])-true_g:+.3f})" for m in out)
            + f" {np.mean(pis):>5.2f} {np.mean(reps):>6.1f}"
        )
