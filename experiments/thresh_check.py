"""Does a pooled reporting threshold beat a per-study one? The truth is known: z = 3.2905."""
import warnings; warnings.simplefilter("ignore")
import numpy as np
from nimare.generate import create_effect_size_coordinate_studyset

TRUE_Z = 3.2905267314919255

print(f"true reporting threshold z = {TRUE_Z:.3f}\n")
print(f"{'foci/study':>11s} {'pooled-min':>11s} {'study-min (median)':>19s} "
      f"{'study-min (mean)':>17s}")
for n_noise in (0, 1, 3, 8, 20):
    pooled, per_study_med, per_study_mean = [], [], []
    for seed in range(25):
        ss = create_effect_size_coordinate_studyset(
            [(0, 0, 0)], effect_sizes=0.8, n_studies=30, sample_size=(20, 40),
            seed=seed, n_noise_foci=n_noise, noise_extent=40.0)
        df = ss.coordinates
        z = np.abs(df["z_stat"].astype(float).to_numpy())
        pooled.append(z.min())
        mins = df.assign(az=z).groupby("id")["az"].min().to_numpy()
        per_study_med.append(np.median(mins))
        per_study_mean.append(mins.mean())
    per = len(df) / df["id"].nunique()
    print(f"{per:11.1f} {np.mean(pooled):11.3f} {np.mean(per_study_med):19.3f} "
          f"{np.mean(per_study_mean):17.3f}")

print("\nA threshold estimated too HIGH makes silence unsurprising, so the censoring term")
print("stops pulling the estimate down and the selection correction under-corrects.")
