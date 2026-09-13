"""Where does a CBES fit actually spend its time? Whole-brain, realistic studyset."""
import cProfile, pstats, io, sys, time, warnings
warnings.simplefilter("ignore")
import numpy as np
from nimare.generate import create_effect_size_coordinate_studyset
from nimare.meta.cbma import CBES
from nimare.utils import get_template

N_ITERS = int(sys.argv[1]) if len(sys.argv) > 1 else 10

mask = get_template(space="mni152_2mm", mask="brain")
print(f"mask voxels: {int(mask.get_fdata().sum())}")

ss = create_effect_size_coordinate_studyset(
    [(0, 0, 0), (-40, -20, 50), (40, 20, -10)],
    effect_sizes=0.6, n_studies=40, sample_size=(20, 40), seed=1,
    n_noise_foci=6, noise_extent=60.0,
    threshold_z=[2.3263, 3.0902, 3.2905, 4.2649],
)
print(f"foci: {len(ss.coordinates)}, studies: {ss.coordinates['id'].nunique()}")

est = CBES(fwhm=10.0, mask=mask, n_iters=N_ITERS, seed=0,
           threshold="study-min-corrected", peak_bias="per-study")

start = time.perf_counter()
profiler = cProfile.Profile()
profiler.enable()
est.fit(ss)
profiler.disable()
elapsed = time.perf_counter() - start
print(f"\ntotal: {elapsed:.1f}s for fit + {N_ITERS} null iterations "
      f"({elapsed / max(N_ITERS, 1):.2f}s per iteration)\n")

stream = io.StringIO()
pstats.Stats(profiler, stream=stream).sort_stats("tottime").print_stats(28)
text = stream.getvalue()
print("\n".join(text.splitlines()[4:40]))
