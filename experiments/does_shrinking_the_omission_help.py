"""Should the omission be shrunk instead of trusted, and does Q say when?

The reporting indicator is asymmetric on purpose. A study that named a voxel contributes the
event `|g| > c`; a study that named nothing within `coverage_radius` contributes `|g| <= c`.
Both are censoring events, but they are informative in opposite regimes: the silence pins mu
down when mu is small, while a report saturates -- once `mu sqrt(n)` is well past `c sqrt(n)`,
`P(|g| > c)` is nearly 1 and its derivative in mu is nearly 0. So at a reported voxel the report
itself carries almost no information about how large the effect is, and the confidence there is
carried by the images. That is a property of the likelihood, not a choice.

What *is* a choice is `coverage_radius`, and its direction is the opposite of what the name
suggests. The indicator is `-1` at the named voxel, `0` at a voxel within the radius of one of
that study's peaks but not named, and `+1` everywhere else the study examined. So the radius is an
*exclusion zone* around reported peaks, not the reach of the silence: **widening it removes
silence**. That trades two errors:

  * **wider**: fewer voxels get an indicator, so the coordinate channel carries less information
    (`var_c` rises), but fewer voxels are called silent that were really above threshold and
    merely sat too far from the reported focus (`b` falls).
  * **narrower**: more information, more misassignment.

`micro_when_coordinates_hurt.py` derives when a biased channel stops paying: the coordinates help
while `b^2 < var_c + var_i`. Since `var_i` falls as images accumulate, less bias can be afforded,
so the error-minimising radius should move *outward* with image count. `coordinate_share` is the
check that the direction is read right at all -- it has to fall as the radius widens.

Run on pain against a held-out reference, at a small image count and a large one.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from load_pain import load_pain
from validate_redesign import build, fit, masker, pooled, score, TABLES

N_SPLITS = int(os.environ.get("NSPLITS", 6))
RADII = [float(v) for v in os.environ.get("RADII", "4,6,8,10,14,20").split(",")]
COUNTS = [int(v) for v in os.environ.get("COUNTS", "1,5").split(",")]
#: Q is only meaningful where the coordinate channel carries weight; see the note at its use.
Q_FLOOR = float(os.environ.get("QFLOOR", 0.10))

if __name__ == "__main__":
    ss = load_pain()
    maps, sizes, study_ids = [], [], []
    for row, n in zip(ss.images.itertuples(), ss.sample_sizes()):
        maps.append(np.nan_to_num(masker.transform(row.z).ravel()))
        sizes.append(int(n if np.ndim(n) == 0 else np.asarray(n).ravel()[0]))
        study_ids.append(row.study_id)
    maps = np.array(maps)
    published = None
    if TABLES == "published":
        by_study = {sid: g[["x", "y", "z"]].astype(float).to_numpy()
                    for sid, g in ss.coordinates.groupby("study_id")}
        published = {i: by_study.get(sid) for i, sid in enumerate(study_ids)}
    total, half = len(maps), len(maps) // 2
    print(f"NIDM pain: {total} studies, tables={TABLES}, {N_SPLITS} splits, "
          f"radii {RADII} mm, image counts {COUNTS}\n", flush=True)

    rng = np.random.default_rng(0)
    orders = [rng.permutation(total) for _ in range(N_SPLITS)]

    out = {}
    for k in COUNTS:
        for radius in RADII:
            for order in orders:
                work, held = order[:half], order[half:]
                truth = np.abs(pooled(held, maps, sizes))
                images, tables = list(work[:k]), list(work[k:])
                collection = build(images, tables, maps, sizes, published)
                shipped = fit(collection, coverage_radius=radius)
                alone = fit(build(images, tables, maps, sizes, published),
                            selection_model="none", coverage_radius=radius)
                if shipped is None:
                    continue
                only = np.abs(pooled(images, maps, sizes))
                use = shipped[1] & np.isfinite(only) & np.isfinite(truth)
                if use.sum() < 100:
                    continue
                s_cbes = score(shipped[0], truth, use)
                s_pool = score(only, truth, use)
                q = np.nan
                # Q is read only where the coordinate channel carries real weight. Algebraically
                # Q -> 0 as w -> 0, because delta = w (g_c - g_i); but delta is a difference of two
                # separately fitted maps, so it carries float noise that does not shrink with w, and
                # dividing that noise by w**2 manufactures values in the thousands. Restricting to
                # w >= Q_FLOOR keeps the statistic measuring disagreement rather than round-off.
                if shipped[4] is not None and alone is not None:
                    w = np.clip(shipped[4], 1e-6, 1 - 1e-6)
                    se = np.maximum(shipped[5], 1e-9)
                    acts = use & (shipped[4] >= Q_FLOOR)
                    if acts.any():
                        q = float(np.nanmedian((((shipped[0] - alone[0]) / se) ** 2
                                                * (1.0 - w) / w)[acts]))
                out.setdefault((k, radius), []).append(
                    (s_cbes[5], s_pool[5], s_cbes[4], float(np.nanmedian(shipped[4][use])), q))
            print(f"  {k} image(s), radius {radius:.0f} mm: "
                  f"{len(out.get((k, radius), []))} splits", flush=True)

    for k in COUNTS:
        print(f"\n{k} image(s)")
        print(f"{'radius':>7} {'rmse CBES':>10} {'rmse pool':>10} {'err at top':>11} "
              f"{'share':>7} {'median Q':>9}")
        best, best_rmse = None, np.inf
        for radius in RADII:
            rows = np.array(out.get((k, radius), []), dtype=float)
            if not len(rows):
                continue
            m = np.nanmean(rows, axis=0)
            if m[0] < best_rmse:
                best, best_rmse = radius, m[0]
            print(f"{radius:7.0f} {m[0]:10.3f} {m[1]:10.3f} {m[2]:+11.3f} {m[3]:7.3f} "
                  f"{m[4]:9.3f}")
        print(f"  best radius by rmse: {best:.0f} mm")

    print("\nIf the statistic means what the algebra says, the best radius shrinks as images "
          "accumulate,\nand it sits near where median Q crosses 1.")
