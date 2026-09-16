"""On pain, at what image count does adding the coordinates stop helping and start hurting?

The headline pain result was measured at one image and at two: CBES `g` had 14% lower rmse than
an inverse-variance pool of the same images, and was far better centred at the top of the map
(+0.005 against +0.137). That is the regime the estimator was built for -- almost no images,
mostly tables.

The opposite regime has to fail eventually. Each coordinate study contributes only its reporting
indicator, read through two assumptions: a coverage radius, and a reporting threshold. Those
assumptions carry bias that does not shrink as studies accumulate. The image channel carries
noise that does. So there must be an image count past which the accumulated assumption bias
exceeds the noise the indicator removes, and CBES is worse than simply pooling the images.

This finds that crossing. The same splits and the same held-out reference are used at every
image count, so every comparison is paired: only the work half's division into image donors and
table donors changes. The table studies are never dropped -- at k images the remaining
`half - k` studies still contribute coordinates -- so what the sweep varies is the mixture, which
is the thing a meta-analyst actually controls.

Reported per image count: rmse and the signed error at the top decile, for CBES `g` and for the
image pool, with the paired difference and its p value. The crossing is where the rmse difference
changes sign.
"""
import logging, os, sys, warnings; warnings.simplefilter("ignore")
logging.getLogger("nimare").setLevel(logging.ERROR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import stats
from load_pain import load_pain
from validate_redesign import build, fit, masker, pooled, score, TABLES

N_SPLITS = int(os.environ.get("NSPLITS", 8))
COUNTS = [int(v) for v in os.environ.get("COUNTS", "1,2,3,4,5,6,7,8,9").split(",")]
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
    total = len(maps)
    half = total // 2
    counts = [k for k in COUNTS if k < half]
    print(f"NIDM pain: {total} studies, half = {half}, tables={TABLES}, "
          f"{N_SPLITS} splits, image counts {counts}\n", flush=True)

    # One permutation per split, reused at every image count, so the reference and the working
    # half are identical across counts and only the mixture changes.
    rng = np.random.default_rng(0)
    orders = [rng.permutation(total) for _ in range(N_SPLITS)]

    results, diag = {}, {}
    for k in counts:
        for split, order in enumerate(orders):
            work, held = order[:half], order[half:]
            truth = np.abs(pooled(held, maps, sizes))
            images, tables = list(work[:k]), list(work[k:])
            shipped = fit(build(images, tables, maps, sizes, published))
            if shipped is None:
                continue
            # The same collection with the selection term switched off: the image channel on its
            # own, fitted by the same code, so differencing the two isolates the coordinates.
            alone = fit(build(images, tables, maps, sizes, published), selection_model="none")
            only = np.abs(pooled(images, maps, sizes))
            use = shipped[1] & np.isfinite(only) & np.isfinite(truth)
            if use.sum() < 100:
                continue
            results.setdefault((k, "CBES g"), []).append(score(shipped[0], truth, use))
            results.setdefault((k, "images only"), []).append(score(only, truth, use))
            # Candidate scale-free predictor of the crossing: how much of the likelihood weight
            # the indicator actually carries, and how many subjects the images bring.
            share = shipped[4]
            strong = shipped[0] >= np.nanpercentile(shipped[0][use], 90)
            # Q, the Hausman statistic, entirely from what the estimator reports:
            #   Q = (g_full - g_images_only)^2 / se_full^2 * (1 - w) / w,   w = coordinate_share
            # which the algebra in micro_when_coordinates_hurt.py says should exceed 1 exactly
            # where the coordinates start costing accuracy. Nothing here needs the truth.
            q = np.full(shipped[0].shape, np.nan)
            acts = np.zeros(shipped[0].shape, dtype=bool)
            # Q is read only where the coordinate channel carries real weight. Algebraically
            # Q -> 0 as w -> 0, because delta = w (g_c - g_i); but delta is a difference of two
            # separately fitted maps, so it carries float noise that does not shrink with w, and
            # dividing that noise by w**2 manufactures values in the thousands. Restricting to
            # w >= Q_FLOOR keeps the statistic measuring disagreement rather than round-off.
            if share is not None and alone is not None:
                w = np.clip(share, 1e-6, 1 - 1e-6)
                delta = shipped[0] - alone[0]
                se = np.maximum(shipped[5], 1e-9)
                q = (delta / se) ** 2 * (1.0 - w) / w
                acts = share >= Q_FLOOR
            diag[k] = diag.get(k, [])
            diag[k].append((
                float(np.nanmedian(share[use])) if share is not None else np.nan,
                float(np.nanmedian(share[use & strong])) if share is not None else np.nan,
                float(sum(sizes[i] for i in images)),
                float(len(tables)),
                float(np.nanmedian(q[use & acts])) if (use & acts).any() else np.nan,
                float(np.nanmedian(q[use & acts & strong]))
                if (use & acts & strong).any() else np.nan,
                float(np.mean(acts[use])),
            ))
        done = len(results.get((k, "CBES g"), []))
        print(f"  {k} image(s): {done}/{N_SPLITS} splits fitted", flush=True)

    print(f"\n{'images':>7} {'tables':>7} | {'rmse CBES':>10} {'rmse pool':>10} "
          f"{'diff':>7} {'p':>7} | {'top CBES':>9} {'top pool':>9} {'diff':>7} {'p':>7}")
    for k in counts:
        a = np.array(results.get((k, "CBES g"), []), dtype=float)
        b = np.array(results.get((k, "images only"), []), dtype=float)
        if not len(a) or len(a) != len(b):
            continue
        row = [k, half - k]
        for j in (5, 4):                       # rmse, then signed error at the top decile
            d = a[:, j] - b[:, j]
            ok = np.isfinite(d)
            p = (stats.ttest_rel(a[ok, j], b[ok, j]).pvalue if ok.sum() > 1 else np.nan)
            row += [np.mean(a[ok, j]), np.mean(b[ok, j]), np.mean(d[ok]), p]
        print(f"{row[0]:7d} {row[1]:7d} | {row[2]:10.3f} {row[3]:10.3f} {row[4]:+7.3f} "
              f"{row[5]:7.3f} | {row[6]:+9.3f} {row[7]:+9.3f} {row[8]:+7.3f} {row[9]:7.3f}")

    print(f"\n{'images':>7} {'tables':>7} {'subjects':>9} {'share':>7} {'share@top':>10} "
          f"{'median Q':>9} {'Q@top':>9} {'w>=floor':>9}")
    for k in counts:
        if k not in diag:
            continue
        m = np.nanmean(np.array(diag[k]), axis=0)
        print(f"{k:7d} {int(m[3]):7d} {int(m[2]):9d} {m[0]:7.3f} {m[1]:10.3f} "
              f"{m[4]:9.3f} {m[5]:9.3f} {m[6]:9.2f}")
    print("\nQ is the Hausman statistic from the estimator's own maps; the algebra says the "
          "coordinates\nstop helping where it passes 1, with no reference to the truth.")

    print("\nA negative rmse difference means CBES beats the image pool; positive means the "
          "coordinates are costing accuracy.")
