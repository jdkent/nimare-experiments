"""Does `reporting.py` reproduce the coordinates a paper actually printed? Never checked.

The standing rule is to extract coordinates "the way papers produce them", and `reporting.py` is
that rule's implementation: correct for multiplicity, keep whole clusters, separate peaks by
8 mm, cap nothing. Every real-data result in this program used it -- including the pain
split-half comparison, which loads the collection's 267 transcribed coordinate rows and then
ignores them in favour of re-extraction.

So the proxy has never been compared against the thing it proxies for, although the NIDM pain
collection carries both for the same 21 studies. This is the descriptive half of that audit,
before any estimator runs:

  * how many peaks each study reports, published against extracted;
  * how far each published peak sits from the nearest extracted one, and the converse, which is
    the quantity that matters for a *silence* -- a voxel counted silent because extraction put
    the peak 15 mm from where the paper put it is a fabricated silence;
  * how much of the published table extraction recovers within the 8 mm separation radius it
    uses, and within the 20 mm coverage radius a silence asserts over.

Reported per scheme, because the scheme is the free choice and the point is which one -- if any
-- lands on what papers print.
"""
import os, sys, warnings; warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import nibabel as nib
from nilearn.datasets import load_mni152_brain_mask
from nilearn.maskers import NiftiMasker
from load_pain import load_pain
from reporting import report_peaks

SCHEMES = os.environ.get("SCHEMES", "cluster,fdr,fwe").split(",")
FOCUS = "max"

mask_img = load_mni152_brain_mask(resolution=4)
masker = NiftiMasker(mask_img).fit()
mask_bool = np.asarray(mask_img.get_fdata() > 0)
shape, affine = mask_img.shape, mask_img.affine
zooms = np.asarray(mask_img.header.get_zooms()[:3], dtype=float)


def nearest(a, b):
    """For each row of `a`, the distance in mm to the closest row of `b`."""
    if not len(a) or not len(b):
        return np.array([])
    return np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min(axis=1)


if __name__ == "__main__":
    ss = load_pain()
    coords = ss.coordinates
    published = {sid: g[["x", "y", "z"]].astype(float).to_numpy()
                 for sid, g in coords.groupby("study_id")}
    maps, sizes, ids = [], [], []
    for row, n in zip(ss.images.itertuples(), ss.sample_sizes()):
        maps.append(np.nan_to_num(masker.transform(row.z).ravel()))
        sizes.append(int(n if np.ndim(n) == 0 else np.asarray(n).ravel()[0]))
        ids.append(row.Index if isinstance(row.Index, str) else row.study_id)
    print(f"NIDM pain: {len(maps)} studies with maps, "
          f"{sum(len(v) for v in published.values())} published peaks "
          f"over {len(published)} studies\n")

    for scheme in SCHEMES:
        rows = []
        for sid, z in zip(ids, maps):
            pub = published.get(sid)
            if pub is None:
                continue
            found, height = report_peaks(z, mask_bool, shape, zooms, scheme=scheme, focus=FOCUS)
            ext = (np.array([nib.affines.apply_affine(affine, np.asarray(ijk, float))
                             for ijk, _ in found]) if found else np.zeros((0, 3)))
            d_pub = nearest(pub, ext)
            rows.append(dict(
                sid=sid, n_pub=len(pub), n_ext=len(ext), height=height,
                med=float(np.median(d_pub)) if len(d_pub) else np.nan,
                within8=float((d_pub <= 8).mean()) if len(d_pub) else np.nan,
                within20=float((d_pub <= 20).mean()) if len(d_pub) else np.nan))
        ok = [r for r in rows if r["n_ext"] > 0]
        silent_studies = sum(1 for r in rows if r["n_ext"] == 0)
        print(f"--- scheme={scheme}: {len(ok)}/{len(rows)} studies extracted anything "
              f"({silent_studies} would drop out entirely)")
        if not ok:
            print()
            continue
        print(f"    published peaks per study   median {np.median([r['n_pub'] for r in rows]):.0f}"
              f"   total {sum(r['n_pub'] for r in rows)}")
        print(f"    extracted peaks per study   median {np.median([r['n_ext'] for r in ok]):.0f}"
              f"   total {sum(r['n_ext'] for r in ok)}")
        print(f"    height threshold (z)        median {np.median([r['height'] for r in ok]):.2f}")
        print(f"    published peak to nearest extracted: median "
              f"{np.median([r['med'] for r in ok]):.1f} mm")
        print(f"    published peaks recovered within  8 mm: "
              f"{np.mean([r['within8'] for r in ok]):.2f}")
        print(f"    published peaks recovered within 20 mm: "
              f"{np.mean([r['within20'] for r in ok]):.2f}\n")
