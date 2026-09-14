"""Harvest every NeuroVault collection NeuroStore knows about, grouped by cognitive paradigm.

The point is a validation set that is not the pain collection: several independent studies of
the *same* paradigm, each contributing a group-level t or z map with a sample size. The maps
give an image-based truth; peaks extracted from them give the coordinates CBES sees.

NeuroStore's base-studies endpoint is the index -- ``data_type=image&level=group`` is 684
studies, which is all of them -- and each resolves to a NeuroVault collection through a
version's metadata URL. NeuroVault then supplies the per-image ``cognitive_paradigm_cogatlas``
that the grouping needs, which NeuroStore does not carry.
"""
import json, os, subprocess, sys, time
from collections import defaultdict

OUT = "/tmp/claude-0/paradigms"
os.makedirs(OUT, exist_ok=True)
CACHE = "/tmp/claude-0/cmp/neurostore_images.json"


def get(url, timeout=30, tries=2):
    for attempt in range(tries):
        raw = subprocess.run(["curl", "-s", "--max-time", str(timeout), url],
                             capture_output=True, text=True).stdout
        try:
            return json.loads(raw)
        except Exception:
            if attempt + 1 < tries:
                time.sleep(1.0)
    return None


base = json.load(open(CACHE))["results"]
candidates = [r for r in base
              if r.get("level") == "group" and (r.get("has_t_maps") or r.get("has_z_maps"))]
print(f"{len(candidates)} group-level base studies with t/z maps", flush=True)

# --- NeuroStore study -> NeuroVault collection id -------------------------------------------
resolved_path = f"{OUT}/resolved.json"
if os.path.exists(resolved_path):
    resolved = json.load(open(resolved_path))
else:
    resolved, seen = [], set()
    for i, rec in enumerate(candidates):
        for vid in (rec.get("versions") or []):
            if not isinstance(vid, str):
                continue
            study = get(f"https://neurostore.org/api/studies/{vid}")
            url = ((study or {}).get("metadata") or {}).get("url") or ""
            if "neurovault.org/collections/" in url:
                cid = url.rstrip("/").split("/")[-1]
                if cid.isdigit() and cid not in seen:
                    seen.add(cid)
                    resolved.append({"study": rec["id"], "collection": cid,
                                     "name": (rec.get("name") or "")[:80]})
                break
        if (i + 1) % 50 == 0:
            print(f"  resolved {len(resolved)} collections from {i + 1} studies", flush=True)
    json.dump(resolved, open(resolved_path, "w"))
print(f"{len(resolved)} distinct NeuroVault collections", flush=True)

# --- NeuroVault collection -> usable group maps, with their paradigm ------------------------
maps_path = f"{OUT}/maps.json"
if os.path.exists(maps_path):
    usable = json.load(open(maps_path))
else:
    usable = []
    for i, rec in enumerate(resolved):
        page = get(f"https://neurovault.org/api/collections/{rec['collection']}/images/?format=json")
        for img in ((page or {}).get("results") or []):
            # NeuroVault spells these "T map" / "Z map", not "t" / "z".
            raw_kind = str(img.get("map_type", "")).strip().lower()
            kind = "t" if raw_kind.startswith("t") else "z" if raw_kind.startswith("z") else ""
            n = img.get("number_of_subjects")
            level = str(img.get("analysis_level") or "group").lower()
            if not kind or not n or int(n) < 10 or level not in ("group", "", "none"):
                continue
            if img.get("not_mni") or img.get("is_thresholded"):
                continue
            paradigm = (img.get("cognitive_paradigm_cogatlas") or "").strip()
            if not paradigm or paradigm.lower() in ("none", "other", "null"):
                continue
            usable.append({
                "collection": rec["collection"], "study": rec["study"],
                "image": img.get("id"), "url": img.get("file"), "map_type": kind,
                "n": int(n), "paradigm": paradigm,
                "paradigm_id": img.get("cognitive_paradigm_cogatlas_id"),
                "contrast": str(img.get("cognitive_contrast_cogatlas") or "")[:60],
                "name": str(img.get("name", ""))[:70],
            })
        if (i + 1) % 25 == 0:
            print(f"  {i + 1}/{len(resolved)} collections scanned, {len(usable)} usable maps",
                  flush=True)
    json.dump(usable, open(maps_path, "w"))

print(f"\n{len(usable)} usable group maps with a named paradigm\n", flush=True)
by = defaultdict(set)
for m in usable:
    by[m["paradigm"]].add(m["collection"])
rows = sorted(by.items(), key=lambda kv: -len(kv[1]))
print(f"{'paradigm':>50} {'collections':>11}")
for paradigm, cols in rows[:20]:
    print(f"{paradigm[:50]:>50} {len(cols):11d}")
print(f"\nparadigms with >=10 collections: {sum(1 for _, c in rows if len(c) >= 10)}")
print(f"paradigms with >=6 collections:  {sum(1 for _, c in rows if len(c) >= 6)}")
