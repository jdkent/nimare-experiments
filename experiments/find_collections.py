"""Page the whole NeuroVault collection index and search it locally.

The API accepts no search or filter parameter -- name, search, q and name__icontains all come
back with the unfiltered count -- so the only way to find a collection by name is to pull the
index and grep it. 17856 collections at a thousand a page is eighteen requests.
"""
import json, subprocess, sys, time

TERMS = ("neuroscout", "hcp", "human connectome", "naturalistic", "narrative", "movie",
         "studyforrest", "budapest", "sherlock")
OUT = "/tmp/claude-0/nv_index.json"


def get(url, tries=3):
    for _ in range(tries):
        raw = subprocess.run(["curl", "-s", "--max-time", "120", url],
                             capture_output=True, text=True).stdout
        try:
            return json.loads(raw)
        except Exception:
            time.sleep(2)
    return None


rows, offset = [], 0
while True:
    page = get(f"https://neurovault.org/api/collections/?format=json&limit=1000&offset={offset}")
    if not page or not page.get("results"):
        break
    for r in page["results"]:
        rows.append({k: r.get(k) for k in
                     ("id", "name", "description", "owner_name", "number_of_images", "DOI",
                      "full_dataset_url", "paper_url")})
    offset += 1000
    print(f"  {len(rows)}/{page.get('count')} indexed", flush=True)
    if offset >= (page.get("count") or 0):
        break
json.dump(rows, open(OUT, "w"))
print(f"\nindexed {len(rows)} collections\n")

for term in TERMS:
    hits = [r for r in rows
            if term in " ".join(str(r.get(k) or "") for k in
                                ("name", "description", "owner_name")).lower()]
    with_images = [h for h in hits if (h.get("number_of_images") or 0) > 0]
    print(f"=== '{term}': {len(hits)} collections, {len(with_images)} with images ===")
    for h in sorted(with_images, key=lambda x: -(x.get("number_of_images") or 0))[:10]:
        print(f"  id={h['id']:>6} imgs={h['number_of_images']:>5} owner={str(h.get('owner_name'))[:18]:>18}"
              f"  {str(h.get('name'))[:60]}")
