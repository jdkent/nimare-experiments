"""Index HCP collection 4337: 18,070 single-subject contrast maps, by subject and contrast.

This is the collection that makes a proper validation possible. Every other reference used
today was derived from the same maps the coordinates came from, so "this study has an effect
here" was never observable independently of the magnitude. With per-subject maps over ~830
subjects per contrast, synthetic studies can be built from one set of subjects and the truth
measured on a disjoint set -- the selection that produces the peaks is then statistically
independent of the quantity being compared against.
"""
import json, subprocess, time

OUT = "/tmp/claude-0/hcp_index.json"


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
    page = get(f"https://neurovault.org/api/collections/4337/images/"
               f"?format=json&limit=1000&offset={offset}")
    if not page or not page.get("results"):
        break
    for r in page["results"]:
        name = str(r.get("name") or "")
        parts = name.split("_")
        if len(parts) < 2:
            continue
        rows.append({"id": r.get("id"), "file": r.get("file"), "subject": parts[0],
                     "contrast": "_".join(parts[1:]), "map_type": r.get("map_type"),
                     "paradigm": r.get("cognitive_paradigm_cogatlas")})
    offset += 1000
    print(f"  {len(rows)}/{page.get('count')} indexed", flush=True)
    if offset >= (page.get("count") or 0):
        break

json.dump(rows, open(OUT, "w"))
from collections import Counter
counts = Counter(r["contrast"] for r in rows)
print(f"\nindexed {len(rows)} maps, {len({r['subject'] for r in rows})} subjects, "
      f"{len(counts)} contrasts\n")
print(f"{'contrast':>28} {'subjects':>9}")
for contrast, n in counts.most_common(24):
    print(f"{contrast:>28} {n:>9}")
