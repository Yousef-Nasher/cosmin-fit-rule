"""Stratified random draw of 100 reviews from the sampling frame, with a replacement queue.

Run from study3/:  python3 code/draw_sample.py [output.csv]
Reads data/study3_frame.csv. Writes, for every record in the frame, its year stratum, whether it
was drawn or placed in the replacement queue, and its rank (draw rank for drawn records, queue
rank within the stratum for the others). Without an argument the table is printed."""
import csv
import random
import sys

SEED = 545749
N_DRAW = 100

rows = list(csv.DictReader(open("data/study3_frame.csv", encoding="utf-8")))
for r in rows:
    r["year_stratum"] = r["year_eff"]
rows.sort(key=lambda r: (r["year_stratum"], r["frame_id"]))

strata = {}
for r in rows:
    strata.setdefault(r["year_stratum"], []).append(r)
order = sorted(strata)
sizes = {s: len(strata[s]) for s in order}
total = sum(sizes.values())

# proportional allocation, largest remainder
raw = {s: sizes[s] / total * N_DRAW for s in order}
alloc = {s: int(raw[s] // 1) for s in order}
for s in sorted(order, key=lambda s: raw[s] - alloc[s], reverse=True)[:N_DRAW - sum(alloc.values())]:
    alloc[s] += 1

# one generator, strata shuffled in a fixed order
rng = random.Random(SEED)
drawn, queue = [], []
for s in order:
    pool = list(strata[s])
    rng.shuffle(pool)
    drawn.extend(pool[:alloc[s]])
    queue.extend((r, rank) for rank, r in enumerate(pool[alloc[s]:], 1))

out = [(r["frame_id"], r["year_stratum"], "drawn", i) for i, r in enumerate(drawn, 1)]
out += [(r["frame_id"], r["year_stratum"], "replacement queue", k) for r, k in queue]
fh = open(sys.argv[1], "w", newline="", encoding="utf-8") if len(sys.argv) > 1 else sys.stdout
w = csv.writer(fh, lineterminator="\n")
w.writerow(["frame_id", "year_stratum", "draw_list", "rank"])
w.writerows(sorted(out, key=lambda t: t[0]))
