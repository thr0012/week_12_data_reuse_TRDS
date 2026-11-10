#!/usr/bin/env python3
"""Synthetic rn6-like inputs with planted signal, to check the pipeline recovers it.
Planted: ~35% of DMRs inside segmental duplications (H2), ~15% on GD18-escaped capture regions (H1)."""
import gzip
import os
import sys

import numpy as np
import pandas as pd

out = sys.argv[1] if len(sys.argv) > 1 else "testdata"
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(7)
CH = {"chr1": 4_000_000, "chr2": 3_000_000, "chrX": 2_500_000, "chrUn_test": 50_000}

# genome: CpG-depleted background + CpG islands; one N gap on chr1
seqs, cgi = {}, []
for c, L in CH.items():
    p = np.array([0.29, 0.21, 0.21, 0.29])
    s = rng.choice(np.frombuffer(b"ACGT", np.uint8), size=L, p=p)
    cg = np.flatnonzero((s[:-1] == 67) & (s[1:] == 71))
    drop = rng.random(cg.size) < 0.75
    s[cg[drop] + 1] = 65  # CpG depletion
    for st in rng.integers(0, L - 1500, size=L // 100_000):
        isl = rng.choice(np.frombuffer(b"CG", np.uint8), size=1000)
        s[st:st + 1000] = isl
        cgi.append((c, st, st + 1000))
    if c == "chr1":
        s[1_000_000:1_200_000] = 78
    seqs[c] = s.tobytes().decode()
with gzip.open(f"{out}/genome.fa.gz", "wt") as fh:
    for c, s in seqs.items():
        fh.write(f">{c}\n")
        for i in range(0, len(s), 60):
            fh.write(s[i:i + 60] + "\n")

def ucsc(rows, path, ncol_after=4):
    with gzip.open(path, "wt") as fh:
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

ucsc([(0, "chr1", 1_000_000, 1_200_000, 1, "N", 200000, "contig", "no")], f"{out}/gap.txt.gz")
ucsc([(0, c, s, e, "CpG: x", 0, 0, 0, 0, 0, 0) for c, s, e in cgi], f"{out}/cpgIslandExt.txt.gz")

# repeats (~35% of genome), classes with families
classes = [("LINE", "L1"), ("SINE", "B2"), ("LTR", "ERVK"), ("LTR", "ERV1"), ("DNA", "hAT"),
           ("Satellite", "Satellite"), ("Simple_repeat", "Simple_repeat"), ("Low_complexity", "Low_complexity")]
rm = []
for c, L in CH.items():
    n = int(L * 0.35 / 600)
    st = np.sort(rng.integers(0, L - 2000, size=n))
    ln = rng.integers(100, 1100, size=n)
    k = rng.choice(len(classes), size=n, p=[.3, .3, .1, .08, .1, .04, .05, .03])
    for s0, l0, kk in zip(st, ln, k):
        cl, fam = classes[kk]
        rm.append((0, 0, 0, 0, 0, c, s0, s0 + l0, 0, "+", "x", cl, fam, 0, 0, 0, 1))
ucsc(rm, f"{out}/rmsk.txt.gz")

# segmental duplications: 30 x 20 kb
sd = []
for _ in range(30):
    c = rng.choice(["chr1", "chr2", "chrX"])
    s0 = int(rng.integers(0, CH[c] - 25_000))
    if c == "chr1" and 980_000 < s0 < 1_200_000:
        continue
    sd.append((0, c, s0, s0 + 20_000))
ucsc(sd, f"{out}/genomicSuperDups.txt.gz")

# capture regions (~8%) and GD18 / sperm methylation; 12% of capture regions escape
cap = []
for c in ["chr1", "chr2", "chrX"]:
    for s0 in rng.integers(0, CH[c] - 4000, size=CH[c] // 40_000):
        if c == "chr1" and 990_000 < s0 < 1_200_000:
            continue
        cap.append((c, int(s0), int(s0 + rng.integers(1500, 4000)), rng.random() < 0.12))
def cov_file(path, kind, rep):
    rows = []
    for c, s0, e0, esc in cap:
        seq = np.frombuffer(seqs[c][s0:e0 + 1].encode(), np.uint8)
        for p in np.flatnonzero((seq[:-1] == 67) & (seq[1:] == 71)) + s0:
            if kind == "gd18":
                mu = rng.uniform(.6, .9) if esc else rng.uniform(0, .1)
            else:
                mu = rng.uniform(.7, .95)
            d = int(rng.integers(4, 15))
            m = rng.binomial(d, mu)
            pos1 = p + 1 if rng.random() < .5 else p + 2  # 1-based, + or - strand
            rows.append((c, pos1, pos1, round(100 * m / d, 2), m, d - m))
    pd.DataFrame(rows).sort_values([0, 1]).to_csv(path, sep="\t", header=False, index=False)
os.makedirs(f"{out}/gd18", exist_ok=True)
os.makedirs(f"{out}/sperm", exist_ok=True)
for r in (1, 2):
    cov_file(f"{out}/gd18/GD18_rep{r}.cov.gz", "gd18", r)
cov_file(f"{out}/sperm/sperm_rep1.cov.gz", "sperm", 1)

# DMRs: 400, 1-based, Skinner-style xlsx with a title row
dm = []
esc_regions = [x for x in cap if x[3]]
for i in range(400):
    u = rng.random()
    L = int(rng.choice([1000, 2000, 3000], p=[.5, .35, .15]))
    if u < .35 and sd:
        _, c, a, b = sd[rng.integers(len(sd))]
        s0 = int(rng.integers(a, b - L))
    elif u < .50:
        c, a, b, _ = esc_regions[rng.integers(len(esc_regions))]
        s0 = int(max(0, (a + b) // 2 - L // 2))
    else:
        c = rng.choice(["chr1", "chr2", "chrX"])
        s0 = int(rng.integers(0, CH[c] - L))
        if c == "chr1" and 995_000 < s0 < 1_200_000:
            s0 = 2_000_000
    dm.append((f"DMR{i+1}", c, s0 + 1, s0 + L, L))
dm.append(("DMR_unplaced", "chrUn_test", 101, 1100, 1000))
df = pd.DataFrame(dm, columns=["DMR Name", "Chromosome", "Start", "Stop", "Length"])
with pd.ExcelWriter(f"{out}/dmrs.xlsx") as xw:
    pd.DataFrame([["S1 Table. Synthetic F3 sperm DMRs (p<1e-4)"]]).to_excel(xw, header=False, index=False)
    df.to_excel(xw, startrow=2, index=False)
print(f"wrote test data to {out}: {len(dm)} DMRs, {len(cap)} capture regions, {len(sd)} segdups")
