#!/usr/bin/env python3
"""
Step 2d: is a DMR table on rn6?

The McBirney tables report the number of CpGs in each DMR ("cpg_num"). If the coordinates are
on rn6, counting CpGs in the rn6 sequence over the same interval should reproduce that number
(exactly or almost). If they are on rn5, counts will not match.

Usage: python3 03_check_dmr_assembly.py data/dmrs/McBirney2017_F3_DMRs.tsv data/genome/rn6.fa.gz
"""
import gzip
import sys

import numpy as np
import pandas as pd

tsv, fasta = sys.argv[1], sys.argv[2]
d = pd.read_csv(tsv, sep="\t")
need = set(d["chr"])
print(f"{len(d)} DMRs on {len(need)} chromosomes; reading those chromosomes from {fasta} ...")

seqs, name, chunks = {}, None, []
def keep(n, ch):
    if n in need:
        seqs[n] = "".join(ch).upper()
with gzip.open(fasta, "rt") as fh:
    for line in fh:
        if line.startswith(">"):
            keep(name, chunks)
            name, chunks = line[1:].split()[0], []
            if len(seqs) == len(need):
                break
        elif name in need:
            chunks.append(line.strip())
keep(name, chunks)
print("missing chromosomes in FASTA:", sorted(need - set(seqs)) or "none")

def count(c, s0, e0):  # CpGs whose C lies in [s0, e0), 0-based
    s = seqs.get(c)
    return np.nan if s is None or e0 >= len(s) else s[s0:e0 + 1].count("CG")

res = {}
for label, off in [("as given (1-based start)", -1), ("shifted +1", 0), ("shifted -1", -2)]:
    got = np.array([count(c, st + off, st + off + L) for c, st, L in zip(d["chr"], d["start"], d["length"])], float)
    ok = ~np.isnan(got)
    diff = np.abs(got[ok] - d["cpg_num"].to_numpy()[ok])
    r = np.corrcoef(got[ok], d["cpg_num"].to_numpy()[ok])[0, 1]
    res[label] = got
    print(f"{label:26s} exact match {100*np.mean(diff == 0):5.1f}%   within 1 CpG {100*np.mean(diff <= 1):5.1f}%   r = {r:.3f}")

got = res["as given (1-based start)"]
ok = ~np.isnan(got)
diff = np.abs(got[ok] - d["cpg_num"].to_numpy()[ok])
if np.mean(diff <= 1) > 0.9:
    print("\n=> Coordinates are on rn6. Use the table as is.")
elif np.corrcoef(got[ok], d["cpg_num"].to_numpy()[ok])[0, 1] < 0.3:
    print("\n=> Counts do not match rn6: probably rn5. LiftOver rn5 -> rn6 needed.")
else:
    print("\n=> Partial match. Send me this output.")
    bad = d[ok][diff > 1].head(5).assign(rn6_count=got[ok][diff > 1][:5])
    print(bad[["name", "chr", "start", "length", "cpg_num", "rn6_count"]].to_string(index=False))
