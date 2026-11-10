#!/usr/bin/env python3
"""
Step 2b: check which assembly a Bismark .cov file was aligned to.

Idea: every line of a .cov file is a cytosine in CpG context. If the file was aligned to rn6,
almost every position should be a CpG in the rn6 sequence (C followed by G on the + strand,
or G preceded by C on the - strand). If it was aligned to another assembly (rn7), positions
are shifted and only ~1-2% will hit a CpG by chance.

Usage: python3 01_check_assembly.py data/genome/rn6.fa.gz data/methyl/gd18/<one file>.cov.gz [chrom]
"""
import gzip
import sys
from collections import Counter

fasta, cov = sys.argv[1], sys.argv[2]
chrom = sys.argv[3] if len(sys.argv) > 3 else "chr1"
want = chrom if chrom.startswith("chr") else "chr" + chrom
max_lines = 200_000

# 1) chromosome names used in the .cov file (first 2 million lines)
names = Counter()
with gzip.open(cov, "rt") as fh:
    for i, line in enumerate(fh):
        names[line.split("\t", 1)[0]] += 1
        if i >= 2_000_000:
            break
print("Chromosome names in .cov (first 2M lines):", dict(names.most_common(8)))

# 2) read one chromosome from the FASTA
seq, keep = [], False
with gzip.open(fasta, "rt") as fh:
    for line in fh:
        if line.startswith(">"):
            if keep:
                break
            name = line[1:].split()[0]
            keep = name == want or "chr" + name == want
            continue
        if keep:
            seq.append(line.strip())
seq = "".join(seq).upper()
if not seq:
    sys.exit(f"{want} not found in {fasta}")
print(f"{want} length in FASTA: {len(seq):,} bp")

# 3) test positions on that chromosome, and classify the ones that are not CpGs
cats = Counter()
cov_by_cat = {}
shift = Counter()
tot = out_of_range = 0
with gzip.open(cov, "rt") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        c, pos = f[0], f[1]
        if c != want and "chr" + c != want:
            continue
        p = int(pos) - 1  # .cov is 1-based
        if p + 2 >= len(seq) or p < 2:
            out_of_range += 1
            continue
        tot += 1
        depth = int(f[4]) + int(f[5]) if len(f) >= 6 else 0
        b = seq[p]
        if seq[p:p + 2] == "CG":
            k = "CpG, + strand (C)"
        elif seq[p - 1:p + 1] == "CG":
            k = "CpG, - strand (G)"
        elif b == "C":
            k = "C not followed by G (CHG/CHH, + strand)"
        elif b == "G":
            k = "G not preceded by C (CHG/CHH, - strand)"
        elif b == "N":
            k = "N in reference"
        else:
            k = f"{b} in reference (not a cytosine: coordinate mismatch or SNP)"
            for d in (-2, -1, 1, 2):  # would a small offset make it a CpG?
                q = p + d
                if seq[q:q + 2] == "CG" or seq[q - 1:q + 1] == "CG":
                    shift[d] += 1
        cats[k] += 1
        cov_by_cat.setdefault(k, []).append(depth)
        if tot >= max_lines:
            break

print(f"Positions tested on {want}: {tot:,}  (beyond chromosome end: {out_of_range:,})\n")
print(f"{'category':62s} {'n':>8s} {'%':>6s} {'median depth':>13s}")
for k, n in cats.most_common():
    d = sorted(cov_by_cat[k])
    print(f"{k:62s} {n:8,d} {100*n/tot:6.1f} {d[len(d)//2]:13d}")
if shift:
    print("\nA/T positions that become CpG with an offset (bp):", dict(shift))

cpg = cats["CpG, + strand (C)"] + cats["CpG, - strand (G)"]
cyt = cpg + cats["C not followed by G (CHG/CHH, + strand)"] + cats["G not preceded by C (CHG/CHH, - strand)"]
print(f"\nCpG: {100*cpg/tot:.1f}%   any cytosine (C or G): {100*cyt/tot:.1f}%")
if cyt / tot > 0.98:
    print("=> SAME assembly. Non-CpG positions are other cytosine contexts; the pipeline keeps only CpGs.")
elif cpg / tot < 0.2:
    print("=> DIFFERENT assembly: liftOver needed.")
else:
    print("=> UNCLEAR: many positions are not cytosines. Check the offset line above before continuing.")
