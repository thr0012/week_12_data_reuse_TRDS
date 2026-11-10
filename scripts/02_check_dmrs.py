#!/usr/bin/env python3
"""
Step 2c: inspect a DMR table before the analysis.

Checks:
  - finds the header row (Skinner tables have title rows above it) and the chr/start/stop columns
  - number of DMRs, chromosomes, length and CpG density
  - coordinates fit inside rn6 chromosomes
  - 1-based or 0-based: MEDIPS uses 1000 bp windows; with 1-based coordinates DMR starts fall on
    1, 1001, 2001... (start mod 1000 == 1) and stops on multiples of 1000

Usage: python3 02_check_dmrs.py data/dmrs/<table>.xlsx data/genome/rn6.chrom.sizes [sheet]
"""
import sys
from collections import Counter

import pandas as pd

path, sizes_path = sys.argv[1], sys.argv[2]
sheet = sys.argv[3] if len(sys.argv) > 3 else 0

if path.endswith((".xlsx", ".xls")):
    xl = pd.ExcelFile(path)
    print("Sheets in file:", xl.sheet_names)
    raw = pd.read_excel(path, sheet_name=sheet, header=None)
else:
    raw = pd.read_csv(path, sep=None, engine="python", header=None)

print("\nFirst 6 rows as read:")
print(raw.head(6).to_string(max_colwidth=25))

hdr = next((i for i in range(min(30, len(raw)))
            if any("chr" in str(x).lower() for x in raw.iloc[i])
            and any("start" in str(x).lower() for x in raw.iloc[i])), None)
if hdr is None:
    sys.exit("\nNo header row with 'chr' and 'start' found. Send me the first rows above.")
df = raw.iloc[hdr + 1:].copy()
df.columns = [str(x).strip() for x in raw.iloc[hdr]]
print(f"\nHeader found on row {hdr}. Columns: {list(df.columns)}")

low = {c.lower(): c for c in df.columns}
cchr = next(low[c] for c in low if c.startswith("chr"))
cs = next(low[c] for c in low if "start" in c)
ce = next(low[c] for c in low if "stop" in c or c == "end" or c.endswith("end"))
print(f"Using: chr = '{cchr}', start = '{cs}', end = '{ce}'")

d = pd.DataFrame({"chr": df[cchr].astype(str).str.strip(),
                  "start": pd.to_numeric(df[cs], errors="coerce"),
                  "end": pd.to_numeric(df[ce], errors="coerce")}).dropna()
d["chr"] = d["chr"].where(d["chr"].str.startswith("chr"), "chr" + d["chr"])
d["len"] = d["end"] - d["start"]
print(f"\nDMRs with valid coordinates: {len(d)} (rows after header: {len(df)})")
print("Chromosomes:", dict(Counter(d['chr']).most_common()))
print("Length (bp):", d["len"].describe()[["min", "25%", "50%", "75%", "max"]].to_dict())

dens_col = next((low[c] for c in low if "density" in c), None)
if dens_col:
    dens = pd.to_numeric(df[dens_col], errors="coerce")
    print(f"CpG density ('{dens_col}'):", dens.describe()[["min", "50%", "max"]].round(2).to_dict())

sizes = dict(l.split()[:2] for l in open(sizes_path))
sizes = {k: int(v) for k, v in sizes.items()}
unknown = sorted(set(d["chr"]) - set(sizes))
over = d[d["chr"].isin(sizes) & (d["end"] > d["chr"].map(sizes))]
print(f"\nChromosomes not in rn6: {unknown or 'none'}")
print(f"DMRs ending past the chromosome end: {len(over)}")

sm, em = (d["start"] % 1000).astype(int), (d["end"] % 1000).astype(int)
print("\nMost common start mod 1000:", dict(Counter(sm).most_common(3)))
print("Most common end   mod 1000:", dict(Counter(em).most_common(3)))
if (sm == 1).mean() > 0.8:
    print("=> 1-BASED coordinates (MEDIPS windows). Default of dmr_escape.py is correct.")
elif (sm == 0).mean() > 0.8:
    print("=> 0-BASED coordinates. Run dmr_escape.py with --zero-based.")
else:
    print("=> Window grid not obvious; keep the 1-based default and tell me what you see.")
