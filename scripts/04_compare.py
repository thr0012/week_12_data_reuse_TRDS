#!/usr/bin/env python3
"""
Step 4: compare all DMR sets side by side.

Reads results/*/enrichment.tsv and report.json and writes:
  results/comparison_table.tsv   one row per set x feature (observed, null, fold, p)
  results/comparison_heatmap.png log2 fold (observed / matched null) for key features

Usage: python3 04_compare.py results
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

root = sys.argv[1] if len(sys.argv) > 1 else "results"
ORDER = ["McBirney_F1", "McBirney_F2", "McBirney_F3",
         "Thorson_S2_lean", "Thorson_S3_kidney", "Thorson_S4_testis", "Thorson_S5_puberty", "Thorson_S6_multiple"]
NICE = {"McBirney_F1": "F1  atrazine vs control", "McBirney_F2": "F2  atrazine vs control",
        "McBirney_F3": "F3  atrazine vs control", "Thorson_S2_lean": "F3  lean vs none",
        "Thorson_S3_kidney": "F3  kidney vs none", "Thorson_S4_testis": "F3  testis vs none",
        "Thorson_S5_puberty": "F3  puberty vs none", "Thorson_S6_multiple": "F3  multiple vs none"}
FEATURES = [("H2", "Artifact-prone (union)"), ("H2", "Tandem repeat (TRF)"), ("H2", "rmsk: Simple_repeat"),
            ("H2", "rmsk: Low_complexity"), ("H2", "rmsk: Satellite"), ("H2", "family: ERVK"),
            ("H2", "rmsk: LINE"), ("H2", "CpG island"),
            ("H1", "Escaped (GD18 >= 30%)"), ("H1", "Mean GD18 methylation"),
            ("H1", "Erased then restored (GD18 < 10%, sperm >= 50%)")]

rows = []
for d in sorted(glob.glob(os.path.join(root, "*", "enrichment.tsv"))):
    name = os.path.basename(os.path.dirname(d))
    e = pd.read_csv(d, sep="\t")
    rep = json.load(open(os.path.join(os.path.dirname(d), "report.json")))
    e["set"] = name
    e["n_dmrs"] = rep.get("n_dmrs_canonical")
    e["pct_covered"] = rep.get("pct_dmrs_covered")
    rows.append(e)
if not rows:
    sys.exit(f"No enrichment.tsv found under {root}/")
allr = pd.concat(rows, ignore_index=True)
allr.to_csv(os.path.join(root, "comparison_table.tsv"), sep="\t", index=False, float_format="%.4g")

sets = [s for s in ORDER if s in set(allr["set"])] + sorted(set(allr["set"]) - set(ORDER))
print(f"Sets found: {sets}\n")

# console summary: the numbers for the slide
def g(s, f, c):
    v = allr[(allr["set"] == s) & (allr["feature"] == f)][c]
    return v.iloc[0] if len(v) else np.nan
print(f"{'set':22s} {'n':>5s} {'%cov':>5s} | {'artifact obs':>12s} {'null':>6s} {'fold':>5s} {'p':>7s} | {'TRF fold':>8s} | {'escape30 obs':>12s} {'null':>6s} {'p':>7s}")
for s in sets:
    a = "Artifact-prone (union)"; e = "Escaped (GD18 >= 30%)"
    r = allr[allr["set"] == s].iloc[0]
    print(f"{s:22s} {int(r['n_dmrs']):5d} {r['pct_covered']:5.1f} | {g(s,a,'observed'):12.3f} {g(s,a,'null_mean'):6.3f} "
          f"{g(s,a,'fold'):5.2f} {g(s,a,'p_two_sided'):7.4f} | {g(s,'Tandem repeat (TRF)','fold'):8.2f} | "
          f"{g(s,e,'observed'):12.3f} {g(s,e,'null_mean'):6.3f} {g(s,e,'p_two_sided'):7.4f}")

# heatmap of log2 fold
feats = [f for t, f in FEATURES if f in set(allr["feature"])]
L = np.full((len(feats), len(sets)), np.nan)
P = np.full_like(L, np.nan)
for j, s in enumerate(sets):
    sub = allr[allr["set"] == s].set_index("feature")
    for i, f in enumerate(feats):
        if f in sub.index and sub.at[f, "null_mean"] > 0 and sub.at[f, "observed"] > 0:
            L[i, j] = np.log2(sub.at[f, "fold"])
            P[i, j] = sub.at[f, "p_two_sided"]
        elif f in sub.index:
            P[i, j] = sub.at[f, "p_two_sided"]
            if sub.at[f, "observed"] == 0 and sub.at[f, "null_mean"] > 0:
                L[i, j] = -np.inf  # no DMR overlaps: maximal depletion

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

# diverging: blue (depleted) - neutral gray - orange (enriched)
cmap = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#e6e5e0", "#eb6834"])
fin = np.isfinite(L)
lim = max(1.0, np.nanmax(np.abs(L[fin])) if fin.any() else 1.0)
zero = np.isneginf(L)
L = np.where(zero, -lim, L)
fig, ax = plt.subplots(figsize=(max(9, 1.05 * len(sets) + 4.2), 0.48 * len(feats) + 2.2))
im = ax.imshow(L, cmap=cmap, norm=TwoSlopeNorm(0, -lim, lim), aspect="auto")
for i in range(len(feats)):
    for j in range(len(sets)):
        if np.isnan(P[i, j]):
            ax.text(j, i, "n/a", ha="center", va="center", fontsize=7.5, color="#8a8984")
            continue
        fold = 0 if zero[i, j] else 2 ** L[i, j]
        star = "***" if P[i, j] < 0.002 else "**" if P[i, j] < 0.01 else "*" if P[i, j] < 0.05 else ""
        ink = "#fcfcfb" if abs(L[i, j]) > 0.75 * lim else "#0b0b0b"
        ax.text(j, i - (0.12 if star else 0), f"{fold:.1f}x", ha="center", va="center", fontsize=8, color=ink)
        if star:
            ax.text(j, i + 0.2, star, ha="center", va="center", fontsize=8, color=ink)
ax.set_xticks(range(len(sets)), [NICE.get(s, s) for s in sets], rotation=35, ha="right", fontsize=8.5)
tag = {f: t for t, f in FEATURES}
ax.set_yticks(range(len(feats)), [f"{f.replace('rmsk: ', '').replace('family: ', '')}  [{tag[f]}]" for f in feats], fontsize=8.5)
n_h2 = sum(1 for t, f in FEATURES if t == "H2" and f in feats)
ax.axhline(n_h2 - 0.5, color="#fcfcfb", lw=4)
n_mcb = sum(1 for s in sets if s.startswith("McBirney"))
if 0 < n_mcb < len(sets):
    ax.axvline(n_mcb - 0.5, color="#fcfcfb", lw=4)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cb.set_label("log2 fold (DMRs / matched background)", fontsize=8, color="#52514e")
cb.ax.tick_params(labelsize=7.5)
fig.suptitle("Atrazine sperm DMRs vs matched background", x=0.01, y=0.985, ha="left", fontsize=11, color="#0b0b0b")
fig.text(0.01, 0.935, "Left: McBirney 2017, atrazine vs control lineage.  Right: Thorson 2020, disease vs no disease within atrazine F3.",
         fontsize=8, color="#52514e")
fig.text(0.01, 0.01, "Cells: fold;  * p<0.05  ** p<0.01  *** p<=0.002 (two-sided empirical, 1,000 matched sets).  "
         "H1 uses DMRs covered in GD18 data only.", fontsize=7, color="#52514e")
fig.tight_layout(rect=(0, 0.03, 1, 0.92))
out = os.path.join(root, "comparison_heatmap.png")
fig.savefig(out, dpi=200)
print(f"\nwrote {out} and {os.path.join(root, 'comparison_table.tsv')}")
