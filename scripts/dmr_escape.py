#!/usr/bin/env python3
"""
Do transgenerational sperm DMRs sit in regions that escape germline reprogramming?

Tests, against matched random backgrounds (same chromosome, same length, same CpG density bin):
  H1  escape: DMRs overlap regions that stay methylated in GD18 gonocytes
      (tested only inside the universe covered by the methylation data)
  H2  artifact: DMRs are enriched in segmental duplications / satellites / simple repeats
      and other repeat classes (tested genome-wide)

Inputs: genome FASTA (rn6), DMR table (csv/tsv/xlsx), UCSC annotation dumps,
        per-CpG methylation calls for GD18 gonocytes (optional: adult sperm).

Example:
  python dmr_escape.py --genome data/genome/rn6.fa.gz \
      --dmrs data/dmrs/mcbirney_F3.xlsx --label McBirney_F3 \
      --rmsk data/annot/rmsk.txt.gz --gap data/annot/gap.txt.gz \
      --segdup data/annot/genomicSuperDups.txt.gz --simple data/annot/simpleRepeat.txt.gz \
      --gd18 data/methyl/gd18/*.cov.gz --sperm data/methyl/sperm/*.cov.gz \
      --out results/McBirney_F3
"""
import argparse
import gzip
import json
import os
import re
import sys
import time

import numpy as np
import pandas as pd

CANON = re.compile(r"^chr([0-9]+|X|Y)$")


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def norm_chr(c):
    c = str(c).strip()
    if c.lower().startswith("chr"):
        c = "chr" + c[3:]
    else:
        c = "chr" + c
    return c


# --------------------------------------------------------------------------- genome

def index_genome(fasta, cache):
    """Per canonical chromosome: length, CpG positions (0-based C), merged N runs."""
    if cache and os.path.exists(cache):
        log(f"loading genome index from {cache}")
        z = np.load(cache, allow_pickle=True)
        return z["index"].item()
    log(f"indexing genome {fasta} (one pass, a few minutes for rn6)")
    index = {}

    def finish(name, chunks):
        if name is None or not CANON.match(name):
            return
        a = np.frombuffer("".join(chunks).encode(), dtype=np.uint8) & 0xDF  # upper-case
        cpg = np.flatnonzero((a[:-1] == 67) & (a[1:] == 71)).astype(np.int64)
        isn = (a == 78).astype(np.int8)
        d = np.diff(np.concatenate(([0], isn, [0])))
        ns, ne = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
        index[name] = {"len": int(a.size), "cpg": cpg, "nruns": (ns.astype(np.int64), ne.astype(np.int64))}
        log(f"  {name}: {a.size:,} bp, {cpg.size:,} CpG")

    name, chunks = None, []
    with opener(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                finish(name, chunks)
                name, chunks = norm_chr(line[1:].split()[0]), []
            else:
                chunks.append(line.strip())
    finish(name, chunks)
    if cache:
        np.savez(cache, index=np.array(index, dtype=object))
    return index


# --------------------------------------------------------------------------- interval tracks

class Track:
    """Merged, sorted intervals per chromosome with O(log n) covered-bp queries."""

    def __init__(self, df=None):
        self.d = {}
        if df is not None and len(df):
            for c, g in df.groupby("chr"):
                s, e = g["start"].to_numpy(np.int64), g["end"].to_numpy(np.int64)
                o = np.argsort(s, kind="stable")
                s, e = s[o], e[o]
                ms, me = [s[0]], [e[0]]
                for a, b in zip(s[1:], e[1:]):
                    if a <= me[-1]:
                        me[-1] = max(me[-1], b)
                    else:
                        ms.append(a)
                        me.append(b)
                ms, me = np.array(ms, np.int64), np.array(me, np.int64)
                cum = np.concatenate(([0], np.cumsum(me - ms)))
                self.d[c] = (ms, me, cum)

    @classmethod
    def from_arrays(cls, per_chr):
        t = cls()
        for c, (ms, me) in per_chr.items():
            cum = np.concatenate(([0], np.cumsum(me - ms)))
            t.d[c] = (ms, me, cum)
        return t

    def _cov(self, c, x):
        ms, me, cum = self.d[c]
        k = np.searchsorted(ms, x, side="left")
        out = np.zeros(x.shape, np.int64)
        m = k > 0
        km = k[m] - 1
        out[m] = cum[km] + np.minimum(me[km], x[m]) - ms[km]
        return out

    def bp(self, c, s, e):
        if c not in self.d:
            return np.zeros(np.shape(s), np.int64)
        return self._cov(c, e) - self._cov(c, s)

    def total(self):
        return int(sum(v[2][-1] for v in self.d.values()))


def read_ucsc(path, chrom_col, start_col, end_col, extra=None):
    """Read a UCSC database dump (tab, no header, leading 'bin' column handled by indices)."""
    usecols = [chrom_col, start_col, end_col] + ([*extra.values()] if extra else [])
    df = pd.read_csv(path, sep="\t", header=None, usecols=usecols, low_memory=False)
    ren = {chrom_col: "chr", start_col: "start", end_col: "end"}
    if extra:
        ren.update({v: k for k, v in extra.items()})
    df = df.rename(columns=ren)
    df["chr"] = df["chr"].map(norm_chr)
    return df[df["chr"].str.match(CANON)]


def build_tracks(args):
    tracks = {}
    if args.rmsk:
        log("reading RepeatMasker")
        rm = read_ucsc(args.rmsk, 5, 6, 7, {"repClass": 11, "repFamily": 12})
        rm["repClass"] = rm["repClass"].astype(str).str.replace("?", "", regex=False)
        rm["repFamily"] = rm["repFamily"].astype(str).str.replace("?", "", regex=False)
        tracks["Any repeat"] = Track(rm[~rm["repClass"].isin(["Unknown"])])
        for cl in ["LINE", "SINE", "LTR", "DNA", "Satellite", "Simple_repeat", "Low_complexity"]:
            sub = rm[rm["repClass"] == cl]
            if len(sub):
                tracks[f"rmsk: {cl}"] = Track(sub)
        for fam in ["ERVK", "ERV1", "ERVL", "ERVL-MaLR", "L1"]:
            sub = rm[rm["repFamily"] == fam]
            if len(sub):
                tracks[f"family: {fam}"] = Track(sub)
        art = [rm[rm["repClass"].isin(["Satellite", "Simple_repeat", "Low_complexity"])][["chr", "start", "end"]]]
    else:
        art = []
    if args.segdup:
        log("reading segmental duplications")
        sd = read_ucsc(args.segdup, 1, 2, 3)
        tracks["Segmental duplication"] = Track(sd)
        art.append(sd[["chr", "start", "end"]])
    if args.simple:
        log("reading simpleRepeat (TRF)")
        tr = read_ucsc(args.simple, 1, 2, 3)
        tracks["Tandem repeat (TRF)"] = Track(tr)
        art.append(tr[["chr", "start", "end"]])
    if args.cpgi:
        tracks["CpG island"] = Track(read_ucsc(args.cpgi, 1, 2, 3))
    if art:
        tracks["Artifact-prone (union)"] = Track(pd.concat(art))
    return tracks


# --------------------------------------------------------------------------- DMRs

def load_dmrs(path, sheet, one_based):
    if path.endswith((".xlsx", ".xls")):
        raw = pd.read_excel(path, sheet_name=sheet or 0, header=None)
    else:
        sep = "\t" if path.endswith((".tsv", ".txt", ".tsv.gz", ".txt.gz")) else ","
        raw = pd.read_csv(path, sep=sep, header=None)
    hdr = None
    for i in range(min(30, len(raw))):
        cells = [str(x).lower() for x in raw.iloc[i].tolist()]
        if any("chr" in x for x in cells) and any("start" in x for x in cells):
            hdr = i
            break
    if hdr is None:
        sys.exit(f"Could not find a header row with 'chr' and 'start' in {path}")
    df = raw.iloc[hdr + 1:].copy()
    df.columns = [str(x).strip() for x in raw.iloc[hdr]]
    low = {c.lower(): c for c in df.columns}
    cchr = next(low[c] for c in low if c.startswith("chr"))
    cs = next(low[c] for c in low if "start" in c)
    ce = next(low[c] for c in low if c in ("stop", "end") or "stop" in c or c.endswith("end"))
    out = pd.DataFrame({
        "chr": df[cchr].map(norm_chr),
        "start": pd.to_numeric(df[cs], errors="coerce"),
        "end": pd.to_numeric(df[ce], errors="coerce"),
    })
    name_col = next((low[c] for c in low if "name" in c or c.startswith("dmr")), None)
    out["name"] = df[name_col].astype(str).values if name_col else [f"DMR{i+1}" for i in range(len(out))]
    out = out.dropna(subset=["start", "end"])
    out["start"] = out["start"].astype(np.int64) - (1 if one_based else 0)
    out["end"] = out["end"].astype(np.int64)
    return out.reset_index(drop=True)


# --------------------------------------------------------------------------- methylation

def read_calls(path):
    """Return DataFrame chr,pos0,meth,unmeth (counts may be NaN for bedGraph)."""
    with opener(path) as fh:
        first = fh.readline()
    has_header = any(k in first for k in ("chr\t", "freqC", "chrBase", "track", "start"))
    df = pd.read_csv(path, sep="\t", header=0 if has_header else None, comment=None, low_memory=False)
    if any(str(c) == "freqC" for c in df.columns):           # methylKit
        cov = df["coverage"].astype(float)
        meth = np.round(cov * df["freqC"].astype(float) / 100)
        return pd.DataFrame({"chr": df["chr"].map(norm_chr), "pos0": df["base"].astype(np.int64) - 1,
                             "meth": meth, "unmeth": cov - meth})
    if str(df.iloc[0, 0]).startswith("track"):
        df = df.iloc[1:]
    n = df.shape[1]
    if n >= 7 and str(df.iloc[0, 5]) in ("CG", "CHG", "CHH"):  # bismark CpG_report
        df = df[df.iloc[:, 5] == "CG"]
        return pd.DataFrame({"chr": df.iloc[:, 0].map(norm_chr), "pos0": df.iloc[:, 1].astype(np.int64) - 1,
                             "meth": df.iloc[:, 3].astype(float), "unmeth": df.iloc[:, 4].astype(float)})
    if n >= 6:                                               # bismark .cov (1-based)
        return pd.DataFrame({"chr": df.iloc[:, 0].map(norm_chr), "pos0": df.iloc[:, 1].astype(np.int64) - 1,
                             "meth": df.iloc[:, 4].astype(float), "unmeth": df.iloc[:, 5].astype(float)})
    if n == 4:                                               # bedGraph (0-based, % only)
        pct = df.iloc[:, 3].astype(float)
        return pd.DataFrame({"chr": df.iloc[:, 0].map(norm_chr), "pos0": df.iloc[:, 1].astype(np.int64),
                             "meth": pct / 100, "unmeth": 1 - pct / 100})
    sys.exit(f"Unrecognized methylation format: {path}")


def load_methylation(files, genome, min_cov, label):
    """Pool replicates, collapse strands onto the CpG C position, filter coverage."""
    per = {}
    for f in files:
        log(f"  reading {label}: {os.path.basename(f)}")
        d = read_calls(f)
        d = d[d["chr"].isin(genome.keys())]
        for c, g in d.groupby("chr"):
            cpg = genome[c]["cpg"]
            p = g["pos0"].to_numpy()
            i = np.searchsorted(cpg, p)
            i = np.clip(i, 0, cpg.size - 1)
            hit_plus = cpg[i] == p
            j = np.clip(np.searchsorted(cpg, p - 1), 0, cpg.size - 1)
            hit_minus = (~hit_plus) & (cpg[j] == p - 1)
            idx = np.where(hit_plus, i, np.where(hit_minus, j, -1))
            keep = idx >= 0
            per.setdefault(c, []).append(pd.DataFrame({"k": idx[keep], "m": g["meth"].to_numpy()[keep],
                                                       "u": g["unmeth"].to_numpy()[keep]}))
    out = {}
    n_tot = 0
    for c, parts in per.items():
        agg = pd.concat(parts).groupby("k", sort=True)[["m", "u"]].sum()
        cov = (agg["m"] + agg["u"]).to_numpy()
        ok = cov >= min_cov
        k = agg.index.to_numpy()[ok]
        frac = agg["m"].to_numpy()[ok] / cov[ok]
        pos = genome[c]["cpg"][k]
        out[c] = (pos, frac, np.concatenate(([0.0], np.cumsum(frac))))
        n_tot += pos.size
    log(f"  {label}: {n_tot:,} CpGs pass coverage >= {min_cov}")
    return out


def region_meth(meth, c, s, e):
    """Number of covered CpGs and mean methylation (fraction) in [s, e)."""
    if c not in meth:
        z = np.zeros(np.shape(s), np.int64)
        return z, np.full(np.shape(s), np.nan)
    pos, _, cum = meth[c]
    a, b = np.searchsorted(pos, s), np.searchsorted(pos, e)
    n = b - a
    with np.errstate(invalid="ignore", divide="ignore"):
        m = np.where(n > 0, (cum[b] - cum[a]) / np.maximum(n, 1), np.nan)
    return n, m


# --------------------------------------------------------------------------- features & sampling

def features(df, genome, ntrack):
    dens = np.zeros(len(df))
    nfrac = np.zeros(len(df))
    for c, g in df.groupby("chr"):
        i = g.index.to_numpy()
        s, e = g["start"].to_numpy(), g["end"].to_numpy()
        cpg = genome[c]["cpg"]
        n = np.searchsorted(cpg, e - 1) - np.searchsorted(cpg, s)
        dens[i] = 100 * n / (e - s)
        nfrac[i] = ntrack.bp(c, s, e) / (e - s)
    return dens, nfrac


def sample_matched(dmr, genome, ntrack, bins, n_perm, rng, meth=None, min_cpg=3, oversample=8,
                   max_rounds=6):
    """For each DMR, draw n_perm random regions on the same chromosome with the same length and
    CpG density bin, <5% N. If meth is given, candidates are anchored on covered CpGs and must
    contain >= min_cpg covered CpGs (covered universe for H1).
    Returns arrays (n_dmr, n_perm) of starts and a count of DMRs that needed resampling."""
    starts = np.full((len(dmr), n_perm), -1, np.int64)
    short = 0
    for c, g in dmr.groupby("chr"):
        idx = g.index.to_numpy()
        L = (g["end"] - g["start"]).to_numpy()
        tb = g["bin"].to_numpy()
        clen = genome[c]["len"]
        cpg = genome[c]["cpg"]
        anchors = meth[c][0] if meth is not None and c in meth else None
        if meth is not None and (anchors is None or anchors.size == 0):
            continue
        filled = np.zeros(len(idx), np.int64)
        for r in range(max_rounds):
            need = filled < n_perm
            if not need.any():
                break
            M = n_perm * oversample * (2 ** r)
            for k in np.flatnonzero(need):
                l = L[k]
                if anchors is None:
                    s = rng.integers(0, max(1, clen - l), size=M)
                else:
                    s = rng.choice(anchors, size=M) - rng.integers(0, l, size=M)
                    s = np.clip(s, 0, max(0, clen - l))
                e = s + l
                n = np.searchsorted(cpg, e - 1) - np.searchsorted(cpg, s)
                b = np.digitize(100 * n / l, bins) - 1
                ok = (b == tb[k]) & (ntrack.bp(c, s, e) / l < 0.05)
                if meth is not None:
                    ok &= region_meth(meth, c, s, e)[0] >= min_cpg
                s = s[ok]
                take = min(n_perm - filled[k], s.size)
                starts[idx[k], filled[k]:filled[k] + take] = s[:take]
                filled[k] += take
        lack = filled < n_perm
        for k in np.flatnonzero(lack & (filled > 0)):  # top up with replacement
            have = starts[idx[k], :filled[k]]
            starts[idx[k], filled[k]:] = rng.choice(have, n_perm - filled[k])
        short += int(lack.sum())
    return starts, short


# --------------------------------------------------------------------------- statistics

def overlap_hits(df_chr, starts, ends, track, min_frac):
    out = np.zeros(starts.shape, bool)
    for c, rows in df_chr.items():
        s, e = starts[rows], ends[rows]
        out[rows] = track.bp(c, s.ravel(), e.ravel()).reshape(s.shape) >= min_frac * (e - s)
    return out


def emp(obs, null):
    null = np.asarray(null)
    mu = np.nanmean(null)
    lo, hi = np.nanpercentile(null, [2.5, 97.5])
    p_hi = (1 + np.sum(null >= obs)) / (null.size + 1)
    p_lo = (1 + np.sum(null <= obs)) / (null.size + 1)
    fold = obs / mu if mu > 0 else np.nan
    return dict(observed=obs, null_mean=mu, null_lo95=lo, null_hi95=hi, fold=fold,
                p_enriched=p_hi, p_depleted=p_lo, p_two_sided=min(1.0, 2 * min(p_hi, p_lo)))


# --------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--genome-cache", default=None, help="npz cache for the genome index (default: <genome>.idx.npz)")
    ap.add_argument("--dmrs", required=True)
    ap.add_argument("--sheet", default=None)
    ap.add_argument("--zero-based", action="store_true", help="DMR starts are already 0-based (default: 1-based)")
    ap.add_argument("--label", default="DMRs")
    ap.add_argument("--rmsk")
    ap.add_argument("--gap")
    ap.add_argument("--segdup")
    ap.add_argument("--simple")
    ap.add_argument("--cpgi")
    ap.add_argument("--gd18", nargs="*", default=[], help="GD18 gonocyte methylation call files (replicates pooled)")
    ap.add_argument("--sperm", nargs="*", default=[], help="adult sperm methylation call files (optional)")
    ap.add_argument("--min-cov", type=int, default=5, help="pooled read depth per CpG")
    ap.add_argument("--min-cpg", type=int, default=3, help="covered CpGs needed to call a region covered")
    ap.add_argument("--escape-thr", type=float, nargs="+", default=[0.3, 0.5, 0.7])
    ap.add_argument("--erased-thr", type=float, default=0.1)
    ap.add_argument("--min-overlap", type=float, default=0.1, help="fraction of region bp overlapping a track to count as a hit")
    ap.add_argument("--n-bins", type=int, default=5, help="CpG density bins (quantiles of DMRs)")
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    genome = index_genome(args.genome, args.genome_cache or (args.genome + ".idx.npz"))

    # N runs from the sequence, plus UCSC gaps if provided
    nr = {c: v["nruns"] for c, v in genome.items()}
    if args.gap:
        g = read_ucsc(args.gap, 1, 2, 3)
        ntrack = Track(pd.concat([g[["chr", "start", "end"]],
                                  pd.DataFrame([(c, s, e) for c, (ss, ee) in nr.items() for s, e in zip(ss, ee)],
                                               columns=["chr", "start", "end"])]))
    else:
        ntrack = Track.from_arrays(nr)
    tracks = build_tracks(args)

    # ---- DMRs
    dmr = load_dmrs(args.dmrs, args.sheet, one_based=not args.zero_based)
    n_in = len(dmr)
    dmr = dmr[dmr["chr"].isin(genome.keys())].reset_index(drop=True)
    dmr["start"] = dmr["start"].clip(lower=0)
    dmr["end"] = np.minimum(dmr["end"], dmr["chr"].map(lambda c: genome[c]["len"]))
    log(f"{args.label}: {n_in} DMRs read, {len(dmr)} on canonical chromosomes")
    dmr["cpg_per_100bp"], dmr["N_frac"] = features(dmr, genome, ntrack)
    q = np.quantile(dmr["cpg_per_100bp"], np.linspace(0, 1, args.n_bins + 1)[1:-1])
    bins = np.unique(np.concatenate(([-np.inf], q, [np.inf])))
    dmr["bin"] = np.digitize(dmr["cpg_per_100bp"], bins) - 1
    by_chr = {c: g.index.to_numpy() for c, g in dmr.groupby("chr")}
    S, E = dmr["start"].to_numpy(), dmr["end"].to_numpy()

    report = {"label": args.label, "n_dmrs_input": n_in, "n_dmrs_canonical": len(dmr),
              "dmr_length_median": float(np.median(E - S)),
              "dmr_cpg_per_100bp_median": float(dmr["cpg_per_100bp"].median()),
              "cpg_density_bin_edges": [float(x) for x in bins[1:-1]], "n_perm": args.n_perm}

    # ---- H2: genome-wide matched background
    rows = []
    if tracks:
        log(f"H2: sampling {args.n_perm} genome-wide matched sets")
        bs, short = sample_matched(dmr, genome, ntrack, bins, args.n_perm, rng)
        report["H2_dmrs_topped_up"] = short
        be = bs + (E - S)[:, None]
        for name, tr in tracks.items():
            obs_hit = overlap_hits(by_chr, S[:, None], E[:, None], tr, args.min_overlap)[:, 0]
            dmr[f"hit: {name}"] = obs_hit
            nh = overlap_hits(by_chr, bs, be, tr, args.min_overlap)
            valid = bs >= 0
            null = np.array([nh[valid[:, j], j].mean() for j in range(args.n_perm)])
            r = emp(obs_hit.mean(), null)
            r.update(test="H2", feature=name, n=len(dmr), genome_fraction=tr.total() / sum(v["len"] for v in genome.values()))
            rows.append(r)
            log(f"  {name:28s} obs {r['observed']:.3f}  null {r['null_mean']:.3f}  fold {r['fold']:.2f}  p_enr {r['p_enriched']:.4f}")

    # ---- Methylation (H1)
    meth_rows, meth_dist = [], {}
    if args.gd18:
        log("loading GD18 methylation")
        gd = load_methylation(args.gd18, genome, args.min_cov, "GD18")
        sp = load_methylation(args.sperm, genome, args.min_cov, "sperm") if args.sperm else None
        n_cov, m_gd = np.zeros(len(dmr), np.int64), np.full(len(dmr), np.nan)
        m_sp = np.full(len(dmr), np.nan)
        for c, i in by_chr.items():
            n_cov[i], m_gd[i] = region_meth(gd, c, S[i], E[i])
            if sp:
                m_sp[i] = region_meth(sp, c, S[i], E[i])[1]
        dmr["GD18_n_cpg"], dmr["GD18_meth"] = n_cov, m_gd
        if sp:
            dmr["sperm_meth"] = m_sp
        covered = n_cov >= args.min_cpg
        report.update(n_dmrs_covered=int(covered.sum()), pct_dmrs_covered=round(100 * covered.mean(), 2),
                      gd18_cpgs=int(sum(v[0].size for v in gd.values())))
        log(f"STEP 0 coverage: {covered.sum()} / {len(dmr)} DMRs ({100*covered.mean():.1f}%) have >= {args.min_cpg} covered CpGs")
        if covered.sum() >= 10:
            sub = dmr[covered].reset_index(drop=True)
            log(f"H1: sampling {args.n_perm} matched sets inside the covered universe")
            cs, short = sample_matched(sub, genome, ntrack, bins, args.n_perm, rng, meth=gd, min_cpg=args.min_cpg)
            report["H1_dmrs_topped_up"] = short
            Ls = (sub["end"] - sub["start"]).to_numpy()
            ce = cs + Ls[:, None]
            nm = np.full(cs.shape, np.nan)
            nsp = np.full(cs.shape, np.nan)
            for c, g in sub.groupby("chr"):
                i = g.index.to_numpy()
                ok = cs[i] >= 0
                s_, e_ = cs[i][ok], ce[i][ok]
                tmp = np.full(cs[i].shape, np.nan)
                tmp[ok] = region_meth(gd, c, s_, e_)[1]
                nm[i] = tmp
                if sp:
                    tmp2 = np.full(cs[i].shape, np.nan)
                    tmp2[ok] = region_meth(sp, c, s_, e_)[1]
                    nsp[i] = tmp2
            om = sub["GD18_meth"].to_numpy()
            r = emp(np.nanmean(om), np.nanmean(nm, axis=0))
            r.update(test="H1", feature="Mean GD18 methylation", n=len(sub))
            meth_rows.append(r)
            for t in args.escape_thr:
                r = emp(np.mean(om >= t), np.nanmean(nm >= t, axis=0))
                r.update(test="H1", feature=f"Escaped (GD18 >= {int(t*100)}%)", n=len(sub))
                meth_rows.append(r)
            if sp:
                osp = sub["sperm_meth"].to_numpy()
                obs = np.mean((om < args.erased_thr) & (osp >= 0.5))
                null = np.nanmean((nm < args.erased_thr) & (nsp >= 0.5), axis=0)
                r = emp(obs, null)
                r.update(test="H1", feature=f"Erased then restored (GD18 < {int(args.erased_thr*100)}%, sperm >= 50%)", n=len(sub))
                meth_rows.append(r)
            for r in meth_rows:
                log(f"  {r['feature']:52s} obs {r['observed']:.3f}  null {r['null_mean']:.3f}  p2 {r['p_two_sided']:.4f}")
            meth_dist = {"obs": om, "null": nm[~np.isnan(nm)]}
            np.savez_compressed(os.path.join(args.out, "gd18_meth_distributions.npz"), **meth_dist)
        else:
            log("Too few covered DMRs for H1 (< 10). See README step 0: switch the escape reference.")

    # ---- write outputs
    res = pd.DataFrame(meth_rows + rows)
    if len(res):
        cols = ["test", "feature", "n", "observed", "null_mean", "null_lo95", "null_hi95", "fold",
                "p_enriched", "p_depleted", "p_two_sided"] + (["genome_fraction"] if "genome_fraction" in res else [])
        res = res[cols]
        res.to_csv(os.path.join(args.out, "enrichment.tsv"), sep="\t", index=False, float_format="%.4g")
    dmr.drop(columns=["bin"]).to_csv(os.path.join(args.out, "dmrs_annotated.tsv"), sep="\t", index=False, float_format="%.4g")
    with open(os.path.join(args.out, "report.json"), "w") as fh:
        json.dump(report, fh, indent=2, default=float)
    plot(res, meth_dist, args)
    log(f"done: {args.out}")


# --------------------------------------------------------------------------- figures

BLUE, GRAY, INK, MUTED = "#2a78d6", "#a9a8a2", "#0b0b0b", "#52514e"


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#d6d5cf")
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.xaxis.label.set_color(MUTED)
    ax.yaxis.label.set_color(MUTED)


def plot(res, meth_dist, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    if len(res):
        panels = [(t, res[res["test"] == t].iloc[::-1].reset_index(drop=True)) for t in ("H1", "H2")]
        panels = [(t, r) for t, r in panels if len(r)]
        titles = {"H1": "H1  Escape: GD18 gonocyte methylation (covered universe)",
                  "H2": "H2  Artifact-prone and repeat regions (genome-wide)"}
        fig, axes = plt.subplots(len(panels), 1, figsize=(8, 0.36 * len(res) + 1.2 * len(panels) + 0.6),
                                 gridspec_kw={"height_ratios": [len(r) + 1 for _, r in panels]}, squeeze=False)
        for ax, (t, r) in zip(axes[:, 0], panels):
            y = np.arange(len(r))
            ax.hlines(y, r["null_lo95"], r["null_hi95"], color=GRAY, lw=6, alpha=0.6, capstyle="round")
            ax.plot(r["null_mean"], y, "|", color=MUTED, ms=12, mew=1.5)
            ax.plot(r["observed"], y, "o", color=BLUE, ms=8, mec="#fcfcfb", mew=1.5, zorder=3)
            for yi, (_, row) in zip(y, r.iterrows()):
                p = row["p_enriched"] if row["observed"] >= row["null_mean"] else row["p_depleted"]
                star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
                ax.text(1.01, yi, f"{row['fold']:.2f}x {star}", transform=ax.get_yaxis_transform(),
                        va="center", fontsize=8.5, color=INK)
            ax.set_yticks(y, r["feature"], fontsize=8.5)
            ax.set_ylim(-0.7, len(r) - 0.3)
            ax.set_xlim(0, 1)
            ax.set_title(f"{titles[t]}  (n = {int(r['n'].iloc[0])} DMRs)", loc="left", fontsize=9.5, color=INK)
            style(ax)
            ax.grid(axis="x", color="#ecebe6", lw=0.8)
            ax.set_axisbelow(True)
        has_h1 = any(t == "H1" for t, _ in panels)
        axes[-1, 0].set_xlabel("Fraction of DMRs" + ("  (top row of H1: mean methylation)" if has_h1 else ""))
        fig.suptitle(f"{args.label}: observed DMRs vs matched background", x=0.02, y=0.995, ha="left", fontsize=11.5, color=INK)
        fig.legend(handles=[Line2D([], [], marker="o", ls="", color=BLUE, label="Observed DMRs"),
                            Line2D([], [], color=GRAY, lw=6, alpha=0.6, label="Matched null, 95% interval"),
                            Line2D([], [], marker="|", ls="", color=MUTED, ms=10, label="Null mean")],
                   frameon=False, fontsize=8, loc="upper left", ncol=3, bbox_to_anchor=(0.01, 0.975))
        fig.text(0.99, 0.005, "fold = observed / null mean;  * p<0.05  ** p<0.01  *** p<0.001 (empirical)",
                 ha="right", fontsize=7.5, color=MUTED)
        fig.tight_layout(rect=(0, 0.02, 1, 0.94))
        fig.savefig(os.path.join(args.out, "fig1_enrichment.png"), dpi=200)
        plt.close(fig)

    if meth_dist:
        fig, ax = plt.subplots(figsize=(6, 3.6))
        b = np.linspace(0, 1, 21)
        ax.hist(meth_dist["null"], bins=b, density=True, color=GRAY, alpha=0.7, label="Matched background",
                edgecolor="#fcfcfb", lw=2)
        ax.hist(meth_dist["obs"], bins=b, density=True, histtype="step", color=BLUE, lw=2, label="DMRs")
        for t in args.escape_thr:
            ax.axvline(t, color=MUTED, lw=0.8, ls=":")
        ax.set_xlabel("Mean CpG methylation in GD18 gonocytes")
        ax.set_ylabel("Density")
        ax.set_title(f"{args.label}: GD18 methylation at DMRs", loc="left", fontsize=11, color=INK)
        ax.legend(frameon=False, fontsize=8.5)
        style(ax)
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, "fig2_gd18_methylation.png"), dpi=200)
        plt.close(fig)


if __name__ == "__main__":
    main()
