#!/usr/bin/env python3
"""
Extract disease-specific DMR tables from Thorson et al. 2020 (PLoS ONE 15:e0239380)
supplementary PDFs (S2 lean, S3 kidney, S4 testis, S5 puberty, S6 multiple; p<1e-4) into TSV.

Repairs: CpG density glued to the gene name (1.2Ust; 0.657SK = 0.65 + 7SK, checked against
cpg_num / length); stop shown in scientific notation (1.06E+08, recomputed as start + length - 1).
Validation printed per table: parsed vs expected counts from the paper, name == chr:start,
duplicates, stop - start + 1 == length, density check, hyper/hypo counts (max_log2fc sign).

Usage: python3 extract_thorson2020.py S2.pdf S3.pdf S4.pdf S5.pdf S6.pdf
Requires: pdfplumber, pandas
"""
import pdfplumber, re, sys, pandas as pd
pat = re.compile(r'^DMR(\S+?):(\d+)\s+(\S+)\s+(\d+)\s+([0-9.]+(?:E\+\d+)?)\s+(\d+)\s+(\d+)\s+([0-9.]+E[-+]?\d+)\s+(-?[0-9.]+(?:E[-+]?\d+)?)\s+(\d+)\s+(\d+(?:\.\d+)?)\s*(.*)$')
expected = {'S2_lean':301,'S3_kidney':261,'S4_testis':693,'S5_puberty':322,'S6_multiple':336}
files = dict(zip(expected, sys.argv[1:]))
for lab, path in files.items():
    rows, bad, n = [], [], 0
    with pdfplumber.open(path) as pdf:
        for pg in pdf.pages:
            for line in (pg.extract_text() or '').splitlines():
                line = line.strip()
                if not (line.startswith('DMR') and ':' in line.split()[0]): continue
                n += 1
                m = pat.match(line)
                if not m:
                    tok = line.split(); c = tok[0][3:].split(':')[0]
                    if len(tok) > 1 and tok[0] == f"DMR{c}:{tok[1]}{c}":
                        m = pat.match(f"DMR{c}:{tok[1]} {c} " + " ".join(tok[1:]))
                (rows.append(m.groups()) if m else bad.append(line))
    df = pd.DataFrame(rows, columns=['name_chr','name_pos','chr','start','stop','length','n_sig_windows','min_p','max_log2fc','cpg_num','cpg_density_pct','gene_annotation_and_category'])
    ok = (df.name_chr == df.chr) & (df.name_pos == df.start)
    df['name'] = 'DMR' + df.name_chr + ':' + df.name_pos
    df['stop_in_pdf'] = df['stop']
    for c in ['start','length','n_sig_windows','cpg_num']: df[c] = df[c].astype(int)
    df['stop'] = df['start'] + df['length'] - 1   # PDF sometimes shows stop as 1.06E+08
    sci = df['stop_in_pdf'].str.contains('E')
    assert (df.loc[~sci,'stop_in_pdf'].astype(int) == df.loc[~sci,'stop']).all()
    df['gene_annotation_and_category'] = df['gene_annotation_and_category'].replace('', pd.NA)
    # gene names starting with a digit glued to density (0.657SK = 0.65 + 7SK): move digits back to the gene
    fixed = 0
    for i in df.index:
        dstr, calc = df.at[i, 'cpg_density_pct'], int(df.at[i, 'cpg_num']) / int(df.at[i, 'length']) * 100
        k = 0
        while abs(float(dstr[:len(dstr) - k]) - calc) > 0.006 and k < len(dstr.split('.')[-1]) - 1:
            k += 1
        if k and abs(float(dstr[:len(dstr) - k]) - calc) <= 0.006:
            df.at[i, 'gene_annotation_and_category'] = dstr[len(dstr) - k:] + df.at[i, 'gene_annotation_and_category']
            df.at[i, 'cpg_density_pct'] = dstr[:len(dstr) - k]
            fixed += 1
    print(f"   density/gene split repaired: {fixed}")
    for c in ['min_p','max_log2fc','cpg_density_pct']: df[c] = df[c].astype(float)
    df['chr'] = 'chr' + df.chr
    df = df[['name','chr','start','stop','length','n_sig_windows','min_p','max_log2fc','cpg_num','cpg_density_pct','gene_annotation_and_category']]
    df.to_csv(f'Thorson2020_{lab}_DMRs.tsv', sep='\t', index=False)
    print(f"== {lab}: DMR lines {n}, parsed {len(df)}, expected {expected[lab]}, unparsed {len(bad)}, name mismatch {(~ok).sum()}, dup {df.name.duplicated().sum()}")
    for b in bad[:4]: print('   UNPARSED:', b)
    print(f"   stops recomputed from scientific notation: {sci.sum()}")
    print(f"   stop-start+1==length: {(df.stop-df.start+1==df.length).mean():.3f} | start mod 1000: {df.start.mod(1000).value_counts().head(2).to_dict()} | max p {df.min_p.max():.2e}")
    print(f"   density check max diff {(df.cpg_num/df.length*100-df.cpg_density_pct).abs().max():.3f} | hyper (LFC>0) {(df.max_log2fc>0).sum()} / hypo {(df.max_log2fc<0).sum()} | chrs {df.chr.nunique()} incl. {sorted(set(df.chr)-{f'chr{i}' for i in range(1,21)})}")
