#!/usr/bin/env python3
"""
Extract DMR tables from the McBirney et al. 2017 (PLoS ONE 12:e0184306) supplementary PDFs
(Tables S3 F1 p<1e-6, S4 F2 p<1e-5, S5 F3 p<1e-9) into TSV.

Tables have no stop column: stop = start + length - 1 (1-based, 100 bp MEDIPS windows).
Repairs rows where the PDF glued the DMR name to the chromosome column
(DMR10:101600401 + 10 -> DMR10:10160040110).
Validation printed per table: rows parsed vs DMR lines, name == chr:start, duplicates,
start mod 100 == 1, CpG density == cpg_num / length.

Usage: python3 extract_mcbirney2017.py F1.pdf F2.pdf F3.pdf   (writes McBirney2017_F*_DMRs.tsv)
Requires: pdfplumber, pandas
"""
import pdfplumber, re, sys, pandas as pd
pat = re.compile(r'^DMR(\S+?):(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([0-9.]+E[-+]?\d+)\s+(\d+)\s+([0-9.]+)(?:\s+(\S+))?(?:\s+(.+))?$')
for path, gen in [(sys.argv[1],'F1'),(sys.argv[2],'F2'),(sys.argv[3],'F3')]:
    rows, bad, n_dmr_like = [], [], 0
    with pdfplumber.open(path) as pdf:
        for pg in pdf.pages:
            for line in (pg.extract_text() or '').splitlines():
                line=line.strip()
                if line.startswith('DMR') and ':' in line.split()[0]:
                    n_dmr_like += 1
                    m = pat.match(line)
                    if not m:  # PDF glued the name to the chr column: DMR10:101600401 + 10 -> DMR10:10160040110
                        tok = line.split()
                        c = tok[0][3:].split(':')[0]
                        if len(tok) > 1 and tok[0] == f"DMR{c}:{tok[1]}{c}":
                            line = f"DMR{c}:{tok[1]} {c} " + " ".join(tok[1:])
                            m = pat.match(line)
                    if m: rows.append(m.groups())
                    else: bad.append(line)
    df = pd.DataFrame(rows, columns=['name_chr','name_pos','chr','start','length','n_sig_windows','min_p','cpg_num','cpg_density_pct','gene','gene_category'])
    df['name']='DMR'+df.name_chr+':'+df.name_pos
    for c in ['start','length','n_sig_windows','cpg_num']: df[c]=df[c].astype(int)
    df['min_p']=df.min_p.astype(float); df['cpg_density_pct']=df.cpg_density_pct.astype(float)
    df['chr']='chr'+df.chr.astype(str)
    df['stop']=df.start+df.length-1
    consistent = (df.name_chr.radd('chr')==df.chr) & (df.name_pos.astype(int)==df.start)
    df=df[['name','chr','start','stop','length','n_sig_windows','min_p','cpg_num','cpg_density_pct','gene','gene_category']]
    out=f'McBirney2017_{gen}_DMRs.tsv'; df.to_csv(out,sep='\t',index=False)
    print(f'== {gen}: lines starting with DMR: {n_dmr_like}, parsed: {len(df)}, unparsed: {len(bad)}, name/coord mismatches: {(~consistent).sum()}, duplicates: {df.name.duplicated().sum()}')
    for b in bad[:5]: print('   UNPARSED:', b)
    print('   max min_p:', df.min_p.max(), '| start mod 100:', df.start.mod(100).value_counts().head(3).to_dict(), '| length mod 100:', df.length.mod(100).value_counts().head(2).to_dict())
    print('   chromosomes:', df.chr.value_counts().to_dict())
    print('   length median/max:', df.length.median(), df.length.max(), '| density/100bp check (cpg_num/length*100 vs column) max abs diff:', round((df.cpg_num/df.length*100-df.cpg_density_pct).abs().max(),2))
