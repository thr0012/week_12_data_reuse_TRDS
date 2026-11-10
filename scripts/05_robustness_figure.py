#!/usr/bin/env python3
"""
Slide figure: McBirney F3 enrichment at >= 10% vs >= 50% repeat overlap, log fold scale, with DMR counts.
Usage: python3 05_robustness_figure.py results/McBirney_F3/enrichment.tsv \
         robustness/McBirney_F3_overlap50/enrichment.tsv robustness/robustness_F3_10vs50.png
"""
import sys
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
e10=pd.read_csv(sys.argv[1],sep='\t').set_index('feature')
e50=pd.read_csv(sys.argv[2],sep='\t').set_index('feature')
N=958
rows=[("Artifact-prone (union)","Artifact-prone (union)"),("Tandem repeat (TRF)","Tandem repeat (TRF)"),
      ("rmsk: Simple_repeat","Simple repeat"),("rmsk: Low_complexity","Low complexity"),("rmsk: Satellite","Satellite"),
      ("family: ERVK","ERVK"),("rmsk: LINE","LINE"),("CpG island","CpG island")]
BLUE,ORANGE,GRAY,INK,MUTED="#2a78d6","#eb6834","#a9a8a2","#0b0b0b","#52514e"
fig,ax=plt.subplots(figsize=(8.6,4.9))
lo_x,hi_x=0.08,80
for i,(f,lab) in enumerate(rows):
    for k,(e,col,dy,name) in enumerate([(e10,BLUE,-0.17,'>= 10%'),(e50,ORANGE,0.17,'>= 50%')]):
        r=e.loc[f]; y=i+dy; mu=r.null_mean
        b=lo_x
        if mu>0:
            a,b=max(r.null_lo95/mu,lo_x),min(r.null_hi95/mu,hi_x)
            ax.hlines(y,a,b,color=col,alpha=0.22,lw=7,capstyle='round')
        n=int(round(r.observed*N))
        if r.observed==0:
            ax.text(b*1.13,y,"0 DMRs",va='center',fontsize=8,color=INK)
        else:
            ax.plot(r.fold,y,'o',color=col,ms=8,mec='#fcfcfb',mew=1.5,zorder=3)
            ax.text(max(min(r.fold,hi_x),b)*1.13,y,f"{r.fold:.1f}x  (n={n})",va='center',fontsize=8,color=INK)
ax.axvline(1,color=MUTED,lw=1,ls=(0,(3,3)))
ax.set_xscale('log'); ax.set_xlim(lo_x,hi_x*4)
ax.set_xticks([0.1,0.25,0.5,1,2,5,10,25,50],["0.1","0.25","0.5","1","2","5","10","25","50"])
ax.set_yticks(range(len(rows)),[l for _,l in rows]); ax.invert_yaxis()
ax.set_xlabel("Fold enrichment vs matched background (log scale)")
for s in ("top","right"): ax.spines[s].set_visible(False)
for s in ("left","bottom"): ax.spines[s].set_color("#d6d5cf")
ax.tick_params(colors=MUTED,labelsize=9); ax.xaxis.label.set_color(MUTED)
ax.tick_params(axis='y',labelcolor=INK,length=0)
ax.grid(axis='x',color='#ecebe6',lw=0.8); ax.set_axisbelow(True)
fig.suptitle("McBirney 2017 F3 DMRs: enrichment grows when repeats must cover more of the DMR",x=0.01,ha='left',fontsize=11,color=INK,y=0.985)
fig.legend(handles=[Line2D([],[],marker='o',ls='',color=BLUE,label='Repeat covers >= 10% of DMR bp'),
                    Line2D([],[],marker='o',ls='',color=ORANGE,label='Repeat covers >= 50% of DMR bp'),
                    Line2D([],[],color=GRAY,lw=7,alpha=0.5,label='Null 95% interval')],
           frameon=False,fontsize=8.5,loc='upper left',ncol=3,bbox_to_anchor=(0.005,0.945))
fig.text(0.01,0.01,"958 DMRs (atrazine vs control lineage, p<1e-9); 1,000 matched random sets (chromosome, length, CpG density). n = DMRs overlapping.",fontsize=7.3,color=MUTED)
fig.tight_layout(rect=(0,0.03,1,0.91))
out=sys.argv[3] if len(sys.argv)>3 else 'robustness_F3_10vs50.png'
fig.savefig(out,dpi=220); print('wrote',out)
