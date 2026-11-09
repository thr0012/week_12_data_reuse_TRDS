# week_12_data_reuse: Do transgenerational sperm DMRs sit in regions that escape germline reprogramming?

**Week 12, Data Reuse project.** Rat, public data.

## Research question

Do published transgenerational sperm DMRs preferentially occur in genomic regions that resist DNA methylation erasure in the germline?

## Background and rationale

In mammals, DNA methylation is erased and re-established twice: in the preimplantation embryo and again in primordial germ cells (PGCs) and gonocytes. For an epimutation induced in F1 germ cells to still be detected in F3 sperm, it has to survive the germline wave at least twice. Published transgenerational DMRs are rarely evaluated against where reprogramming is actually incomplete.

The rat is the right system for a pilot because it has both:
- transgenerational sperm DMRs (Skinner lab, atrazine model), and
- a stage-resolved methylome of male rat germ cells that includes the point of lowest methylation (gestational day 18).

## Hypotheses

| | Hypothesis | Prediction |
|---|---|---|
| **H1** | Escape (biological persistence) | DMRs overlap regions that stay methylated in GD18 gonocytes more than matched background |
| **H2** | Technical or genomic artifact | DMRs are enriched in segmental duplications, satellites, simple repeats or low mappability regions, independent of GD18 status |
| **H0** | No association | DMR distribution matches matched background for both annotation sets |

H1 and H2 are not mutually exclusive. Some retrotransposon families escape reprogramming for real (IAPs in mouse), so enrichment in repeats alone does not discriminate. What discriminates is overlap with an **independent** measure of methylation in GD18 germ cells.

## Potential Datasets

| # | Role | Study | What will used |
|---|---|---|---|
| 1 | Primary DMR set | McBirney et al. 2017, *PLoS ONE* | F3 atrazine vs control lineage sperm DMRs (MeDIP-seq) | 
| 2 | Secondary DMR set | Thorson et al. 2020, *PLoS ONE* | Disease-specific F3 sperm DMRs, Tables S2 to S6 (chr, start, stop, CpG density). rn6, 1 kb windows, edgeR p < 1e-4 | 
| 3 | Escape reference (independent method) | El Omri-Charai et al. 2023, *Development* | Methyl-seq capture of 9 male germ cell stages (n = 4 per stage); GD18 gonocytes = methylation minimum | 
| 4 | Escape reference (same lab, sensitivity only) | Ben Maamar et al. 2022, *iScience* | MeDIP-seq of E13 PGCs, E16 prospermatogonia, P10 spermatogonia to sperm. Tables S1 to S6 = top 1,000 DMRs per stage comparison | 
| 5 | Genome annotation | UCSC rn6 | RepeatMasker, genomicSuperDups, gap, CpG islands, RefSeq genes | 

**OBS** Dataset 1 compares atrazine vs control lineages, which is the transgenerational exposure signature. Dataset 2 compares F3 males with vs without a disease **within** the atrazine lineage, so it answers a slightly different question (disease biomarkers). Dataset 1 is the main test; dataset 2 checks whether the pattern holds for disease DMRs.

## What will be evaluate
