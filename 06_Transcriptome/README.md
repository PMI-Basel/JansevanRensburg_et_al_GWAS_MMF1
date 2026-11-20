# RNA-seq analysis — Janse van Rensburg et al.

Two stages, run in order:

- `01_Bioinformatics/` — read processing and quantification (Snakemake: STAR, FastQC,
  MultiQC, SAMtools, Qualimap), producing the gene-level count matrix
  `01_Bioinformatics/OUTPUT/ReadsPerGene_all.txt`. One sample pair (`WTm2`/`WTp2`) was
  switched during the experiment and is corrected here — see
  `01_Bioinformatics/OUTPUT/swap_samples/`.
- `02_Statistics/` — differential expression, co-expression clustering and GO enrichment
  in R, taking the corrected count matrix as input.
