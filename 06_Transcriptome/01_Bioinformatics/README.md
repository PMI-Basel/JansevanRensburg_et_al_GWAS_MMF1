# RNA-seq bioinformatics pipeline

Read processing and quantification for the root transcriptome of *Arabidopsis thaliana*
Col-0 (WT) and *mmf1-2* grown on BXplus and BXminus soil microbiomes
(Janse van Rensburg et al.). The downstream statistical analysis is not here — it lives in
the manuscript code repository, document `06_Transcriptome.qmd`, which takes
`OUTPUT/ReadsPerGene_all.txt` as its input.

## What the pipeline does

Snakemake workflow (`Snakefile`), one rule per step, submitted to SLURM:

1. `Index_REFgenome` — STAR 2.7.9a genome index from the TAIR10 assembly and the
   Araport11 annotation (`--sjdbOverhang 50`, `--genomeSAindexNbases 12`).
2. `FastQC` / `MultiQC` — per-fastq quality control (FastQC 0.11.8, MultiQC 1.7).
3. `alignment` — STAR, two lanes per sample given as comma-separated read lists,
   `--quantMode GeneCounts`, `--alignIntronMin 50`, `--alignIntronMax 10223`.
4. `bam_sort_byname` — SAMtools 1.9, name-sorted BAM for Qualimap.
5. `alignment_quality` — Qualimap 2.2.1 `rnaseq`, strand-specific-reverse, paired.
6. `alignment_quality_report` — renders `scripts/alignment_quality.Rmd` into
   `alignment_quality.report.html` (STAR mapping rates plus genomic origin of reads).
7. `make_countmatrix` — `scripts/merge_readcounts.R` collects column 4 of each
   `*ReadsPerGene.out.tab` into `07_counts/ReadsPerGene_all.txt`.

20 samples: WT and *mmf1-2*, BXminus (`m`) and BXplus (`p`), five biological
replicates each, two lanes and paired reads per sample.

## The sample switch (read this before using the counts)

`WTm2` and `WTp2` were exchanged at some point between sowing and sequencing. The
analysis therefore uses a corrected count matrix in which those two column labels are
swapped back; `OUTPUT/swap_samples/` holds the whole operation:

| file | content |
|---|---|
| `ReadsPerGene_all_original.txt` | matrix as the pipeline produced it |
| `swap_samples.R` | exchanges the `WTm2` and `WTp2` column names |
| `ReadsPerGene_all_swapped.txt` | corrected matrix — identical to `OUTPUT/ReadsPerGene_all.txt` |
| `README.txt` | reason for the correction |

Only those two columns differ between the two matrices, and they are exchanged exactly
(`original$WTm2 == swapped$WTp2` and `original$WTp2 == swapped$WTm2`); all 18 other
samples and all 36,799 gene rows are untouched. `ReadsPerGene_all_swapped.txt` is
byte-identical (MD5 `cde4e9685ce3640662dc408fb935f2b8`) to the file the statistical
analysis reads, so the published results are based on the corrected labelling.

Both matrices are kept deliberately: a reader can reproduce the correction, or redo the
analysis without it.

## Reference genome

Not included (1.5 GB). Download before running:

- Assembly: TAIR10 chromosome fasta (`TAIR10.genome.fa`, TAIR / Araport)
- Annotation: `Araport11_GTF_genes_transposons.Mar202021.gtf`

Place both under `01_REFgenome/` with the names given in `config.yaml`; the
`Index_REFgenome` rule builds `01_REFgenome/REFindex/` from them.

## Running it

```bash
# 1. fetch and rename the raw reads (writes 00_rawdat/ and samples.txt)
cd get_samples && Rscript sample_renaming.R && bash get_samples.sh && cd ..

# 2. run the workflow on the cluster
sbatch job_exec.sh
```

`job_exec.sh` and `cluster_config.yaml` are written for a SLURM cluster with
environment modules and a conda environment named `snakemake`; adapt the account,
partition, module names and the conda activation path to your site.

## Raw data

Raw fastq files are deposited at European Nucleotide Archive (http://www.ebi.ac.uk/ena) PRJEB109937. `get_samples/get_samples.tab`
records which sequencing-facility file became which sample name, so the deposited files
can be matched to the sample names used throughout the analysis.

## Software versions

STAR 2.7.9a · FastQC 0.11.8 · MultiQC 1.7 · SAMtools 1.9 · Qualimap 2.2.1 ·
R 4.1.0 · Snakemake (conda environment)
