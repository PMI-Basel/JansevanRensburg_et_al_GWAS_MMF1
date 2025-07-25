#!/bin/bash

#SBATCH --time=1-00:00:00
#SBATCH --qos=1day
#SBATCH --mem=2g
#SBATCH --output=log.out
#SBATCH --error=log.error
#SBATCH --job-name=blastn
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=jan.waelchli@unibas.ch
#SBATCH --mail-type=ALL

#load modules

#get the taxa with blast
module load gompi/2021a
#module load BLAST+/2.11.0-gompi-2021a
module load BLAST+/2.12.0-gompi-2021a

#blast
blastn -db /scicore/data/managed/BLAST_FASTA/latest/nt -query ../../01_dadapipe/13_output/fungi/100/fungi_SEQ100.fasta \
-evalue 1e-5 \
-max_target_seqs 5 \
-num_threads $SLURM_CPUS_PER_TASK \
-outfmt "6 qseqid pident staxids" \
-out taxid.tab

#-outfmt "7 qseqid qframe sseqid sframe qstart qend sstart send evalue bitscore pident stitle" \
#-outfmt "6 qseqid pident staxids" \
#-evalue 1e-5 \
#-max_target_seqs 5 \
