#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --qos=6hours
#SBATCH --mem=2g
#SBATCH --output=RNASeq.out
#SBATCH --error=RNASeq.error
#SBATCH --job-name=RNASeq
#SBATCH --cpus-per-task=1
## set your own address, or delete this line and the next
#SBATCH --mail-user=USER@INSTITUTION
#SBATCH --mail-type=ALL

#create folder for log files
mkdir -p log/jobs 2> /dev/null

#activate minconda3 and snakemake environment
## path to your own conda installation
source ${CONDA_BASE:-$HOME/miniconda3}/bin/activate
conda activate snakemake

#unlock snakemake folder (necessary when last Snakemake job failed)
snakemake --unlock

#start snakemake
snakemake \
--jobs 100 \
--cluster-config cluster_config.yaml \
--cluster "sbatch \
--time={cluster.time} \
--qos={cluster.qos} \
--cpus-per-task={cluster.cpus-per-task} \
--mem={cluster.mem-per-cpu} \
--output={cluster.output} \
--error={cluster.error}" \
-p \
--jobname RNASeq_{jobid} \
--latency-wait 60

#mv log files
mv RNASeq.error log
mv RNASeq.out log
