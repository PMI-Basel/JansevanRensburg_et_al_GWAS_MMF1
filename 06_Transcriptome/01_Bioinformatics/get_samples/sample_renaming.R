## Build get_samples.tab: sequencing-facility filename -> analysis sample name.
## Run from this directory:  Rscript sample_renaming.R

## directory holding the raw fastq files (override with RNASEQ_RAW_DIR)
RAW_DIR <- Sys.getenv("RNASEQ_RAW_DIR", unset = "/path/to/raw_data/RNAseq_Henry")

#import sample table
samples <- read.table("all_RNASeq_samples.txt")
colnames(samples) <- "path"

#add full path
samples$path <- file.path(RAW_DIR, samples$path)

#remove path
library(tidyverse)
samples$ori_name <- samples$path %>% basename()

#subset to samples from Henry
samples <- samples %>% filter(grepl("(WT|mmf1)", ori_name))

#get name and lane
samples$name <- str_extract(samples$ori_name, "(WT|mmf1)(p|m)[0-9]+")
samples$lane <- str_extract(samples$ori_name, "L00(1|2)")
samples$dir <- str_extract(samples$ori_name, "R(1|2)")
samples$new_name <- paste0(samples$name, ".", samples$lane, "_", samples$dir, ".fastq.gz")

#make table
to_write <- samples %>% dplyr::select(path, new_name)
write_tsv(to_write, "get_samples.tab", col_names = F, eol = "\n")
