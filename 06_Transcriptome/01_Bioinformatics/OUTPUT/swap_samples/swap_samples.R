## Run from this directory:  Rscript swap_samples.R

RPG_original <- read.table("ReadsPerGene_all_original.txt")

#swap WTm2 and WTp2
RPG_swapped <- RPG_original
colnames(RPG_swapped)[colnames(RPG_original)=="WTm2"] <- "WTp2"
colnames(RPG_swapped)[colnames(RPG_original)=="WTp2"] <- "WTm2"
RPG_swapped <- RPG_swapped[,colnames(RPG_original)]

write.table(RPG_swapped, "ReadsPerGene_all_swapped.txt")
