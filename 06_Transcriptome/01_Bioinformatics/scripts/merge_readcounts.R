samples <- list.files("04_mapping/")

files <- list.files(path= paste0("04_mapping/", samples), pattern = "*ReadsPerGene.out.tab", full.names = TRUE)
counts.files <- lapply( files, read.table, skip = 4 )

counts <- as.data.frame( sapply( counts.files, function(x) x[ , 4 ] ) ) 

colnames(counts) <- samples
row.names(counts) <- counts.files[[1]]$V1

write.table(counts, "07_counts/ReadsPerGene_all.txt", row.names = T, col.names = T, sep = "\t", quote = F)


