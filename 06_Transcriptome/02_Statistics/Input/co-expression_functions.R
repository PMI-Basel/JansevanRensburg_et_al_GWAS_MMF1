DESeq2_analysis <- function(cts, col, shrinkage, funct) {
  
  require(DESeq2)
  
  ## merge count data with coldata
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = col,
                                design = funct)
  
  ## set factor levels to which level to compare against
  dds$conditioning <- relevel(dds$conditioning, ref = "minus")
  
  ## run analysis
  dds <- DESeq(dds, quiet = T, parallel=T)
  ## pre-filter the data to remove rows with very few reads (necessary if computation takes too long)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  ## this has to be repeated to return dds as variable in environment. otherwise DESeq doesn't recognize that DESeq() has been run and won't proceed with results
  
  # data shrinkage (optional)
  if(shrinkage){
    res.shrink <- lfcShrink(dds, coef ="conditioning_plus_vs_minus" , type="apeglm", quiet = T, parallel=T)
    DESeq2.results <- list(dds, res.shrink) 
    names(DESeq2.results) <- c("DESeqDataSet", "DESeqResults_shrink")
    (DESeq2.results)
  }
  else{
    res <- results(dds, contrast=c("conditioning","minus", "plus"), parallel=T)
    DESeq2.results <- list(dds, res) 
    names(DESeq2.results) <- c("DESeqDataSet", "DESeqResults")
    (DESeq2.results)
  }
  
}


################################################################
### function for making heatmaps and running GO ORA analysis ###
################################################################

GOanalysis <- function(x, deseq, dotplot) {
 
  require(clusterProfiler)
  require(DOSE)
  golist <- split(x$gene,x$Cluster)
  clustercomp <- compareCluster(golist,
                                OrgDb = org.At.tair.db,
                                universe = row.names(deseq[[2]]),
                                ont = "BP",
                                keyType = "TAIR",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.01,
                                fun = "enrichGO")
  clustercomp <- clusterProfiler::simplify(clustercomp, cutoff=0.7, by="p.adjust") ##simplify GO terms
  #clustercomp@compareClusterResult$Description <- gsub("biosynthetic process", "biosyn. process", clustercomp@compareClusterResult$Description )
  #clustercomp@compareClusterResult$Description <- gsub("isopentenyl diphosphate biosyn. process, methylerythritol 4-phosphate pathway", "methylerythritol 4-phosphate pathway", clustercomp@compareClusterResult$Description )
  clustercomp.df <- data.frame(clustercomp) %>% 
    dplyr::select(c(Cluster, Description, ID, geneID, p.adjust)) %>% 
    cSplit("geneID", sep="/") %>% 
    melt(id.vars=c("Cluster", "Description", "ID", "p.adjust"), value.name = "geneID", variable.name="rm", na.rm = T) %>% 
    dplyr::select(-rm)
  
  if(dotplot) {
  GO <- dotplot(clustercomp, font.size = 8, showCategory=5) +
    scale_color_gradient(low="#9b2226", high="#e9d8a6") +
    theme(legend.title = element_text(size=8),
          legend.text = element_text(size=8),
          legend.spacing.x = unit(0.2, 'cm'),
          axis.text.x = element_text(angle = 45, hjust=0)) +
    coord_flip() +
    scale_x_discrete(limits = rev) +
    scale_y_discrete(limits = rev) +
    scale_y_discrete(position = "right")

  }
  else {return(list(clustercomp.df, clustercomp))}
  
}




############################################################################
### make custom enrichment plots, functions adapted from clusterprofiler ###
############################################################################

default_labeller <- function(n) {
  function(str){
    str <- gsub("_", " ", str)
    ep_str_wrap(str, n)
  }
}

ep_str_wrap <- function(string, width) {
  x <- gregexpr(' ', string)
  vapply(seq_along(x),
         FUN = function(i) {
           y <- x[[i]]
           n <- nchar(string[i])
           len <- (c(y,n) - c(0, y)) ## length + 1
           idx <- len > width
           j <- which(!idx)
           if (length(j) && max(j) == length(len)) {
             j <- j[-length(j)]
           }
           if (length(j)) {
             idx[j] <- len[j] + len[j+1] > width
           }
           idx <- idx[-length(idx)] ## length - 1
           start <- c(1, y[idx] + 1)
           end <- c(y[idx] - 1, n)
           words <- substring(string[i], start, end)
           paste0(words, collapse="\n")
         },
         FUN.VALUE = character(1)
  )
}

dotplot_internal <- function(object, x = "geneRatio", color = "p.adjust",
                             showCategory=5, size=NULL, split = NULL,
                             font.size=8, title = "", orderBy="x",
                             decreasing=TRUE, label_format = 30) {
  
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
    if (is.null(size))
      size <- "Count"
  } else if (x == "count" || x == "Count") {
    x <- "Count"
    if (is.null(size))
      size <- "GeneRatio"
  } else if (is(x, "formula")) {
    x <- as.character(x)[2]
    if (is.null(size))
      size <- "Count"
  } else {
    ## message("invalid x, setting to 'GeneRatio' by default")
    ## x <- "GeneRatio"
    ## size <- "Count"
    if (is.null(size))
      size  <- "Count"
  }
  
  df <- fortify(object, showCategory = showCategory, split=split)
  ## already parsed in fortify
  ## df$GeneRatio <- parse_ratio(df$GeneRatio)
  
  if (orderBy !=  'x' && !orderBy %in% colnames(df)) {
    message('wrong orderBy parameter; set to default `orderBy = "x"`')
    orderBy <- "x"
  }
  
  if (orderBy == "x") {
    df <- dplyr::mutate(df, x = eval(parse(text=x)))
  }
  
  label_func <- default_labeller(label_format)
  if(is.function(label_format)) {
    label_func <- label_format
  }
  
  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description,
                           levels=rev(unique(df$Description)), ordered = T)
  
  
  ggplot(df, aes_string(x="Cluster", y="Description", size=size, color=colorBy)) +
    geom_point() +
    scale_color_continuous(low="red", high="blue", name = color,
                           guide=guide_colorbar(reverse=T)) +
    scale_y_discrete(labels = label_func) +
    ylab(NULL) + ggtitle(title) + theme_dose(font.size) +
    scale_size(range=c(3, 8))
  
  #return(df)
  
}


clustersort <- function(CL, deseq) {
  
  ### add LFC
  CL$LFC <- deseq[[2]]$log2FoldChange[match(CL$gene, row.names(deseq[[2]]))]
  ### make new Cluster order from highest median LFC to lowest
  medianexpr <- CL %>% 
    group_by(Cluster) %>% 
    summarise(medianLFC = median(LFC)) %>% 
    arrange(desc(medianLFC)) %>% 
    mutate(newCluster = c(1:nrow(.))) 
  ### replace old Cluster with newCluster
  CL$Cluster <- medianexpr$newCluster[match(CL$Cluster, medianexpr$Cluster)]
  ### relevel
  CL$Cluster <- factor(CL$Cluster, levels=c(1:nrow(medianexpr)), ordered = T)
  
  return(CL)
}




###########################################
######### make annotated heatmap ##########
###########################################


clusterheatmap <- function(CL, deseq, GOdf, GOnumber, cluster_order) {
  require(tidyr)
  require(viridis)
  require(ComplexHeatmap)
  require(org.At.tair.db)
  require(viridis)
  require(gridExtra)
  require(viridis)
  require(clusterProfiler)
  require(ggsci)
  
  x <- CL[order(CL$Cluster),]
  x <- filter(x, Cluster %in% base::unique(GOdf[[1]]$Cluster))
  mat <- pivot_wider(x, names_from = Sample, values_from = profile, id_cols = gene)
  
  # make matrix
  mat$gene <- factor(mat$gene, levels(as.factor(x$gene)), ordered = T)
  mat <- mat %>% tibble::column_to_rownames("gene")
  mat <- as.matrix(mat)
  # crop outliers
  max <- quantile(mat , .99)
  idx <- mat > max
  mat[idx] <- max
  
  min <- quantile(mat, .01)
  idx <- mat < min
  mat[idx] <- min
  
  ### make matrix for LFC 
  LFC <- data.frame(deseq[[2]]) %>% 
    tibble::rownames_to_column("gene") %>% 
    dplyr::select(gene, log2FoldChange) %>% 
    dplyr::filter(gene %in% c(row.names(mat)))
  
  LFC$Cluster <- x$Cluster[match(LFC$gene, x$gene)]
  LFC <- LFC[order(LFC$Cluster),]
  
  clus <- sort(LFC$Cluster, decreasing = F) ## cluster vector
  
  ### make matrix for padj value
  padj <- data.frame(deseq[[2]]) %>% 
    tibble::rownames_to_column("gene") %>% 
    dplyr::select(gene, padj) %>% 
    dplyr::filter(gene %in% c(row.names(mat)))
  
  padj$Cluster <- x$Cluster[match(padj$gene, x$gene)]
  padj <- padj[order(padj$Cluster),]
  
  ## make heatmap
  #conditioning = gsub("[0-9].*|.*_", "", base::colnames(mat))
  conditioning = ifelse(grepl("m[0-9]",base::colnames(mat)),"minus","plus")
  ha = HeatmapAnnotation(df = data.frame(conditioning = conditioning),
                         col = list(conditioning = c("plus" = "#000000", "minus" = "#DCDCDC")),
                         show_annotation_name = F,
                         annotation_legend_param = list(title_gp = gpar(fontsize = 8, fontface=2), labels_gp = gpar(fontsize = 8), conditioning=list(nrow=1)),
                         simple_anno_size = unit(3, "mm"))
  
  
  if (cluster_order) {
    column_dend <- as.dendrogram(hclust(dist(t(mat)))) 
  }
  else {
    column_dend <- rev(as.dendrogram(hclust(dist(t(mat))))) 
  }
 
  heatmap <- Heatmap(mat,
                     name = "expression\nprofile",
                     cluster_rows = F,
                     cluster_columns = column_dend,
                     show_row_names = F,
                     split = LFC$Cluster,
                     col = rev(rocket(256)),
                     bottom_annotation = ha,
                     column_names_gp = grid::gpar(fontsize = 8),
                     row_names_gp = grid::gpar(fontsize = 8),
                     row_title_gp = gpar(fontsize = 8),
                     width = unit(5, "cm"),
                     heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface=2), labels_gp = gpar(fontsize = 8), direction = "vertical")) +
    
    Heatmap(LFC$log2FoldChange, 
            name = "log2FoldChange", 
            width = unit(6, "mm"),
            column_names_gp = gpar(fontsize = 8),
            heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface=2), labels_gp = gpar(fontsize = 8), direction = "vertical")) +
    
    Heatmap(padj$padj,
            name = "p.adjust", 
            width = unit(6, "mm"),
            column_names_gp = gpar(fontsize = 8),
            col = cividis(256),
            heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface=2), labels_gp = gpar(fontsize = 8), direction = "vertical")) 
  
  
  ########## annotations ###########
  ## make heatmap annotations (adapted from https://jokergoo.github.io/2020/05/31/word-cloud-as-heatmap-annotation/)
  
  GO <-  GOdf[[1]][!base::duplicated(GOdf[[1]][ , c("Cluster", "ID", "geneID", "p.adjust")]),]  %>% ## filter for unique
    dplyr::select(c(Cluster, ID, Description, p.adjust)) %>% 
    base::unique() %>% 
    dplyr::select(!ID)
  
  ## select the first 5 GO terms with the lowest p.adjust for each cluster
  ## I converted the p.adjust to -log10(p.adjust) to get positive values - in case the GO terms should be displayed with different font sizes (this is   later defined with the variable gbl)
  ## I also simplified GO descriptions (e.g. biosyn. process = biosynthesis)
  
  df <- setDT(GO)[, .SD[which.minn(p.adjust, GOnumber)], by=Cluster] %>% 
    dplyr::mutate(freq=-log10(p.adjust), Description = gsub("biosyn. process", "biosynthesis",  Description),  Description = gsub("metabolic process", "metabolism", Description)) %>% 
    dplyr::rename(words=Description) %>% 
    dplyr::select(!p.adjust) %>% 
    dplyr::mutate(words = paste(words, "\t\t", sep = "")) %>% 
    base::unique()
  
  ## split df into list by Cluster
  d <- df %>% group_split(Cluster)
  d <- lapply(d, function(x) {dplyr::select(x, -Cluster)})
  names(d) <- unique(GO$Cluster)
  d <- lapply(d, function(x) {data.frame(x)})
  
  ## defining the alignment variable per cluster
  align_to <- split(seq_len(nrow(mat)), clus)
  align_to <- align_to[names(align_to) %in% names(d)]
  
  
  ## specify font size range for annotations (change fontsize_range to different values if annotations should have different font sizes according to   p.adjust values)
  fontsize_range = c(8, 8)
  gbl <-  lapply(names(align_to), function(nm) {
    kw = d[[nm]][, 1]
    freq = d[[nm]][, 2]
    fontsize = scale_fontsize(freq, rg = c(1, max(10, freq)), fs = fontsize_range)
    
    word_cloud_grob(text = kw, fontsize = fontsize, col = c("#561643", "#6C0E23", "#8F5C38", "#7180AC", "grey30"), word_space = 4)
  })
  names(gbl) <-  names(align_to)
  
  ## specify annotation margins (height + width)
  margin <- unit(8, "pt")
  gbl_h <- lapply(gbl, function(x) convertHeight(grobHeight(x), "cm") + margin)
  gbl_h <-  do.call(unit.c, gbl_h)
  
  gbl_w <- lapply(gbl, function(x) convertWidth(grobWidth(x), "cm"))
  gbl_w <- do.call(unit.c, gbl_w)
  gbl_w <- max(gbl_w) + margin
  
  
  ## function for annotation panels
  panel_fun <- function(index, nm) {
    # background
    grid.rect(gp = gpar(fill = "#EEF1EF", col = NA))
    # border
    grid.lines(c(0, 1, 1, 0), c(0, 0, 1, 1), gp = gpar(col = "#AAAAAA"),
               default.units = "npc")
    gb <- gbl[[nm]]
    # a viewport within the margins
    pushViewport(viewport(x = margin/2, y = margin/2,
                          width = grobWidth(gb), height = grobHeight(gb),
                          just = c("left", "bottom")))
    grid.draw(gb)
    popViewport()
  }
  
  annotatedhm <- heatmap + rowAnnotation(d = anno_link(align_to = align_to,
                                                       which="row",
                                                       panel_fun = panel_fun,
                                                       size = gbl_h,
                                                       gap = unit(2, "mm"),
                                                       width = gbl_w + unit(5, "mm"),
                                                       link_gp = gpar(fill = "#EEF1EF", col = "#AAAAAA"),
                                                       internal_line = F))
  
  draw(annotatedhm, heatmap_legend_side = "left", annotation_legend_side = "left", merge_legend = TRUE)
  
  
}

### overrepresentation analysis (hypergeometic test) for TF in co-expressed clusters
TF_ORA <- function(CL, TFdf, deseq) {
  
  ## merge TF families and cluster info
  TFclusters <- data.frame(merge(CL, TFdf, by="gene", all = F)) %>%
    dplyr::select(-c(Sample, profile)) %>% 
    base::unique() 
  
  ## run loop to check of overrepresented TF families in individual clusters
  clusters <- as.numeric(max(CL$Cluster))
  results <- list()
  
  for (i in 1:clusters) {
    
    ## specify universal genes
    universe <- base::rownames(deseq[[2]])
    
    ## specify gene assay with TFs
    assayed <- data.frame(universe) %>% 
      dplyr::rename(Gene = universe) %>% 
      dplyr::mutate(TFfamily = TFdf$family[base::match(Gene, TFdf$gene)]) %>% 
      na.omit()
    assayed <- split(assayed$Gene,assayed$TFfamily)
    
    ## specify significant genes (i.e. genes within a cluster)
    CL <- dplyr::filter(TFclusters, Cluster == i) %>% 
      dplyr::select(gene, family)
    siggenes <- CL$gene
    significant <- Map(function(x, y) x[x %in% y], assayed, MoreArgs=list(y=siggenes))
    
    result <- hyperg(assayed, significant, universe) %>% ## hypergeometric test
      tibble::rownames_to_column("TFamily")
    result <- apply(result,2,as.character) %>% data.frame()
    result$p <- as.numeric(result$p)
    result$padj <- p.adjust(result$p, method = "BH") ## calculate padj values
    result <- filter(result, padj < 0.05)
    results <- append(results, list(result))
    
  }
  
  names(results) <- c(1:length(results))
  results <-Filter(function(x) nrow(x) > 0, results) ## remove empty df from list
  results <- Map(cbind, results, cluster = names(results))
  results <- do.call(rbind, results)
  
  return(results)
}



### GRF normalized count plot
GRFplot <- function(deseq, DE_grfs, grfs, plottitle) {
  ncounts <- DESeq2::estimateSizeFactors(deseq[[1]]) %>% 
    counts(., normalized=TRUE) %>% 
    reshape2::melt(d.vars = "X", value.name = "Ncounts") %>% 
    mutate(conditioning = gsub("R_|\\d+", "", Var2)) %>% 
    dplyr::rename(Gene=Var1) %>% 
    filter(Gene %in% grfs$gene & Gene %in% DE_grfs$gene) %>% 
    mutate(Gene = grfs$GRF[match(Gene, grfs$gene)])
  
  ggplot(ncounts, aes(x=conditioning, y=Ncounts, col=Gene)) +
    geom_boxplot(col="grey40", show.legend = F, outlier.color = NA) +
    facet_wrap(.~Gene, scales="free") +
    theme_classic() +
    geom_jitter(width = 0.1, show.legend = F) +
    scale_color_uchicago() +
    theme(plot.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.background = element_rect(fill="grey90")) +
  labs(x="conditioning", y="normalized counts") +
    ggtitle(plottitle)
  
}

################################################################
### Function for Heatmap with legend at the bottom of heatmap ###
################################################################

clusterheatmap2 <- function(CL, deseq, GOdf, GOnumber, cluster_order) {
  require(tidyr)
  require(viridis)
  require(ComplexHeatmap)
  require(org.At.tair.db)
  require(viridis)
  require(gridExtra)
  require(viridis)
  require(clusterProfiler)
  require(ggsci)
  
  x <- CL[order(CL$Cluster),]
  x <- filter(x, Cluster %in% base::unique(GOdf[[1]]$Cluster))
  mat <- pivot_wider(x, names_from = Sample, values_from = profile, id_cols = gene)
  
  # make matrix
  mat$gene <- factor(mat$gene, levels(as.factor(x$gene)), ordered = T)
  mat <- mat %>% tibble::column_to_rownames("gene")
  mat <- as.matrix(mat)
  # crop outliers
  max <- quantile(mat , .99)
  idx <- mat > max
  mat[idx] <- max
  
  min <- quantile(mat, .01)
  idx <- mat < min
  mat[idx] <- min
  
  ### make matrix for LFC 
  LFC <- data.frame(deseq[[2]]) %>% 
    tibble::rownames_to_column("gene") %>% 
    dplyr::select(gene, log2FoldChange) %>% 
    dplyr::filter(gene %in% c(row.names(mat)))
  
  LFC$Cluster <- x$Cluster[match(LFC$gene, x$gene)]
  LFC <- LFC[order(LFC$Cluster),]
  
  clus <- sort(LFC$Cluster, decreasing = F) ## cluster vector
  
  ### make matrix for padj value
  padj <- data.frame(deseq[[2]]) %>% 
    tibble::rownames_to_column("gene") %>% 
    dplyr::select(gene, padj) %>% 
    dplyr::filter(gene %in% c(row.names(mat)))
  
  padj$Cluster <- x$Cluster[match(padj$gene, x$gene)]
  padj <- padj[order(padj$Cluster),]
  
  ## make heatmap
  #conditioning = gsub("[0-9].*|.*_", "", base::colnames(mat))
  conditioning = ifelse(grepl("m[0-9]",base::colnames(mat)),"minus","plus")
  ha = HeatmapAnnotation(df = data.frame(conditioning = conditioning),
                         col = list(conditioning = c("plus" = "gold2", "minus" = "seagreen")),
                         show_annotation_name = F,
                         annotation_legend_param = list(title_gp = gpar(fontsize = 8, fontface=2), labels_gp = gpar(fontsize = 8), conditioning=list(nrow=1)),
                         simple_anno_size = unit(3, "mm"))
  
  
  if (cluster_order) {
    column_dend <- as.dendrogram(hclust(dist(t(mat)))) 
  }
  else {
    column_dend <- rev(as.dendrogram(hclust(dist(t(mat))))) 
  }
  
  heatmap <- Heatmap(mat,
                     name = "expression\nprofile",
                     cluster_rows = F,
                     cluster_columns = column_dend,
                     show_row_names = F,
                     split = LFC$Cluster,
                     col = rev(rocket(256)),
                     bottom_annotation = ha,
                     column_names_gp = grid::gpar(fontsize = 8),
                     row_names_gp = grid::gpar(fontsize = 8),
                     row_title_gp = gpar(fontsize = 8),
                     width = unit(5, "cm"),
                     heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface=2), labels_gp = gpar(fontsize = 8), direction = "horizontal")) +
    
    Heatmap(LFC$log2FoldChange, 
            name = "log2FoldChange", 
            width = unit(6, "mm"),
            column_names_gp = gpar(fontsize = 8),
            heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface=2), labels_gp = gpar(fontsize = 8), direction = "horizontal")) +
    
    Heatmap(padj$padj,
            name = "p.adjust", 
            width = unit(6, "mm"),
            column_names_gp = gpar(fontsize = 8),
            col = cividis(256),
            heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface=2), labels_gp = gpar(fontsize = 8), direction = "horizontal")) 
  
  
  ########## annotations ###########
  ## make heatmap annotations (adapted from https://jokergoo.github.io/2020/05/31/word-cloud-as-heatmap-annotation/)
  
  GO <-  GOdf[[1]][!base::duplicated(GOdf[[1]][ , c("Cluster", "ID", "geneID", "p.adjust")]),]  %>% ## filter for unique
    dplyr::select(c(Cluster, ID, Description, p.adjust)) %>% 
    base::unique() %>% 
    dplyr::select(!ID)
  
  ## select the first 5 GO terms with the lowest p.adjust for each cluster
  ## I converted the p.adjust to -log10(p.adjust) to get positive values - in case the GO terms should be displayed with different font sizes (this is   later defined with the variable gbl)
  ## I also simplified GO descriptions (e.g. biosyn. process = biosynthesis)
  
  df <- setDT(GO)[, .SD[which.minn(p.adjust, GOnumber)], by=Cluster] %>% 
    dplyr::mutate(freq=-log10(p.adjust), Description = gsub("biosyn. process", "biosynthesis",  Description),  Description = gsub("metabolic process", "metabolism", Description)) %>% 
    dplyr::rename(words=Description) %>% 
    dplyr::select(!p.adjust) %>% 
    dplyr::mutate(words = paste(words, "\t\t", sep = "")) %>% 
    base::unique()
  
  ## split df into list by Cluster
  d <- df %>% group_split(Cluster)
  d <- lapply(d, function(x) {dplyr::select(x, -Cluster)})
  names(d) <- unique(GO$Cluster)
  d <- lapply(d, function(x) {data.frame(x)})
  
  ## defining the alignment variable per cluster
  align_to <- split(seq_len(nrow(mat)), clus)
  align_to <- align_to[names(align_to) %in% names(d)]
  
  
  ## specify font size range for annotations (change fontsize_range to different values if annotations should have different font sizes according to   p.adjust values)
  fontsize_range = c(8, 8)
  gbl <-  lapply(names(align_to), function(nm) {
    kw = d[[nm]][, 1]
    freq = d[[nm]][, 2]
    fontsize = scale_fontsize(freq, rg = c(1, max(10, freq)), fs = fontsize_range)
    
    word_cloud_grob(text = kw, fontsize = fontsize, col = c("#561643", "#6C0E23", "#8F5C38", "#7180AC", "grey30"), word_space = 4)
  })
  names(gbl) <-  names(align_to)
  
  ## specify annotation margins (height + width)
  margin <- unit(8, "pt")
  gbl_h <- lapply(gbl, function(x) convertHeight(grobHeight(x), "cm") + margin)
  gbl_h <-  do.call(unit.c, gbl_h)
  
  gbl_w <- lapply(gbl, function(x) convertWidth(grobWidth(x), "cm"))
  gbl_w <- do.call(unit.c, gbl_w)
  gbl_w <- max(gbl_w) + margin
  
  
  ## function for annotation panels
  panel_fun <- function(index, nm) {
    # background
    grid.rect(gp = gpar(fill = "#EEF1EF", col = NA))
    # border
    grid.lines(c(0, 1, 1, 0), c(0, 0, 1, 1), gp = gpar(col = "#AAAAAA"),
               default.units = "npc")
    gb <- gbl[[nm]]
    # a viewport within the margins
    pushViewport(viewport(x = margin/2, y = margin/2,
                          width = grobWidth(gb), height = grobHeight(gb),
                          just = c("left", "bottom")))
    grid.draw(gb)
    popViewport()
  }
  
  annotatedhm <- heatmap + rowAnnotation(d = anno_link(align_to = align_to,
                                                       which="row",
                                                       panel_fun = panel_fun,
                                                       size = gbl_h,
                                                       gap = unit(2, "mm"),
                                                       width = gbl_w + unit(5, "mm"),
                                                       link_gp = gpar(fill = "#EEF1EF", col = "#AAAAAA"),
                                                       internal_line = F))
  
  draw(annotatedhm, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE)
  
  
}

