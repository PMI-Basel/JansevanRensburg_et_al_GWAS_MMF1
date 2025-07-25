libraries <- function(){
  #BiocManager
  if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
  if (!require(phyloseq)){BiocManager::install("phyloseq")}; library(phyloseq)
  if (!require(vegan)){BiocManager::install("vegan")}; library(vegan)
  if (!require(Biostrings)){BiocManager::install("Biostrings")}; library(Biostrings)
  if (!require(DECIPHER)){BiocManager::install("DECIPHER")}; library(DECIPHER)

  #data manipulation
  if (!require(readxl)){install.packages("readxl")}; library(readxl)
  if (!require(tidyr)){install.packages("tidyr")}; library(tidyr)
  if (!require(dplyr)){install.packages("dplyr")}; library(dplyr)
  if (!require(tibble)){install.packages("tibble")}; library(tibble)

  #figures and tables
  if (!require(ggplot2)){install.packages("ggplot2")}; library(ggplot2)
  if (!require(cowplot)){install.packages("cowplot")}; library(cowplot)
  if (!require(pander)){install.packages("pander")}; library(pander)
  if (!require(RColorBrewer)){install.packages("RColorBrewer")}; library(RColorBrewer)
  if (!require(ggh4x)){install.packages("ggh4x")}; library(ggh4x)
  if (!require(ggpubr)){install.packages("ggpubr")}; library(ggpubr)
  if (!require(ggnewscale)){install.packages("ggnewscale")}; library(ggnewscale)
}
