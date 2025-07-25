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
  if (!require(forcats)){install.packages("forcats")}; library(forcats)
  if (!require(stringr)){install.packages("stringr")}; library(stringr)

  #figures and tables
  if (!require(ggplot2)){install.packages("ggplot2")}; library(ggplot2)
  if (!require(ggpattern)){install.packages("ggpattern")}; library(ggpattern)
  if (!require(cowplot)){install.packages("cowplot")}; library(cowplot)
  if (!require(pander)){install.packages("pander")}; library(pander)
  if (!require(RColorBrewer)){install.packages("RColorBrewer")}; library(RColorBrewer)
  if (!require(ggh4x)){install.packages("ggh4x")}; library(ggh4x)
  
  #stats
  if (!require(hilldiv)){install.packages("hilldiv")}; library(hilldiv)
  if (!require(emmeans)){install.packages("emmeans")}; library(emmeans)
  if (!require(multcomp)){install.packages("multcomp")}; library(multcomp)
  if (!require(ALDEx2)){BiocManager::install("ALDEx2")}; library(ALDEx2)
  if (!require(ANCOMBC)){BiocManager::install("ANCOMBC")}; library(ANCOMBC)
  if (!require(Maaslin2)){BiocManager::install("Maaslin2")}; library(Maaslin2)
  if (!require(metagenomeSeq)){BiocManager::install("metagenomeSeq")}; library(metagenomeSeq)
  if (!require(sjstats)){install.packages("sjstats")}; library(sjstats)
  }
