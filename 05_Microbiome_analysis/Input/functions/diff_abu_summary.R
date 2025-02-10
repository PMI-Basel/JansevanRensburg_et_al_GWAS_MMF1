
# function to create a summary data frame for the differential abundance analysis

fun_diff_abu_summary <- function(summ, means, SEs, TAX){
  
  summary <- pivot_longer(means, !(ASV), names_to = "conditioning", values_to = "means")
  
  summary$se <- pivot_longer(SEs, !(ASV), names_to = "conditioning", values_to = "se") %>% pull(se)
  summary$ASV <- factor(summary$ASV, levels = unique(summary$ASV)) # keep the sorted order of ASVs
  summary$conditioning <- as.factor(summary$conditioning)
  
  # get score
  summary$score <- 0
  for (p in unique(summary$conditioning)) {
    for (a in unique(summary$ASV)) {
      score <- summ %>%
        filter(ASV == a) %>%
        pull(score)
      if (length(score) == 1) {
        summary$score[(summary$conditioning == p & summary$ASV == a)] <- score
      }
    }
  }
  
  
  # add if different abundant (score > 2)
  summary$diff_abundant <- ifelse(summary$score > 2, "differently_abundant", "unchanged")
  summary$diff_abundant <- as.factor(summary$diff_abundant)
  
  
  # add abundance in BX+ compared to BX-
  summary$abu_BX <- "unchanged"
  for (p in unique(summary$conditioning)) {
    for (a in unique(summary$ASV)) {
      mean_BXm <- summary$means[(summary$ASV == a) & (summary$conditioning == "BXm")]
      mean_BXp <- summary$means[(summary$ASV == a) & (summary$conditioning == "BXp")]
      abu_diff <- summary$diff_abundant[(summary$ASV == a) & (summary$conditioning == p)]
      if (abu_diff == "differently_abundant") {
        abu_BX_value <- ifelse(mean_BXp > mean_BXm, "higher in BX+", "lower in BX+")
        summary$abu_BX[(summary$ASV == a) & (summary$conditioning == p)] <- abu_BX_value
      }
    }
  }
  
  # add phyla and class
  
  TAX_temp <- TAX
  TAX_temp$ASV_ID <- rownames(TAX_temp)
  
  
  summary$phylum <- NA
  summary$class <- NA
  summary$order <- NA
  summary$family <- NA
  summary$genus <- NA
  
  for (a in unique(summary$ASV)) {
    if (a %in% TAX_temp$ASV_ID) {
      summary$phylum[summary$ASV == a] <- as.character(TAX_temp$Phylum[TAX_temp$ASV_ID == a])
      summary$class[summary$ASV == a] <- as.character(TAX_temp$Class[TAX_temp$ASV_ID == a])
      summary$order[summary$ASV == a] <- as.character(TAX_temp$Order[TAX_temp$ASV_ID == a])
      summary$family[summary$ASV == a] <- as.character(TAX_temp$Family[TAX_temp$ASV_ID == a])
      summary$genus[summary$ASV == a] <- as.character(TAX_temp$Genus[TAX_temp$ASV_ID == a])
    }
  }
  
  return(summary)
  
}