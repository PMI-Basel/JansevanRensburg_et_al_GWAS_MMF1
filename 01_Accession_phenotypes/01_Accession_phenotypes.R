## ----r------------------------------------------------------------------------
cat("\\small")
knitr::knit_hooks$set(purl = knitr::hook_purl)
knitr::opts_chunk$set(echo = TRUE)

##clear the objects from memory
rm(list=ls())
# Purpose: Clears all objects from the R environment to 
#ensure a clean workspace.

## ----r------------------------------------------------------------------------
#| label: Libraries
#| message: false
#| warning: false
#| paged-print: false
library(tidyverse)
#devtools::install_github("Gregor-Mendel-Institute/gwaR"); 
#https://github.com/Gregor-Mendel-Institute/gwaR
library(gwaR)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DescTools)
library(ggpubr)
library(MESS)
library(cowplot)
library(knitr)
library(ggdendro)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(openxlsx)
# Purpose: Loads the necessary R packages for data manipulation, visualization, 
#and GWAS analysis.

## ----r tidy=TRUE--------------------------------------------------------------
library("readxl")
tot_plant_specific_snps <- read_excel("Input/Dataset1.xlsx", sheet = "SNP_table")

# Change the orde of the table to show highest log10_p values first
tot_plant_specific_snps <- tot_plant_specific_snps[order(tot_plant_specific_snps$log10_p, decreasing = TRUE),]
kable(tot_plant_specific_snps[1:3, ], caption = "SNP table.")

#Import plant areas of accessions
plant_region_areas <- read_excel("Input/Dataset1.xlsx", sheet = "plant_region_area")

# Purpose: Imports phenotype data files into R as data frames.


## ----r------------------------------------------------------------------------
# Read phenotypes from excel file
Col0_plant_region_areas <- read_excel("Input/Dataset1.xlsx", sheet = "Col0_growth")

# --- Plot growth over time (area vs. age) ---

Col0_growth <- ggplot(Col0_plant_region_areas, aes(x = age, y = plant_region_area, colour = BX, fill = BX)) +
  # Line for group mean
  stat_summary(fun = mean, geom = "line", size = 1, alpha = 0.8, show.legend = FALSE) +
  
  # Error bars for SE
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  
  # Points for group mean
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2.5, alpha = 0.4, color = "black") +
  
  # Axis labels
  labs(
    title = "Daily rosette area: Col-0",
    x = "Days after germination",
    y = expression("Rosette area (10"^4*")")
  ) +
  
  # Rescale Y axis to match IP_Hue_3 formatting
  scale_y_continuous(labels = function(x) x / 10000) +
  
  # Optional: Custom breaks if you want consistency
  # scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  
  # Color & fill definitions
  scale_colour_manual(
    values = c("seagreen", "gold2"),
    guide = "none"  # hide color legend
  ) +
  scale_fill_manual(
    values = c("seagreen", "gold2"),
    name = "Microbiome",
    labels = c("BX-", "BX+"),
    guide = guide_legend(
      override.aes = list(shape = 21, size = 2.5, alpha = 0.4)
    )
  ) +
  
  # Theme formatting
  theme_bw() +
  theme(
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    plot.title = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.title.x = element_text(size = 7),
    legend.position = c(0.2, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "black")
  )

Col0_growth

## ----r------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(stringr)

# Ensure age_bx is constructed consistently
Col0_plant_region_areas <- Col0_plant_region_areas %>%
  mutate(age_bx = paste0("DAG", age, "_", BX))

# Extract numeric age and make ordering with BX- before BX+
bx_levels <- Col0_plant_region_areas %>%
  mutate(
    age_num = as.numeric(str_extract(age_bx, "\\d+")),
    bx_order = ifelse(BX == "minus", 0, 1)
  ) %>%
  arrange(age_num, bx_order) %>%
  distinct(age_bx) %>%
  pull(age_bx)

# Apply the custom order
Col0_plant_region_areas <- Col0_plant_region_areas %>%
  mutate(age_bx = factor(age_bx, levels = bx_levels))

# Plot
Col0_plant_growth <- ggplot(Col0_plant_region_areas, aes(x = factor(age), y = plant_region_area, fill = BX)) +
  geom_boxplot(
    position = position_dodge(width = 1),
    outlier.shape = NA,
    width = 0.7,
    alpha = 0.6,
    color = "black"
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 1),
    shape = 21,
    size = 1.0,
    alpha = 0.8,
    color = "black"
  ) +
  scale_fill_manual(
    values = c("minus" = "seagreen", "plus" = "gold2"),
    name = "Microbiome",
    labels = c("BX-", "BX+"),
    guide = guide_legend(override.aes = list(shape = 21, size = 1, alpha = 0.4))
  ) +
  scale_y_continuous(
    labels = function(x) x / 10000,
    limits = c(0, 1500000)
  ) +
  labs(
    title = "Rosette Area of Col-0 Over Time (BX Comparison)",
    x = "Days After Germination",
    y = expression("Rosette area (10"^4*")")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 7),
    plot.title = element_text(size = 7),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.position = c(0.10, 0.95),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "black")
  )

Col0_plant_growth

## ----r------------------------------------------------------------------------
library(dplyr)
library(broom)

# Run t-tests per age
t_test_results <- Col0_plant_region_areas %>%
  group_by(age) %>%
  filter(BX %in% c("minus", "plus")) %>%
  summarise(
    t_test = list(t.test(plant_region_area ~ BX))
  ) %>%
  mutate(tidy_result = purrr::map(t_test, broom::tidy)) %>%
  unnest(tidy_result) %>%
  select(
    age,
    estimate1,  # mean of "minus"
    estimate2,  # mean of "plus"
    statistic,
    p.value,
    conf.low,
    conf.high
  ) %>%
  mutate(p.value = signif(p.value, 3))

# View results
print(t_test_results)

# Optional: Save to Excel
# openxlsx::write.xlsx(t_test_results, "Col0_BX_t_test_results.xlsx", rowNames = FALSE)

## ----r------------------------------------------------------------------------
library(dplyr)
library(purrr)
library(tidyr)
library(broom)

# Make sure purrr::map is used instead of maps::map
anova_results <- Col0_plant_region_areas %>%
  group_by(age) %>%
  nest() %>%
  mutate(
    anova_model = purrr::map(data, ~ aov(plant_region_area ~ BX, data = .x)),
    anova_tidy  = purrr::map(anova_model, broom::tidy),
    tukey_hsd   = purrr::map(anova_model, ~ TukeyHSD(.x, "BX") %>% broom::tidy())
  ) %>%
  select(age, anova_tidy, tukey_hsd) %>%
  unnest(cols = c(anova_tidy)) %>%
  filter(term == "BX")

# Print ANOVA summary results
print(anova_results)

## ----r------------------------------------------------------------------------
# Tukey HSD post-hoc summary
tukey_results <- Col0_plant_region_areas %>%
  group_by(age) %>%
  nest() %>%
  mutate(
    anova_model = purrr::map(data, ~ aov(plant_region_area ~ BX, data = .x)),
    tukey = purrr::map(anova_model, ~ TukeyHSD(.x, "BX") %>% broom::tidy())
  ) %>%
  unnest(tukey)

# View Tukey results
print(tukey_results)

## ----r------------------------------------------------------------------------
library(lme4)
library(lmerTest)   # for p-values
library(emmeans)    # for post-hoc

# --- Add unique plant ID ---
Col0_plant_region_areas <- Col0_plant_region_areas %>%
  mutate(
    plant_id = paste(line, BX, number, sep = "_"),
    BX = factor(BX, levels = c("minus", "plus")),
    age = as.numeric(age)  # or factor(age)
  )

## ----r------------------------------------------------------------------------
# Fit linear mixed-effects model
model_lmm <- lmer(plant_region_area ~ age * BX + (1 | plant_id), data = Col0_plant_region_areas)

# Summary with p-values
summary(model_lmm)


## ----r------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(stringr)


# --- Calculate mean per age and BX (match new x-axis structure) ---
mean_df <- Col0_plant_region_areas %>%
  group_by(age, BX) %>%
  summarise(mean_area = mean(plant_region_area), .groups = "drop")

# --- Create connecting lines between BX-minus and BX-plus means for each age ---
line_df <- mean_df %>%
  pivot_wider(names_from = BX, values_from = mean_area) %>%
  filter(!is.na(minus) & !is.na(plus))  # Remove any incomplete rows

# Plot with means and connecting lines
Col0_plant_growth1 <- Col0_plant_growth +
  # Add line segments between means
  geom_segment(
    data = line_df,
    aes(
      x = as.numeric(factor(age)) - 0.25,  # position for BX-minus
      xend = as.numeric(factor(age)) + 0.25,  # position for BX-plus
      y = minus,
      yend = plus
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 0.5
  ) +
  # Add triangle mean points
  geom_point(
    data = mean_df,
    aes(
      x = as.numeric(factor(age)) + ifelse(BX == "minus", -0.25, 0.25),
      y = mean_area,
      fill = BX
    ),
    shape = 24,
    size = 1.3,
    color = "black",
    inherit.aes = FALSE
  )

Col0_plant_growth1

ggsave("Col0_growth_feedback.pdf", plot = Col0_plant_growth1, width = 12, height = 6, units = "cm", path = "Output/Figures/")

## ----r tidy=TRUE--------------------------------------------------------------
# Extract plant area
temp1 <- tot_plant_specific_snps %>%
  {gwaR::phenotype_by_snp(plant_region_areas %>% mutate(ACC_ID = genotype_id),
                          "plant_region_area",
                          gwas_table = .,
                          specific = "bx",
                          SNPrank = 7,
                          plot = FALSE, 
                          SNPmatrix = "Input/SNPmatrix_full.fst", # Change this path!
                          nobees = TRUE)}

# Extract dag
temp2 <- tot_plant_specific_snps %>% 
  {gwaR::phenotype_by_snp(plant_region_areas %>% mutate(ACC_ID = genotype_id),
                          "dag",
                          gwas_table = .,
                          specific = "bx",
                          SNPrank = 7,
                          plot = FALSE, 
                          SNPmatrix = "Input/SNPmatrix_full.fst", # Change this path!
                          nobees = TRUE)}

# Combine to one df
plant_region_areas_SNP <- data.frame(
  ACC_ID = temp1$ACC_ID, 
  bx = temp1$bx, 
  plant_region_area = temp1$phenotype_value,
  dag = as.integer(temp2$phenotype_value),  # Convert dag to integers
  SNP = temp1$SNP
)

# Remove temporary objects
rm(temp1, temp2)

# Purpose: Extracts and combines phenotype data based on SNP rank 7 and prepares it for further analysis.
kable(plant_region_areas_SNP[1:3, ], caption = "Phenotype table for SNP rank 7.")

## ----r------------------------------------------------------------------------
# Filter the plant region areas for DAG ≤ 20
plant_region_areas_DAG20 <- subset(plant_region_areas, subset = dag <= 21)

# Calculate mean plant region area for each combination of bx, dag, and genotype_id
plant_region_areas_DAG20 <- plant_region_areas_DAG20 %>%
  group_by(bx, dag, genotype_id) %>%
  mutate(mean = mean(plant_region_area, na.rm = TRUE)) %>%
  ungroup()

# Remove duplicates and unnecessary columns
plant_region_areas_DAG20 <- plant_region_areas_DAG20 %>%
  select(-c(1, 6)) %>%
  distinct()

# Pivot data to spread 'bx' treatment values into separate columns and calculate the difference
plant_region_areas_DAG20 <- plant_region_areas_DAG20 %>%
  spread(key = "bx", value = "mean") %>%
  mutate(area_diff = `+` - `-`) %>%
  arrange(desc(area_diff))

# Rank the entries by area difference
plant_region_areas_DAG20 <- plant_region_areas_DAG20 %>%
  mutate(rank = row_number(desc(area_diff)))

# Purpose: Ranks the genotypes by the calculated area difference.
kable(plant_region_areas_DAG20[1:3, ], caption = "Phenotype table for DAG 20.")

## ----r------------------------------------------------------------------------
# Step 1: Create a dataframe with preferences based on "dag" 15
preferences <- plant_region_areas_DAG20 %>%
  filter(dag == 10) %>%
  mutate(preference = case_when(
    area_diff >= 0 ~ "plus",
    TRUE ~ "minus"
  )) %>%
  select(genotype, preference)
# Purpose: Assigns a preference ("plus" or "minus") to each genotype based on the area difference at day 10.

# Step 2: Merge preferences back to the original dataframe
plant_region_areas_DAG20 <- plant_region_areas_DAG20 %>%
  left_join(preferences, by = "genotype")
# Purpose: Merges the preference data back into the main dataframe.

# Order by rank (if needed, depending on the requirement)
plant_region_areas_DAG20 <- plant_region_areas_DAG20 %>%
  arrange(desc(rank))
# Purpose: Orders the dataframe by rank for consistent visualization.

## ----r------------------------------------------------------------------------
genotype_DAG20 <- plant_region_areas_DAG20 %>% filter(dag == 20)# dag 20 is the 20th day after germination was picked up
genotype_order <- genotype_DAG20$genotype
time_order <- 0:20


## ----r------------------------------------------------------------------------
# Ensure 'genotype' is a factor with levels in the order specified by 'genotype_order'
plant_region_areas_DAG20_M <- plant_region_areas_DAG20 %>%
  mutate(genotype = factor(genotype, levels = genotype_order))
# Purpose: Converts the genotype column to a factor with specified levels to maintain order in the heatmap.

# Create a heatmap of area differences by day and genotype, separated by preference
ACC_heatmap <- plant_region_areas_DAG20_M %>%
  filter(dag <= 19) %>%
  ggplot(aes(x = dag, y = genotype, fill = area_diff)) +
    geom_tile() +  # Draws the tiles for the heatmap
    scale_fill_gradient2(high = "gold2", low = "seagreen", midpoint = 0) +  # Color gradient
    labs(x = "Days after germination", y = "Accessions") +  # Axis labels
    facet_wrap(~preference, nrow = 2) +  # Separates the heatmap by 'preference'
    coord_cartesian(xlim = c(5, 19)) +  # Limits the x-axis without removing data outside the limits
    theme(axis.text = element_text(size = 6))  # Customizes axis text size

# Display the plot
ACC_heatmap
# Purpose: Creates a heatmap showing the area differences over time for all accessions, grouped by preference.
#ggsave("ACC_rel_growth_time.pdf", plot = ACC_heatmap, width = 20, height = 12, units = "cm", path = "Output/Figures/")

## ----r------------------------------------------------------------------------
# Filter for specific genotypes with the highest opposing feedbacks
#and days less than or equal to 20

# Filter top 20 positive and negative accessions based on area difference
plant_region_areas_DAG20_extreme20 <- plant_region_areas_DAG20 %>%
  filter(genotype_id %in% c("9851", "7350", "5867", "5921", "9972", "10006", 
                            "9823", "7161", "6088", "9084", "6909", "9392", 
                            "8365", "6136", "7217", "506", "7520", "6108", 
                            "9653", "9821", "9591"),
         dag <= 20) %>%
  mutate(group = factor(genotype, 
                        levels = c('IP-Hue-3', 'Tac-0', 'DraIV 3-7', 'Dra2-1', 
                                   'Bozen-1.2', 'Kastel-1', 'IP-Bae-0', 'Gd-1', 
                                   'Stu1-1', 'Lerik2-6', "Col-0", "Hadd-3", 
                                   "Rak-2", "T840", "Lm-2", "BRR60", "Lp2-2", 
                                   "T480", "Giffo-1", "IP-Aru-0", "IP-Vad-0")))

# Pseudo-log10 transformation function
pseudoLog10 <- function(x) { 
  log((x / 2) + sqrt((x / 2)^2 + 1)) / log(10)
}

# Apply transformation
plant_region_areas_DAG20_extreme20 <- plant_region_areas_DAG20_extreme20 %>%
  mutate(pseudo_log_area_diff = pseudoLog10(area_diff))

## ----r------------------------------------------------------------------------
# Load required packages
if (!requireNamespace("dtw", quietly = TRUE)) install.packages("dtw")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) install.packages("ComplexHeatmap")
library(dtw)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

## ----r------------------------------------------------------------------------
# Transform data into wide format
wide_data <- plant_region_areas_DAG20_extreme20 %>%
  select(dag, group, pseudo_log_area_diff) %>%
  pivot_wider(names_from = dag, values_from = pseudo_log_area_diff)

# Ensure columns (days after germination) are sorted numerically
col_order <- order(as.numeric(colnames(wide_data)[-1]))  # Exclude 'group' column
wide_data <- wide_data %>%
  select(group, colnames(wide_data)[-1][col_order])

# Convert to matrix
heatmap_matrix <- as.matrix(wide_data[, -1])  # Exclude 'group'
rownames(heatmap_matrix) <- wide_data$group   # Set row names to genotypes

# **2. Calculate DTW Distance Matrix**
# Using vectorized computation
dtw_dist <- proxy::dist(heatmap_matrix, method = function(x, y) dtw(x, y)$distance)


## ----r------------------------------------------------------------------------
fit <- hclust(dtw_dist, method = "ward.D2")

# Reverse row order for visualization (optional)
row_order <- rev(1:nrow(heatmap_matrix))

## ----r------------------------------------------------------------------------
# **4. Create and Customize the Heatmap**
ht <- Heatmap(
  heatmap_matrix,
  name = "Phenotype Change",
  col = colorRamp2(
    c(min(heatmap_matrix), 0, max(heatmap_matrix)),
    c("orange2", "white", "#0072B2"),
  ),
  cluster_rows = fit,
  cluster_columns = FALSE,
  row_order = row_order,
  show_row_dend = TRUE,
  row_dend_side = "right",
  row_names_side = "left",
  row_title = "",
  row_dend_width = unit(30, "mm"),
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  column_title = "Days after germination",
  column_title_gp = gpar(fontsize = 7),
  column_title_side = "bottom",
  column_names_rot = 0,
  heatmap_legend_param = list(title = "Log Area Diff", legend_direction = "horizontal"),
  show_heatmap_legend = TRUE
)

# Draw the heatmap in the document
draw(ht, heatmap_legend_side = "bottom")

# **5. Export Heatmap to PDF**
output_dir <- "Output/Figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
pdf(file.path(output_dir, "ACC20_rel_growth_time_DTW.pdf"), width = 6, height = 5)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

## ----r------------------------------------------------------------------------
# Filter for observations on day 15, assign preference, and sort
plant_region_areas_DAG15 <- plant_region_areas_DAG20 %>%
  filter(dag == 15) %>%
  mutate(preference = case_when(
    area_diff >= 0 ~ "plus",
    TRUE ~ "minus"
  )) %>%
  arrange(desc(area_diff)) %>%
  mutate(rank = row_number())

## ----r------------------------------------------------------------------------
# Prepare data for specific genotypes to label
labels_data <- plant_region_areas_DAG15 %>%
  filter(genotype %in% c("IP-Hue-3", "Col-0", "Hadd-3")) %>%
  mutate(
    text_position = ifelse(area_diff > 0, area_diff + 4000, area_diff - 5000) # Offset for labels
  )

# Create the plot
ACC_ratio <- ggplot(plant_region_areas_DAG15, aes(x = rank, y = area_diff, fill = preference)) +
  
  # Add bars with original fill color based on "preference"
  geom_bar(stat = "identity", linewidth = 0.001, alpha = 0.8, width = 1) +
  
  # Add a second layer to highlight specific genotypes in black
  geom_bar(
    data = plant_region_areas_DAG15 %>% filter(genotype %in% c("IP-Hue-3", "Col-0", "Hadd-3")),
    aes(x = rank, y = area_diff),
    stat = "identity", fill = "black", linewidth = 0.001, width = 1
  ) +
  
  # Customize the fill scale for the legend
  scale_fill_manual(
    values = c("orange2", "#0072B2"),
    name = "BX feedback",  # Update legend title
    labels = c("negative", "positive")  # Update legend labels
  ) +
  
  # Add points at the exact value of area_diff
  geom_point(data = labels_data,
             aes(y = area_diff), 
             color = "black", size = 1.2) +
  
  # Add text labels dynamically positioned
  geom_text(data = labels_data,
            aes(y = text_position, label = genotype), 
            color = "black", size = 3, hjust = 0.5, vjust = -0.2) +  # Increased text size
  
  xlim(-2, 412) +
  
  # Modify the y-axis label and scale
  labs(
    y = expression("Growth feedback (10"^4*")"), 
    x = "Accession rank"
  ) +
  
  # Scale y-axis to divide values by 10,000 for display
  scale_y_continuous(labels = function(x) x / 10000) +
  
  ggtitle("Growth feedback distribution of accessions") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8),  # Slightly increased font size
    axis.title = element_text(size = 9), 
    plot.title = element_text(size = 9, face = "bold"),
    legend.position = c(0.95, 0.95),  # Custom position for the top-right corner
    legend.justification = c(1, 1),   # Align legend box to top-right
    legend.background = element_rect(fill = "white", color = "black"),  # Optional: Add background
    legend.title = element_text(size = 7),  # Increased legend title size
    legend.text = element_text(size = 6)    # Increased legend items size
  )

# Display the plot
ACC_ratio

# Save the plot
ggsave("ACC_distr_rel_feedback.pdf", plot = ACC_ratio, width = 16, height = 7, units = "cm", path = "Output/Figures/")

## ----r------------------------------------------------------------------------
# Filter data for genotype_id 9851 (IP-Hue-3) and DAG ≤ 20
plant_region_areas_IPHue3 <- plant_region_areas %>%
  filter(dag <= 20, genotype_id == 9851)

# Calculate mean plant region area for each bx and dag combination
plant_region_areas_IPHue3 <- plant_region_areas_IPHue3 %>%
  group_by(bx, dag) %>%
  summarise(mean = mean(plant_region_area, na.rm = TRUE), .groups = "drop")

# Remove unnecessary columns and ensure unique rows
plant_region_areas_IPHue3 <- distinct(plant_region_areas_IPHue3)

# Transform data: separate columns for each treatment
shaded_data_IP_Hue <- plant_region_areas_IPHue3 %>%
  pivot_wider(names_from = bx, values_from = mean, names_prefix = "mean_")

# Ensure correct column names (assuming `+` and `-` are original bx values)
shaded_data_IP_Hue <- shaded_data_IP_Hue %>%
  rename_with(~ gsub("mean_\\+", "mean_group1", .x)) %>%
  rename_with(~ gsub("mean_\\-", "mean_group2", .x)) %>%
  mutate(
    lower = pmin(mean_group1, mean_group2, na.rm = TRUE),
    upper = pmax(mean_group1, mean_group2, na.rm = TRUE)
  )

# Define x-axis breaks and labels
x_breaks <- seq(floor(min(plant_region_areas_IPHue3$dag)), 
                ceiling(max(plant_region_areas_IPHue3$dag)), by = 5)
x_labels <- as.character(x_breaks)

# Display final data structure
#head(shaded_data_IP_Hue)

## ----r------------------------------------------------------------------------
# Create the plot
IP_Hue_3 <- ggplot() +
  # Add ribbon without legend
  geom_ribbon(
    data = shaded_data_IP_Hue,
    aes(x = dag, ymin = lower, ymax = upper),
    fill = "#0072B2", alpha = 0.2,
    show.legend = FALSE
  ) +
  # Add points with legend
  geom_point(
    data = plant_region_areas_IPHue3,
    aes(x = dag, y = mean, colour = bx, fill = bx),
    shape = 21, size = 2.5, alpha = 0.4
  ) +
  # Add line without legend
  geom_line(
    data = plant_region_areas_IPHue3,
    aes(x = dag, y = mean, colour = bx),
    show.legend = FALSE
  ) +
  # Modify the y-axis label and scale
  ggtitle("Daily rosette area: IP-Hue-3") +
  labs(
    y = expression("Rosette area (10"^4*")"), 
    x = "Days after germination"
  ) +
  scale_y_continuous(labels = function(x) x / 10000) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  # Define manual color and fill scales with legend settings
  scale_colour_manual(
    values = c("seagreen", "gold2"),
    guide = "none" # Suppress colour legend
  ) +
  scale_fill_manual(
    values = c("seagreen", "gold2"),
    name = "Microbiome",  # Update legend title
    labels = c("BX-", "BX+"),  # Update legend labels
    guide = guide_legend(
      override.aes = list(shape = 21, size = 2.5, alpha = 0.4) # Match point aesthetics
    )
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    plot.title = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.title.x = element_text(size = 7),
    legend.position = c(0.2, 0.95),  # Position legend in top-right
    legend.justification = c(1, 1),   # Align legend box to top-right
    legend.background = element_rect(fill = "white", color = "black")
  )

# Display the plot
print(IP_Hue_3)

# Purpose: Displays the plot.

# Save the plot
ggsave("IP_hue_3_growth.pdf", plot = IP_Hue_3, width = 10, height = 10, units = "cm", path = "Output/Figures/")
# Purpose: Saves the plot as a pdf file.

## ----r------------------------------------------------------------------------
# Compute the difference between the two groups
shaded_data_IP_Hue <- shaded_data_IP_Hue %>%
  mutate(difference = mean_group1 - mean_group2)
# Purpose: Calculates the difference in mean plant region area between the two treatments for each day.

# Plot with shaded area and differences as bars
IP_Hue_3_response <- ggplot() +
  geom_bar(data = shaded_data_IP_Hue, aes(x = dag, y = difference), stat = "identity", fill = "#0072B2", alpha = 0.8) +
  ggtitle("Daily growth feeedbacks: IP-Hue-3") +
  labs(
    y = expression("Growth feedback (10"^4*")"), 
    x = "Days after germination"
  ) +
  scale_y_continuous(labels = function(x) x / 10000) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  theme_bw() +
  theme(axis.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        plot.title = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        legend.position = "none")       # Remove legend
# Purpose: Creates a bar plot showing the daily growth feedback for IP-Hue-3, calculated as the difference between the two treatments.

# Display the plot
print(IP_Hue_3_response)
# Purpose: Displays the plot.

# Save the plot
ggsave("IP_Hue_3_response.pdf", plot = IP_Hue_3_response, width = 10, height = 10, units = "cm", path = "Output/Figures/")
# Purpose: Saves the plot as a PNG file.

## ----r------------------------------------------------------------------------
# Combine plots for IP-Hue-3 using cowplot

# Load the cowplot library
library(cowplot)

# Combine the growth plot and the growth feedback plot
IP_Hue_3_combined <- plot_grid(IP_Hue_3, IP_Hue_3_response, nrow = 1, align = "h")

# Display the combined plot

print(IP_Hue_3_combined)

# Save the combined plot as a PNG file

ggsave("IP_Hue_3_combined.pdf", plot = IP_Hue_3_combined, width = 18, height = 8, units = "cm", path = "Output/Figures/")


## ----r------------------------------------------------------------------------
# Filter data for genotype_id 9392
plant_region_areas_Hadd3 <- plant_region_areas %>% 
  filter(dag < 21, genotype_id == 9392)
# Purpose: Filters the dataset to include only data for the accession with genotype_id 9392 (Hadd-3) and days less than 20.

# Calculate mean plant region area for each combination of bx and dag
plant_region_areas_Hadd3 <- plant_region_areas_Hadd3 %>%
  group_by(bx, dag) %>%
  mutate(mean = mean(plant_region_area, na.rm = TRUE)) %>%
  ungroup()
# Purpose: Calculates the mean plant region area for each treatment (bx) and day (dag) combination.

# Remove unnecessary columns and keep only unique entries
plant_region_areas_Hadd3 <- plant_region_areas_Hadd3 %>%
  select(-c(1, 6)) %>%
  distinct()
# Purpose: Removes unnecessary columns and ensures that only unique rows remain in the dataset.

# Create a new data frame for the shaded area between curves
# Spread the data to have separate columns for each group
shaded_data_Hadd3 <- plant_region_areas_Hadd3 %>%
  pivot_wider(names_from = bx, values_from = mean)
# Purpose: Transforms the data to have separate columns for each treatment, preparing it for shaded area plotting.

# Ensure the correct names for the new columns
# Assuming the groups are "group1" and "group2"
shaded_data_Hadd3 <- shaded_data_Hadd3 %>%
  rename(mean_group1 = `+`, mean_group2 = `-`) %>%
  mutate(lower = pmin(mean_group1, mean_group2, na.rm = TRUE),
         upper = pmax(mean_group1, mean_group2, na.rm = TRUE))
# Purpose: Renames the columns for clarity and calculates the lower and upper bounds for the shaded area.

## ----r------------------------------------------------------------------------


# Create the plot
Hadd_3 <- ggplot() +
  # Add ribbon without legend
  geom_ribbon(
    data = shaded_data_Hadd3,
    aes(x = dag, ymin = lower, ymax = upper),
    fill = "orange2", alpha = 0.2,
    show.legend = FALSE
  ) +
  # Add points with legend
  geom_point(
    data = plant_region_areas_Hadd3,
    aes(x = dag, y = mean, colour = bx, fill = bx),
    shape = 21, size = 2.5, alpha = 0.4
  ) +
  # Add line without legend
  geom_line(
    data = plant_region_areas_Hadd3,
    aes(x = dag, y = mean, colour = bx),
    show.legend = FALSE
  ) +
  # Title and axis labels
  ggtitle("Daily rosette area: Hadd-3") +
  labs(
    y = expression("Rosette area (10"^4*")"), 
    x = "Days after germination"
  ) +
  # Scales for y-axis and x-axis
  scale_y_continuous(labels = function(x) x / 10000) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  # Define manual color and fill scales with legend settings
  scale_colour_manual(
    values = c("seagreen", "gold2"),
    guide = "none" # Suppress colour legend
  ) +
  scale_fill_manual(
    values = c("seagreen", "gold2"),
    name = "Microbiome",  # Update legend title
    labels = c("BX-", "BX+"),  # Update legend labels
    guide = guide_legend(
      override.aes = list(shape = 21, size = 2.5, alpha = 0.4) # Match point aesthetics
    )
  ) +
  # Theme adjustments
  theme_bw() +
  theme(
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    plot.title = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.title.x = element_text(size = 7),
    legend.position = c(0.2, 0.95),  # Position legend in top-right
    legend.justification = c(1, 1),   # Align legend box to top-right
    legend.background = element_rect(fill = "white", color = "black")
  )

# Display the plot
Hadd_3

ggsave("Hadd_3_growth.pdf", plot = Hadd_3, width = 10, height = 10, units = "cm", path = "Output/Figures/")

## ----r------------------------------------------------------------------------
# Compute the difference between the two groups
shaded_data_Hadd3 <- shaded_data_Hadd3 %>%
  mutate(difference = mean_group1 - mean_group2)
# Purpose: Calculates the difference in mean plant region area between the two treatments for each day.

# Plot with shaded area and differences as bars
Hadd3_response <- ggplot() +
  geom_bar(data = shaded_data_Hadd3, aes(x = dag, y = difference), stat = "identity", fill = "orange2", alpha = 0.8) +
  ggtitle("Daily growth feeedbacks: Hadd-3") +
  labs(
    y = expression("Growth feedback (10"^4*")"), 
    x = "Days after germination"
  ) +
  scale_y_continuous(labels = function(x) x / 10000) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  theme_bw() +
  theme(axis.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        plot.title = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        #axis.text.y = element_blank(),  # Remove y-axis text
        legend.position = "none")       # Remove legend
# Purpose: Creates a bar plot showing the daily growth feedback for Hadd-3, calculated as the difference between the two treatments.

# Display the plot
print(Hadd3_response)
# Purpose: Displays the plot.

# Save the plot
ggsave("Hadd_3_response.pdf", plot = Hadd3_response, width = 10, height = 10, units = "cm", path = "Output/Figures/")

## ----r------------------------------------------------------------------------
# Combine plots for IP-Hue-3 using cowplot


# Combine the growth plot and the growth feedback plot
Hadd_3_combined <- plot_grid(Hadd_3, Hadd3_response, nrow = 1, align = "h")

# Display the combined plot

print(Hadd_3_combined)

# Save the combined plot as a PNG file

ggsave("Hadd_3_combined.pdf", plot = Hadd_3_combined, width = 18, height = 8, units = "cm", path = "Output/Figures/")

## ----r------------------------------------------------------------------------
# Combine the growth plot and the growth feedback plot
Extreme_phenotypes <- plot_grid(IP_Hue_3_combined, Hadd_3_combined, nrow = 2, align = "h")

# Display the combined plot

print(Extreme_phenotypes)
# Save the plot
ggsave("Extreme_ACC_phenotypes.pdf", plot = Extreme_phenotypes, width = 14, height = 13, units = "cm", path = "Output/Figures/")
# Purpose: Saves the combined growth feedback plot as a PDF file.

## ----r------------------------------------------------------------------------
ACC_coordinates <- read_excel("Input/Dataset1.xlsx", sheet = "accession_coordinates")
kable(ACC_coordinates[1:3, ], caption = "Coordinate table.")

## ----r------------------------------------------------------------------------
# Filter for observations on day 15, assign preference, and sort
plant_region_areas_DAG10 <- plant_region_areas_DAG20 %>%
  filter(dag == 15) %>%
  mutate(preference = case_when(
    area_diff >= 0 ~ "plus",
    TRUE ~ "minus"
  )) %>%
  arrange(desc(area_diff)) %>%
  mutate(rank = row_number())

## ----r------------------------------------------------------------------------
# Merge the data tables and include only Lat and Long columns from ACC_coordinates
coordinates_merged_data <- plant_region_areas_DAG10 %>%
  left_join(
    ACC_coordinates %>% select(ACC_ID, Lat, Long),  # Select only relevant columns
    by = c("genotype_id" = "ACC_ID")               # Specify the join key
  )

## ----r------------------------------------------------------------------------
# Perform Pearson correlation test for Latitude
cor_test_lat <- cor.test(coordinates_merged_data$Lat, coordinates_merged_data$area_diff, 
                         method = "pearson", use = "complete.obs")

# Perform Pearson correlation test for Longitude
cor_test_long <- cor.test(coordinates_merged_data$Long, coordinates_merged_data$area_diff, 
                          method = "pearson", use = "complete.obs")

# Function to extract correlation test results safely
extract_cor_results <- function(test_result, variable_name) {
  data.frame(
    Variable = variable_name,
    Correlation = round(test_result$estimate, 3),  # Rounded for readability
    p_value = round(test_result$p.value, 4),
    CI_lower = if (!is.null(test_result$conf.int)) round(test_result$conf.int[1], 3) else NA,
    CI_upper = if (!is.null(test_result$conf.int)) round(test_result$conf.int[2], 3) else NA
  )
}

# Extract results
lat_results <- extract_cor_results(cor_test_lat, "Latitude")
long_results <- extract_cor_results(cor_test_long, "Longitude")

# Combine both results into a single table
cor_results <- bind_rows(lat_results, long_results)

# Display the results as a neat table
kable(cor_results, caption = "Correlation Test Results for Latitude and Longitude", digits = 4)

## ----r------------------------------------------------------------------------
# Linear regression with Latitude
model_lat <- lm(area_diff ~ Lat + genotype_id, data = coordinates_merged_data)
summary_lat <- summary(model_lat)

# Extract coefficients and R-squared
intercept_lat <- coef(summary_lat)[1, 1]
slope_lat <- coef(summary_lat)[2, 1]
r_squared_lat <- summary_lat$r.squared
p_value_lat <- coef(summary_lat)[2, 4]
p_value_f_lat <- summary_lat$fstatistic[1] # F-statistic value
print(summary_lat)
if (p_value_lat < 0.05) {
    print("Latitude significantly affects area_diff.")
} else {
    print("Latitude does not significantly affect area_diff.")
}
# Purpose: Performs a linear regression to examine the effect of latitude on area difference, extracting relevant statistics and checking significance.

# Linear regression with Longitude
model_long <- lm(area_diff ~ Long + genotype_id, data = coordinates_merged_data)
summary_long <- summary(model_long)

# Extract coefficients and R-squared
intercept_long <- coef(summary_long)[1, 1]
slope_long <- coef(summary_long)[2, 1]
r_squared_long <- summary_long$r.squared
p_value_long <- coef(summary_long)[2, 4]
p_value_f_long <- summary_long$fstatistic[1] # F-statistic value
print(summary_long)
if (p_value_lat < 0.05) {
    print("Longitute significantly affects area_diff.")
} else {
    print("Longitue does not significantly affect area_diff.")
}

## ----r------------------------------------------------------------------------
# Linear regression with Latitude
model_lat <- lm(area_diff ~ Lat, data = coordinates_merged_data)
summary_lat <- summary(model_lat)

# Extract coefficients and R-squared for Latitude
lat_results <- data.frame(
  Variable = "Latitude",
  Intercept = coef(summary_lat)[1, 1],
  Slope = coef(summary_lat)[2, 1],
  R_squared = summary_lat$r.squared,
  p_value_slope = coef(summary_lat)[2, 4],
  F_statistic = summary_lat$fstatistic[1]
)

# Linear regression with Longitude
model_long <- lm(area_diff ~ Long, data = coordinates_merged_data)
summary_long <- summary(model_long)

# Extract coefficients and R-squared for Longitude
long_results <- data.frame(
  Variable = "Longitude",
  Intercept = coef(summary_long)[1, 1],
  Slope = coef(summary_long)[2, 1],
  R_squared = summary_long$r.squared,
  p_value_slope = coef(summary_long)[2, 4],
  F_statistic = summary_long$fstatistic[1]
)

# Combine the results for Latitude and Longitude
regression_results <- rbind(lat_results, long_results)

# Display the results as a neat table
kable(regression_results, caption = "Linear Regression Results for Latitude and Longitude", digits = 4)

## ----r------------------------------------------------------------------------
# Compute linear regression model
lm_lat <- lm(area_diff ~ Lat, data = coordinates_merged_data, na.action = na.omit)

# Extract R² and p-value
r_squared_lat <- summary(lm_lat)$r.squared
p_value_lat <- summary(lm_lat)$fstatistic[1] %>% 
  pf(df1 = 1, df2 = lm_lat$df.residual, lower.tail = FALSE)

# Define y-axis limits dynamically (prevents cut-off issues)
y_limits <- c(min(coordinates_merged_data$area_diff, na.rm = TRUE) * 1.1, 
              max(coordinates_merged_data$area_diff, na.rm = TRUE) * 1.1)

# Create the scatter plot
corr_plot_lat <- ggplot(coordinates_merged_data, 
                        aes(x = Lat, y = area_diff, color = preference)) +
  
  # Scatter points
  geom_point(size = 1, alpha = 0.8) +
  
  # Custom color scale
  scale_colour_manual(
    values = c("orange2", "#0072B2"),
    name = "BX feedback",
    labels = c("negative", "positive")
  ) +
  
  # Regression line
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", 
              color = "red4", linewidth = 1) +
  
  # Labels
  labs(
    title = "Phenotype vs. Latitude with Regression Line",
    x = "Latitude",
    y = "Phenotype (Area Diff)"
  ) +
  
  # Theme modifications for A4 compatibility
  theme_classic() +
  ylim(y_limits) +
  theme(
    legend.position = c(0.97, 0.03),  # Bottom-right inside plot
    legend.justification = c(1, 0),  
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.title = element_text(size = 8, face = "bold"),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6)
  ) +
  
  # Add R² and p-value annotation
  annotate(
    "text", x = max(coordinates_merged_data$Lat, na.rm = TRUE), 
    y = min(y_limits) + (max(y_limits) - min(y_limits)) * 0.85,  # Position dynamically
    label = sprintf("R² = %.2f\np (F) = %.3g", r_squared_lat, p_value_lat),
    hjust = 1, vjust = 0, size = 3, color = "black"
  )

# Display the plot
print(corr_plot_lat)


## ----r------------------------------------------------------------------------
# Compute linear regression model
lm_long <- lm(area_diff ~ Long, data = coordinates_merged_data, 
              na.action = na.omit)

# Extract R² and p-value
r_squared_long <- summary(lm_long)$r.squared
p_value_long <- summary(lm_long)$fstatistic[1] %>% 
  pf(df1 = 1, df2 = lm_long$df.residual, lower.tail = FALSE)

# Define y-axis limits dynamically (prevents cut-off issues)
y_limits <- c(min(coordinates_merged_data$area_diff, na.rm = TRUE) * 1.1, 
              max(coordinates_merged_data$area_diff, na.rm = TRUE) * 1.1)

# Create the scatter plot
corr_plot_long <- ggplot(coordinates_merged_data, aes(x = Long, y = area_diff, 
                                                      color = preference)) +
  
  # Scatter points
  geom_point(size = 1, alpha = 0.8) +
  
  # Custom color scale
  scale_colour_manual(
    values = c("orange2", "#0072B2"),
    name = "BX feedback",
    labels = c("negative", "positive")
  ) +
  
  # Regression line
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", 
              color = "red4", linewidth = 1) +
  
  # Labels
  labs(
    title = "Phenotype vs. Longitude with Regression Line",
    x = "Longitude",
    y = "Phenotype (Area Diff)"
  ) +
  
  # Theme modifications for A4 compatibility
  theme_classic() +
  ylim(y_limits) +
  theme(
    legend.position = c(0.97, 0.03),  # Bottom-right inside plot (same as lat plot)
    legend.justification = c(1, 0),  
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.title = element_text(size = 8, face = "bold"),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6)
  ) +
  
  # Add R² and p-value annotation
  annotate(
    "text", x = max(coordinates_merged_data$Long, na.rm = TRUE), 
    y = min(y_limits) + (max(y_limits) - min(y_limits)) * 0.85,  # Position dynamically
    label = sprintf("R² = %.2f\np (F) = %.3g", r_squared_long, p_value_long),
    hjust = 1, vjust = 0, size = 3, color = "black"
  )

# Display the plot
print(corr_plot_long)

## ----r------------------------------------------------------------------------
# Add simple data points to the map
# Combine the two plots on top of each other
ACC_corr <- plot_grid(
  corr_plot_lat, 
  corr_plot_long, 
  labels = c("B", ""), 
  label_size = 12, 
  nrow = 1, 
  ncol = 2, 
  rel_heights = c(1, 1)
  # Adjust relative heights if necessary
)
# Purpose: Combines the scatter plots for latitude and longitude into a single figure, displayed side by side.

# Display the combined plot
print(ACC_corr)
# Purpose: Displays the combined scatter plot.

# Save the plot
ggsave("ACC_corr.pdf", plot = ACC_corr, width = 20, height = 10, units = "cm", path = "Output/Figures/")
# Purpose: Saves the combined scatter plot as a PDF file.

## ----r------------------------------------------------------------------------
# Load necessary libraries
library(maps)
# Purpose: Loads the `maps` package for map plotting.

# Obtain world map data
world_map <- map_data("world")
# Purpose: Retrieves the data necessary to create a world map.

# Subset the world map data for the region of interest
World_map_small <- subset(world_map, long > -130 & long < 90)
World_map_small <- subset(World_map_small, lat > -0 & lat < 110)
# Purpose: Limits the map data to the specific region of interest by filtering longitude and latitude ranges.

# Create a base plot with ggplot2
p <- ggplot() + coord_fixed() +
  xlab("Longitude") + ylab("Latitude")
# Purpose: Sets up the basic structure for the map plot with fixed coordinates.

# Add map to base plot
base_world_messy <- p + geom_polygon(data = World_map_small, aes(x = long, 
                                                                 y = lat, 
                                                                 group = group), 
                                     colour = "grey", fill = "grey", alpha = 0.9)
# Purpose: Adds the filtered map data to the plot as polygons representing countries.

# Display the messy base world map
#print(base_world_messy)
# Purpose: Displays the initial, unrefined world map.

# Strip the map down to make it clean and beautiful
cleanup <- theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(), 
                 panel.background = element_rect(fill = 'white', 
                                                 colour = 'white'), 
                 axis.line = element_line(colour = "grey"))
# Purpose: Cleans up the map by removing unnecessary elements and setting a minimalist theme.

base_world <- base_world_messy + cleanup
# Purpose: Applies the cleaned-up theme to the world map.

# Display the clean base world map
#print(base_world)
# Purpose: Displays the refined world map.

# Separate data for layers
df_layer_1 <- coordinates_merged_data[coordinates_merged_data$preference == "plus", ]
df_layer_2 <- coordinates_merged_data[coordinates_merged_data$preference == "minus", ]
# Purpose: Separates the coordinates data into two layers based on whether the response is positive or negative.

## ----r------------------------------------------------------------------------
# Add simple data points to the map
# Add a new column 'growth_type' to both data frames if it doesn't exist
# df_layer_1$growth_type <- "plus"
# df_layer_2$growth_type <- "minus"

ACC_map <- base_world +
  # All points plotted with color based on growth_type
  geom_point(
    data = df_layer_2,
    aes(x = Long, y = Lat, color = "minus"),
    pch = 19,
    size = 1.4,
    alpha = 0.5
  ) +
  geom_point(
    data = df_layer_1,
    aes(x = Long, y = Lat, color = "plus"),
    pch = 19,
    size = 1.4,
    alpha = 0.5
  ) +
  scale_color_manual(
    values = c("plus" = "#0072B2", "minus" = "orange2"),
    name = "Growth feedback",
    labels = c("plus" = "positive", "minus" = "negative")
  ) +
  theme(
    legend.position = c(0.01, 0.01),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7, face = "bold")
  )

ACC_map
# Purpose: Displays the world map with the positive and negative accession data points.

# Save the plot as a PDF
ggsave("ACC_response_map_white.pdf", plot = ACC_map, height = 15, width = 20, units = "cm", path = "Output/Figures/")
# Purpose: Saves the map with data points as a PDF file.

## ----r------------------------------------------------------------------------
# Get session information
SI <- sessionInfo()

# Extract the attached packages and their versions
attached_packages <- SI$otherPkgs

# Print the R version
cat("R version:", SI$R.version$version.string, "\n\n")

# Print the attached packages and their versions
cat("Attached packages:\n")
invisible(lapply(attached_packages, function(pkg) {
  cat(pkg$Package, "-", pkg$Version, "\n")
}))


## ----r create package library, echo = FALSE-----------------------------------
knitr::write_bib(file = 'packages.bib')

