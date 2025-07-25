seq_nr_plot <- function(df, title = "") {
  # Ensure 'line' is a factor matching the colors
  df$line <- factor(df$line, levels = names(level_cols_line))
  
  ggplot(df, aes(x = condition, y = col_sum, fill = line)) + 
    geom_boxplot(
      position = position_dodge2(width = 0.75, preserve = "single"),
      outlier.colour = NA
    ) +
    geom_jitter(
      aes(color = line, group = line),
      size = 1,
      stroke = 0.2,
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)
    ) +
    facet_grid(. ~ condition, scales = "free_x", space = "free") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = 5)
    ) +
    ylab("log2(seq)") +
    scale_fill_manual(values = level_cols_line) +
    scale_color_manual(values = level_cols_line) +
    scale_y_continuous(trans = "log2") +
    ggtitle(title)
}


