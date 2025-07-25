samples_plot <- function(df) {
  # Ensure 'line' is a factor matching the colors
  df$line <- factor(df$line, levels = names(level_cols_line))
  
  ggplot(df, aes(x = x, y = col_sum, fill = line, color = line)) +
    geom_bar(
      stat = "identity",
      width = 1,
      color = "black",
      size = 0.1
    ) +
    scale_fill_manual(values = level_cols_line) +
    scale_colour_manual(values = level_cols_line) +
    xlab("samples") +
    ylab("seq / samples") +
    facet_grid(. ~ condition, scales = "free_x", space = "free") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      legend.position = "none"
    )
}


