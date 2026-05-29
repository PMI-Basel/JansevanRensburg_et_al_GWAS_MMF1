#' A dark theme to plot genespace results via gggenomes
#'
#' @param base_size Base text size
#' @param base_family Base text family
#' @param base_line_size Base line size
#' @param base_rect_size Base rectangle size
#'
#' @returns
#' @export
#'
#' @examples
#' @author Niklas Schandry
theme_ggg_dark <- function(base_size = 14,
                           base_family = "",
                           base_line_size = 14 / 30,
                           base_rect_size = 14 / 30) {
  theme(
    line = element_line(color = "white"),
    text = element_text(color = "white"),
    rect = element_rect(fill = "black", colour = NA),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "black", colour = NA),
    # x-axis
    axis.line.x = element_line(color = "white", size = .4),
    axis.title.x = element_blank(),
    axis.text.x = element_text(color = "white", size = 7),
    axis.ticks.length.x = unit(.7, "mm"),
    # y-axis
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    legend.background = element_rect(fill = "black", colour = NA),
    legend.box.background = element_rect(fill = "black", colour = NA),
    legend.key = element_rect(fill = "black", colour = NA),
    strip.background = element_rect(fill = "black", colour = NA),
    strip.text = element_text(color = "white"),
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    axis.ticks = element_line(color = "white"),
    axis.line = element_line(color = "white")
  )
}

#' A light theme to plot genespace results via gggenomes
#'
#' @param base_size Base text size
#' @param base_family Base text family
#' @param base_line_size Base line size
#' @param base_rect_size Base rectangle size
#'
#' @returns
#' @export
#'
#' @examples
#' @author Niklas Schandry
theme_ggg_light <- function(base_size = 14,
                            base_family = "",
                            base_line_size = base_size/ 30,
                            base_rect_size = base_size / 30) {
  theme(
    line = element_line(color = "black"),
    text = element_text(color = "black",size = base_size),
    rect = element_rect(fill = "white", colour = NA),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    # x-axis
    axis.line.x = element_line(color = "black", size = base_line_size),
    axis.title.x = element_blank(),
    axis.text.x = element_text(color = "black", size = base_size),
    axis.ticks.length.x = unit(.7, "mm"),
    # y-axis
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.box.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black")
  )
}