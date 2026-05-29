#' Subset Plotting data for gggenomes
#'
#' This function takes a ploting data list for gggenomes created via `gs_to_gggenomes` and performs adjustments,
#' effectively cutting off parts of the sequence and subsetting the resulting dataframe.
#' It will also create polygon links and store those in $polygons
#' Currently this is intended to cut of parts of the sequence on the left (5')
#'
#' @param plotdat a list of dataframes, containing dataframes named `seqs`, `genes` and `links`
#' @param seqs_to_adjust name of the sequences to adjust, corresponding to `$seq_id` in `seqs`
#' @param shift length to shift, default is 300000
#' @param seqs_to_keep list of sequnences (seq_ids) to keep, if `NULL` all genomes will be kept
#' 
#'
#' @return a list of dataframes, identical to the list provided as `plotdat` with adjustments to `seqs`, `genes` and `links`,
#' and an extra dataframe named `polygons` that can be used to draw curved polygons as links via `geom_polygon`.
#'  
#' @export
#'
#' @examples
#' 
#' @author Niklas Schandry
subset_plotdat <- function(plotdat = NULL , seqs_to_adjust = NULL, shift = 300000, seqs_to_keep = NULL, debug = F) {
    if(is.null(plotdat)) {
      stop("Please provide plotting data")
    }
    if(!is.list(plotdat)) {
      stop("plotdat needs to be a list of dataframes")
    }
    if(!is.character(seqs_to_adjust)) {
      stop("seqs_to_adjust needs to be a character vector")
    }
    if(sum(!(seqs_to_adjust %in% plotdat$seqs$seq_id)) > 0 ) {
      stop("sequences to adjust need to all be part of plotdat$seqs$seq_id")
    }
    if(is.null(seqs_to_keep)) {
      seqs_to_keep <- unique(plotdat$seqs$seq_id)
    }
  pdat <- plotdat
  message("Adjusting sequences")
  pdat$seqs <- pdat$seqs %>%
    filter(seq_id %in% seqs_to_keep) %>% 
    mutate(length = case_when(seq_id %in% seqs_to_adjust ~ length - shift,
                              TRUE ~ length))
  message("Adjusting genes")
  pdat$genes <- pdat$genes %>%
    filter(seq_id %in% seqs_to_keep) %>% 
    mutate(
      start = case_when(seq_id %in% seqs_to_adjust ~ start - shift,
                        TRUE ~ start),
      end = case_when(seq_id %in% seqs_to_adjust ~ end - shift,
                      TRUE ~ end)
    ) %>%
    filter(start > 0 & end > 0)
  message("Adjusting links")
  
  pdat$links <- pdat$links %>%
    filter(seq_id %in% seqs_to_keep & seq_id2 %in% seqs_to_keep) %>% 
    mutate(
      start = case_when(seq_id %in% seqs_to_adjust ~ start - shift,
                        TRUE ~ start),
      end = case_when(seq_id %in% seqs_to_adjust ~ end - shift,
                      TRUE ~ end),
      length = case_when(seq_id %in% seqs_to_adjust ~ length - shift,
                         TRUE ~ length),
      start2 = case_when(seq_id2 %in% seqs_to_adjust ~ start2 - shift,
                         TRUE ~ start2),
      end2 = case_when(seq_id2 %in% seqs_to_adjust ~ end2 - shift,
                       TRUE ~ end2),
      length2 = case_when(seq_id2 %in% seqs_to_adjust ~ length2 - shift,
                          TRUE ~ length2)
    ) %>%
    filter(start > 0 & end > 0 & start2 > 0 & end2 > 0)
  if(debug) {print(pdat)}
  message("Computing polygons")
  pdat$polygons <- compute_polygons(plotdat = pdat,
                                    as.character(pdat$seqs$seq_id))
  return(pdat)
}

# Make a list length containing numbers, names according to genome order
#' Get genomes
#'
#' @param plotdat list of dataframes, as returned by `gs_to_gggenomes()`
#' @param genome_order desired order of genomes
#'
#' @return Inverted order of sequences to arrange in gggenomes.
#' @export
#'
#' @examples
#' 
#' @author Niklas Schandry
get_genomes <- function(plotdat, genome_order) {
  genomes_numbers <- 
    c(1:nrow(plotdat$seqs)) %>% 
    as.list
  names(genomes_numbers) <- rev(genome_order)
  genomes_numbers
}

# Compute polygons via GENESPACE::calc_curvePolygon
#' Compute polygons
#'
#' Compute polygons to draw links in gggenomes via geom_poly. 
#' This function wraps `GENESPACE::calc_curvePolygon`
#' @param plotdat list of dataframes, as returned by `gs_to_gggenomes()`
#' @param genome_order desired order of genomes in the plot (top to bottom).
#'
#' @return dataframe containing polygons. 
#' Column `direct` indicates if this is direct link between two sequences drawn next to each other.
#' I recommend to subset for `direct` for larger figures to avoid drawing many links across the plot.
#' 
#' @export
#'
#' @examples
#' 
#' @author Niklas Schandry
compute_polygons <- function(plotdat, genome_order) {
  # Get the order of genomes
  genomes <- get_genomes(plotdat, genome_order)
  # Compute all polygons
  parallel::mclapply(1:nrow(plotdat$links),
                     \(row) {
                       plotdat$links[row,] %$%
                         GENESPACE::calc_curvePolygon(
                           start1 = start,
                           end1 = end,
                           start2 = start2,
                           end2 = end2,
                           y1 = genomes[[paste(seq_id)]],
                           y2 = genomes[[paste(seq_id2)]]
                         ) %>%
                         as.data.frame() %>%
                         mutate(
                           link = as.integer(row),
                           og = plotdat$links[row,]$og,
                           group = plotdat$links[row,]$group,
                           group = case_when(is.na(group) ~ "Other",
                                           TRUE ~ group)
                         )
                     },
                     mc.cores = 8
                    ) %>%
    bind_rows() %>%
    group_by(link) %>%
    mutate(direct = case_when(max(y) - min(y)  == 1 ~ TRUE,
                              TRUE ~ FALSE)) %>%
    ungroup()
}

