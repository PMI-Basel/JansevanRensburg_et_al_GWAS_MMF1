#' gs_to_ggenomes
#'
#' Convert genespace output to a format that can be plotted by gggenomes
#'
#' @param genespace_output genespace output
#' @param bed bed file describing the region of interest in the reference genome
#' @param genomes genomes used for genespace
#' @param pad add padding to non-reference genomes
#' @param annotations a dataframe of gene annotations.
#' @param ncores number of cores to use (default 24)
#' @param pangenes_only Return pangenes, if this is TRUE this function is a wrapper for GENESPACE::query_pangenes
#'
#' @return
#' Returns a list of dataframes
#'
#'    $genes: a table of genes for each genome
#'    $seqs: list of sequences and their lengths
#'    $links = table containing all vs all synteny relationshops
#'    $query = as_tibble(query_features), the features of all non-reference genomes
#'    $ref = as_tibble(ref_features), features of reference genome
#'    $dat = combination of query and ref
#'
#' @export
#'
#' @examples
#'
#' @author Niklas Schandry

gs_to_gggenomes <- function(genespace_output,
                            bed,
                            genomes,
                            annotations,
                            pad = 0,
                            ncores = 24,
                            pangenes_only = F) {
  start_fn_time <- now()
  start_time <- now()
  # Get pangenens
  pangenes <- query_pangenes(genespace_output,
                             bed,
                             transform = FALSE) %>%
    .[[1]] %>%
    dplyr::select(genome, feat_id = id, og, pgID) %>%
    dplyr::mutate(feat_id = glue::glue("{genome}_{feat_id}")) %>%
    mutate(og = as.factor(og))
  
  message(glue::glue("# Obtained pangenes in {as.numeric(difftime(now(),start_time, units = 'secs'))  %>% round(2)} seconds"))
  if (pangenes_only) {
    return(pangenes)
  }
  # This part takes long,
  # I forgot why it is necessary to do it that way
  # There was a reason though
  message(glue::glue("# Obtaining syntenic relationships"))
  start_time <- now()
  dat <-  prepare_genespace_synteny(
    genespace_output,
    bed_table = bed,
    genomes = genomes,
    padding = pad
  )
  message(glue::glue("# Obtained syntenic relationships in {as.numeric(difftime(now(),start_time, units = 'secs'))  %>% round(2)} seconds"))
  start_time <- now()
  
  ref_features <- dtplyr::lazy_dt(dat) %>%
    dplyr::select(
      seq_id = chr1,
      start = start1,
      end = end1,
      feat_id = id1,
      genome = genome1
    ) %>%
    dplyr::mutate(
      seq_id = glue::glue_data(., "{genome}_{seq_id}") %>% as.character,
      feat_id = glue::glue_data(., "{genome}_{feat_id}") %>% as.character
    ) %>%
    dplyr::distinct() %>%
    as.data.frame()
  
  query_features <- dtplyr::lazy_dt(dat) %>%
    dplyr::select(
      seq_id = chr2,
      start = start2,
      end = end2,
      feat_id = id2,
      genome = genome2
    ) %>%
    dplyr::mutate(
      seq_id = glue::glue_data(., "{genome}_{seq_id}") %>% as.character,
      feat_id = glue::glue_data(., "{genome}_{feat_id}") %>% as.character
    ) %>%
    dplyr::distinct() %>%
    as.data.frame()
  
  gs_features <- bind_rows(ref_features, query_features) %>%
    distinct()
  
  message(glue::glue("# Constructed features data in {as.numeric(difftime(now(),start_time, units = 'secs'))  %>% round(2)} seconds"))
  start_time <- now()
  
  syntenic_genes <- dtplyr::lazy_dt(gs_features) %>%
    dplyr::distinct() %>%
    dplyr::group_by(seq_id) %>%
    dplyr::mutate(
      offset = min(start),
      start = start - offset + 1,
      end = end - offset + 1
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-offset) %>%
    dplyr::left_join(pangenes, by = c("genome", "feat_id")) %>%
    as.data.frame()
  
  message(glue::glue("# Constructed syntenic genes in {as.numeric(difftime(now(),start_time, units = 'secs'))  %>% round(2)} seconds"))
  start_time <- now()
  ## Annotations
  
  message("# Subsetting annotations")
  start_time <- now()
  roi_annotations <-
    gs_features %>%
    dplyr::distinct() %>%
    {
      parallel::mclapply(unique(.$seq_id),
                         \(contig) {
                           annotations %>%
                             dplyr::filter(
                               seq_id == contig,
                               start >= min(gs_features[gs_features$seq_id == contig, ]$start) - 1,
                               end <= max(gs_features[gs_features$seq_id == contig, ]$end) +
                                 1
                             ) %>%
                             dplyr::mutate(
                               offset = min(start),
                               start = start - offset + 1,
                               end = end - offset + 1
                             )
                         },
                         mc.cores = ncores)
    } %>%
    dplyr::bind_rows() %>%
    dplyr::left_join(syntenic_genes %>%
                       dplyr::select(seq_id, feat_id, og),
                     by = c("seq_id", "feat_id")) %>%
    dplyr::select(seq_id,
                  start,
                  end,
                  strand,
                  type,
                  feat_id,
                  name,
                  og,
                  bin_id,
                  type,
                  introns) %>%
    distinct()
  message(glue::glue("# Finished subsetting annotations in {as.numeric(difftime(now(),start_time, units = 'secs')) %>% round(2)} seconds"))
  
  message("# Creating sequence table")
  start_time <- now()
  
  seqs <- dtplyr::lazy_dt(syntenic_genes) %>%
    dplyr::select(bin_id = genome, seq_id, start, end) %>%
    group_by(seq_id) %>%
    summarize(
      length = as.integer(max(end) - min(start)),
      bin_id = unique(bin_id),
      .groups = "drop"
    ) %>%
    distinct() %>%
    as.data.frame()
  
  message(glue::glue("# Finished creating sequences in {as.numeric(difftime(now(),start_time, units = 'secs')) %>% round(2)} seconds"))
  
  message("# Creating ranges")
  start_time <- now()
  # Take the combined output
  
  gs_ranges <- rbind(
    dat %>%
      dplyr::select(
        seq_id = chr1,
        start = start1,
        end = end1,
        feat_id = id1,
        genome = genome1
      ),
    dat %>%
      dplyr::select(
        seq_id = chr2,
        start = start2,
        end = end2,
        feat_id = id2,
        genome = genome2
      )
  ) %>%
    dtplyr::lazy_dt() %>%
    distinct() %>%
    # Get the start and end of the sequence
    group_by(genome, seq_id) %>%
    summarize(
      start2 = min(c(start, end)) - pad,
      end2 = max(c(start, end)) + pad,
      .groups = "drop"
    ) %>%
    dplyr::select(
      genome = genome,
      chr = seq_id,
      start = start2,
      end = end2
    ) %>%
    as.data.frame() %>%
    mutate(seq_id = glue::glue("{genome}_{chr}"))
  message(glue::glue("# Finished creating ranges in {as.numeric(difftime(now(),start_time, units = 'secs')) %>% round(2)} seconds"))
  
  # Run nested lapply.
  #  Outer runs genespace per range from initial genome
  #  Inner subsets results to only keep those that fall into initial range
  
  message("# Collecting synteny for ranges")
  start_time <- now()
  
  gs_ranged_links <- lapply(unique(gs_ranges$seq_id),
                            \(contig) {
                              tmp_bed <- gs_ranges %>%
                                dplyr::filter(seq_id == contig) %>%
                                dplyr::select(-seq_id)
                              
                              tmp_hits <-
                                GENESPACE::query_hits(genespace_output,
                                                      bed = tmp_bed,
                                                      synOnly = T)[[1]] %>%
                                as.data.frame()
                              
                              # Subset hits to only keep those that fall into the pre-existing intervals.
                              parallel::mclapply(unique(gs_ranges$seq_id)[unique(gs_ranges$seq_id) != contig],
                                                 \(hit_contig) {
                                                   # Create subset of target ranges
                                                   hit_range <-
                                                     gs_ranges[gs_ranges$seq_id == hit_contig,]
                                                   
                                                   # Filter tmp_hits for those in target range
                                                   tmp_hits %>%
                                                     dplyr::filter(
                                                       genome2 == unique(hit_range$genome) &
                                                         chr2 == unique(hit_range$chr) &
                                                         start2 >= min(hit_range$start) &
                                                         end2 <= max(hit_range$end)
                                                     )
                                                 },
                                                 mc.cores = ncores) %>%
                                bind_rows
                            }) %>%
    bind_rows()
  
  message(glue::glue("# Updated synteny in {as.numeric(difftime(now(),start_time, units = 'secs')) %>% round(2)} seconds"))
  
  message("# Converting synteny to links")
  start_time <- now()
  # This code converts output of query_genes to a table that can be used with gggenomes geom_link
  links <- dtplyr::lazy_dt(gs_ranged_links) %>%
    dplyr::mutate(id = 1:nrow(.)) %>%
    dplyr::mutate(
      seq_id1 = glue::glue_data(., "{genome1}_{chr1}") %>% as.character(),
      seq_id2 = glue::glue_data(., "{genome2}_{chr2}") %>% as.character()
    ) %>%
    tidyr::pivot_longer(starts_with("seq_id"),
                        values_to = "seq_id",
                        names_to = "seq_id_type") %>%
    # Compute offset per seq_id
    dplyr::group_by(seq_id) %>%
    dplyr::mutate(offset = min(start1[seq_id_type == "seq_id1"],
                               start2[seq_id_type == "seq_id2"]) - 1) %>%
    # For each comparison, adjust the starts and ends
    dplyr::group_by(id) %>%
    dplyr::mutate(
      start1 = start1 - offset[seq_id_type == "seq_id1"],
      end1 = end1 - offset[seq_id_type == "seq_id1"] - 1,
      start2 = start2 - offset[seq_id_type == "seq_id2"],
      end2 = end2 - offset[seq_id_type == "seq_id2"] - 1,
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-offset, -id) %>%
    # Pivot back to wide
    tidyr::pivot_wider(names_from = seq_id_type, values_from = seq_id) %>%
    dplyr::left_join(seqs %>%
                       dplyr::select(seq_id1 = seq_id, length1 = length), by = "seq_id1") %>%
    dplyr::left_join(seqs %>%
                       dplyr::select(seq_id2 = seq_id, length2 = length),
                     by = "seq_id2") %>%
    dplyr::select(
      seq_id = seq_id1,
      length = length1,
      start = start1,
      end = end1,
      seq_id2,
      length2,
      start2,
      end2
    ) %>%
    # Dereplicate table, keeping only one entry per genome-start-end pair
    # This is done by sorting the seq_id-start-end strings and
    # keeping only one per combination
    as_tibble() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(distinction =
                    paste(sort(c(
                      paste(seq_id, start, end, sep = "_"),
                      paste(seq_id2, start2, end2, sep = "_")
                    )),
                    collapse = "_")) %>%
    dplyr::group_by(distinction) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-distinction)
  
  message(glue::glue("# Converted synteny to links in {as.numeric(difftime(now(),start_time, units = 'secs')) %>% round(2)} seconds"))
  
  message(glue::glue("# Finished. (Total time: {as.numeric(difftime(now(),start_fn_time, units = 'secs')) %>% round(2)} seconds)"))
  
  list(
    genes = as_tibble(roi_annotations),
    seqs = as_tibble(seqs),
    links = as_tibble(links),
    syntenic_genes = as_tibble(syntenic_genes),
    gs_features = as_tibble(gs_features),
    query = as_tibble(query_features),
    ref = as_tibble(ref_features),
    dat = as_tibble(dat)
  )
}


#' Prepare genespace output for synteny analysis
#'
#' @param genespace_output genespace output object
#' @param bed_table a bed file containing a region of interest in the "reference"
#' @param genomes a list of genomes (should be identical to those used for genespace)
#' @param padding padding to add around regions to avoid cutting
#' @param verbose Some extra messages for debugging
#'
#' This functions loops through genomes to obtain all gene-vs-gene relationships for the region of interest
#'
#' @return A data frame containing the all-vs-all gene level synteny relationships.
#' @export
#'
#' @examples
#'
#' @author Niklas Schandry
#'
prepare_genespace_synteny <- function(genespace_output,
                                      bed_table,
                                      genomes,
                                      padding = 0,
                                      verbose = 0) {
  ref_hits <- GENESPACE::query_hits(genespace_output,
                                    bed = bed_table,
                                    synOnly = T)
  
  # Get reference hits
  ref_hits <- ref_hits[[1]] %>%
    filter(!genome1 == genome2) %>%
    filter(str_detect(chr1, bed_table$chr) &
             str_detect(chr2, bed_table$chr))
  if (verbose > 0) {
    message(glue::glue("Obtained reference hits"))
  }
  result <- ref_hits
  
  # Remove reference from genomes
  genomes <- genomes[!genomes == bed_table$genome]
  
  extra_hits <- lapply(genomes, \(cur_genome) {
    if (verbose > 0) {
      message(glue::glue("Obtaining hits for {cur_genome}"))
    }
    # Use the reference hits to construct a new bed file for the first genome
    tmp_bed <- dtplyr::lazy_dt(ref_hits) %>%
      dplyr::filter(genome2 == cur_genome) %>%
      dplyr::select(
        genome = genome2,
        chr = chr2,
        start = start2,
        end = end2
      ) %>%
      dplyr::group_by(genome, chr) %>%
      # Add padding around start and end of bed file
      dplyr::summarize(
        start2 = min(c(start, end)) - padding,
        end2 = max(c(start, end)) + padding,
        .groups = "drop"
      ) %>%
      dplyr::select(
        genome = genome,
        chr = chr,
        start = start2,
        end = end2
      ) %>%
      as.data.frame()
    if (verbose > 0) {
      message(glue::glue("Created bed for {cur_genome}"))
    }
    if (nrow(tmp_bed) == 0) {
      message(glue::glue("{cur_genome} returned empty bed file"))
    }
    # Extract all syntenic hits, if there are rows in the tmp_bed
    if (nrow(tmp_bed) > 0) {
      if (verbose > 0) {
        message(glue::glue("Obtaining hits for {cur_genome}"))
      }
      tmp_hits <- GENESPACE::query_hits(genespace_output,
                                        bed = tmp_bed,
                                        synOnly = T)
      tmp_hits <- tmp_hits[[1]] %>%
        dplyr::filter(!genome1 == genome2) %>%
        dplyr::filter(str_detect(chr1, bed_table$chr) &
                        str_detect(chr2, bed_table$chr))
      if (verbose > 0) {
        message(glue::glue("Obtained {nrow(tmp_hits)} hits for {cur_genome}"))
      }
      tmp_hits
    }
    
  })
  result <- bind_rows(ref_hits, bind_rows(extra_hits)) %>%
    dplyr::filter(str_detect(chr1, bed_table$chr) &
                    str_detect(chr2, bed_table$chr))
  return(result)
}
