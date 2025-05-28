#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(readr)
library(glue)
library(Rsamtools)
library(tools)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript bam_to_tsv.R <bam_file> <output_dir>")
}
bam_file <- args[1]
output_dir <- args[2]

# Create output directory if it doesn't exist
create_output_dir <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
}

# Prepare and index BAM file if needed
prepare_bam_file <- function(bam_file) {
  sorted_bam <- sub("\\.bam$", "_sorted.bam", bam_file)
  if (!file.exists(sorted_bam)) {
    system(paste("samtools sort -o", shQuote(sorted_bam), shQuote(bam_file)))
  }
  if (!file.exists(paste0(sorted_bam, ".bai"))) {
    indexBam(sorted_bam)
  }
  return(sorted_bam)
}

# Normalize lengths of columns from scanBam
normalize_length <- function(x, target_len) {
  if (is.null(x)) {
    return(rep(NA, target_len))
  } else if (length(x) < target_len) {
    return(c(x, rep(NA, target_len - length(x))))
  } else {
    return(head(x, target_len))
  }
}

process_bam_file <- function(bam_file, output_dir, mapq_threshold = 0) {
  bam <- Rsamtools::BamFile(bam_file)
  count <- Rsamtools::countBam(bam)$records
  if (count == 0) return("no_alignments")

  param <- Rsamtools::ScanBamParam(
    what = c("qname", "flag", "rname", "strand", "pos", "mapq", "cigar", "seq", "qual"),
    tag = c("AS", "XS", "NM", "MD")
  )
  bam_data <- Rsamtools::scanBam(bam, param = param)[[1]]

  # Find actual number of alignments
  lengths <- sapply(bam_data, function(x) if (is.null(x)) 0 else length(x))
  actual_length <- max(lengths)

  # Build the data frame
  bam_df <- data.frame(
    qname = normalize_length(bam_data$qname, actual_length),
    flag  = normalize_length(bam_data$flag, actual_length),
    rname = normalize_length(bam_data$rname, actual_length),
    strand = normalize_length(bam_data$strand, actual_length),
    pos = normalize_length(bam_data$pos, actual_length),
    mapq = normalize_length(bam_data$mapq, actual_length),
    cigar = normalize_length(bam_data$cigar, actual_length),
    seq = normalize_length(bam_data$seq, actual_length),
    qual = normalize_length(bam_data$qual, actual_length),
    AS = normalize_length(bam_data$tag$AS, actual_length),
    XS = normalize_length(bam_data$tag$XS, actual_length),
    NM = normalize_length(bam_data$tag$NM, actual_length),
    MD = normalize_length(bam_data$tag$MD, actual_length),
    stringsAsFactors = FALSE
  )

  bam_df$is_secondary <- bitwAnd(bam_df$flag, 256) == 256
  bam_df$is_primary <- !bam_df$is_secondary
  bam_df$unmapped <- bitwAnd(bam_df$flag, 0x4) != 0

  # Tally unmapped reads
  unmapped_count <- sum(bam_df$unmapped, na.rm = TRUE)

  filtered_df <- bam_df[
  (!bam_df$unmapped) & (!is.na(bam_df$mapq)) & (bam_df$mapq >= mapq_threshold),
  ]

  # Add a summary row with unmapped count
  if (unmapped_count > 0) {
    if (nrow(filtered_df) == 0) {
      # create empty dataframe with correct columns + note
      filtered_df <- data.frame(
        qname = character(),
        flag = integer(),
        rname = character(),
        strand = character(),
        pos = integer(),
        mapq = numeric(),
        cigar = character(),
        seq = character(),
        qual = character(),
        AS = integer(),
        XS = integer(),
        NM = integer(),
        MD = character(),
        is_secondary = logical(),
        is_primary = logical(),
        unmapped = logical(),
        note = character(),
        stringsAsFactors = FALSE
      )
    } else {
      filtered_df$note <- NA_character_
    }

    summary_row <- data.frame(
      qname = "UNMAPPED_SUMMARY",
      flag = 4,
      rname = NA,
      strand = NA,
      pos = NA,
      mapq = NA,
      cigar = NA,
      seq = NA,
      qual = NA,
      AS = NA,
      XS = NA,
      NM = NA,
      MD = NA,
      is_secondary = NA,
      is_primary = NA,
      unmapped = TRUE,
      note = paste("Total unmapped reads:", unmapped_count),
      stringsAsFactors = FALSE
    )

    filtered_df <- rbind(filtered_df, summary_row)
  }
  df <- filtered_df %>%
    mutate(
      flag = as.integer(flag),
      mapq = as.numeric(mapq),
      is_secondary = is_secondary == "TRUE",
      is_primary = is_primary == "TRUE"
    )

  df_old <- df

  # Adding primary read sequence and protein id list to full sv
  df_primary <- df %>% 
    dplyr::select(qname, rname, is_primary, seq) %>%
    dplyr::filter(is_primary) %>%
    dplyr::rename(primary_rname = rname) %>%
    dplyr::select(-is_primary)

  df_final <- df_old %>%
    select(-seq) %>% 
    left_join(df_primary)
  return(df_final)
}

# Main pipeline runner
run_bam_to_tsv_pipeline <- function(bam_file, output_dir, mapq_threshold = 0) {
  create_output_dir(output_dir)
  sorted_bam <- prepare_bam_file(bam_file)
  processed_df <- process_bam_file(sorted_bam, output_dir, mapq_threshold)
  out_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(bam_file)), ".tsv"))
  utils::write.table(processed_df, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

run_bam_to_tsv_pipeline(bam_file, output_dir, mapq_threshold = 0)