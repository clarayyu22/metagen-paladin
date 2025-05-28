library(dplyr)
library(ggplot2)
library(readr)
library(glue)


# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript merge_and_summarize.R <input_directory> <output_summary_file>")
}
input_dir <- args[1]
output_dir <- args[2]

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load data
# Find all _filtered.tsv files recursively
tsv_files <- list.files(path = input_dir, pattern = ".tsv$", recursive = TRUE, full.names = TRUE)
if (length(tsv_files) == 0) {
  stop("No TSV files found in the specified directory.")
}

# Load and combine all TSV files
df_list <- lapply(tsv_files, function(file) {
  read_tsv(file, col_types = cols(.default = "c"))
})
df <- bind_rows(df_list)

df <- df %>%
  mutate(
    flag = as.integer(flag),
    mapq = as.numeric(mapq),
    is_secondary = is_secondary == "TRUE",
    is_primary = is_primary == "TRUE"
  )

# Extract and sum unmapped counts from UNMAPPED_SUMMARY rows
unmapped_total <- df %>%
  filter(qname == "UNMAPPED_SUMMARY") %>%
  pull(note) %>%
  gsub("Total unmapped reads: ", "", .) %>%
  as.integer() %>%
  sum(na.rm = TRUE)

# Create combined UNMAPPED_SUMMARY row with count and label in rname
unmapped_row <- tibble(
  rname = "Total unmapped reads",
  is_primary = NA,
  count = unmapped_total,
  primary_rname_list = NA
)

primary_map_df <- df %>% 
  filter(is_secondary) %>% 
  select(rname, primary_rname) %>% 
  group_by(rname) %>% 
  dplyr::summarise(
    primary_rname_list = paste(unique(primary_rname), collapse = ","),
    .groups = "drop"
  ) %>%
  mutate(is_secondary = TRUE)

# What primary reads are above the mapq threshold?
quality_qnames <- df %>% 
  filter(is_primary) %>%
  filter(mapq >= 10) %>%
  pull(qname)

summary_stats <- df %>%
  filter(!is.na(rname)) %>%
  filter(qname  %in% quality_qnames) %>% 
  # filter(mapq >= 10) %>%
  group_by(rname, is_primary, is_secondary) %>%
  summarise(count = n(), .groups = "drop") %>%
  # group_by(is_primary) %>%
  # mutate(rel_abundance = count / sum(count)) %>% 
  left_join(primary_map_df)

summary_stats$count %>% summary

summary_final <- summary_stats %>%
  select(-is_secondary) %>%
  bind_rows(unmapped_row)

output_file <- file.path(output_dir, "summary.tsv")
write.table(summary_final, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
cat("Summary relative abundance statistics written to", output_file, "\n")