
create_matrix <- function(tsv_files, output_csv) {
    source("notebooks/R_scripts/helpers_general.R")
    # Load libraries
    suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(Matrix)
    library(stringr)
    library(tidyr)
    library(purrr)
    library(glue)
    library(utils)
    })
    base::cat(glue("{get_time()} Starting function\n"))
    # ---- Load metadata ----
    metadata_file <- "/central/scratch/clarayu/sample_meta_full_2025-05-20.tsv"
    metadata <- readr::read_tsv(metadata_file, show_col_types = FALSE) %>%
    dplyr::select(sample, accession) %>%
    dplyr::rename(sample_id = sample, run_id = accession)

    # ---- Process files ----
    protein_sample_counts <- purrr::map_dfr(tsv_files, function(file) {
    if (!is_valid_alignment_file(file)) return(NULL)

    run_id <- extract_run_id(file)
    df <- base::tryCatch(readr::read_tsv(file, col_types = cols(.default = "c"), progress = FALSE), error = function(e) return(NULL))

    if (is.null(df)) return(NULL)

    # Extract unmapped count
    unmapped_count <- df %>%
        dplyr::filter(qname == "UNMAPPED_SUMMARY") %>%
        dplyr::pull(note) %>%
        base::gsub("Total unmapped reads: ", "", .) %>%
        as.integer() %>%
        sum(na.rm = TRUE)

    # Extract primary mappings
    primary_counts <- df %>%
        dplyr::filter(is_primary == "TRUE") %>%
        dplyr::count(rname, name = "abundance")

    # If nothing mapped, we still want to return the unmapped
    out <- dplyr::bind_rows(primary_counts, tibble(rname = "UNMAPPED", abundance = unmapped_count))
    out$run_id <- run_id
    return(out)
    })

    if (base::nrow(protein_sample_counts) == 0) {
    base::stop("❌ No valid alignments found in any TSV files.")
    }

    # ---- Map to sample IDs ----
    joined <- dplyr::left_join(protein_sample_counts, metadata, by = "run_id") %>%
    dplyr::filter(!is.na(sample_id))

    # ---- Create sparse matrix ----
    proteins <- base::unique(joined$rname)
    samples <- base::unique(joined$sample_id)

    # Indices for sparse matrix
    i <- base::match(joined$rname, proteins)
    j <- base::match(joined$sample_id, samples)
    x <- joined$abundance

    mat <- Matrix::sparseMatrix(i = i, j = j, x = x,
                        dims = c(length(proteins), length(samples)),
                        dimnames = list(proteins, samples))

    # ---- Save matrix ----
    dense_df <- as.matrix(mat) %>% as.data.frame()
    utils::write.table(dense_df, file = output_csv,
        sep = ",", quote = FALSE, col.names = NA)
    base::cat(glue("✅ {get_time()} Sparse matrix saved to: {output_csv}\n"))
}

