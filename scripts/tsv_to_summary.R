

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript tsv_to_summary.R <tsv_file> <output_dir>")
}
tsv_file <- args[1]
output_dir <- args[2]


summarize_with_threshold <- function(file_path) {
  library(data.table)
  library(fs)
  library(glue)
  library(stringr)
  
  thresholds_list <- list(
    MAPQ_min = 1.69,
    AS_min = 15.6,
    XS_max = 29.0,
    NM_max = 5.70
  )

  summary_output_folder <- output_dir
  log_output_folder <- output_dir

  base_name <- fs::path_file(file_path)
  file_sample_name <- sub("\\.tsv$", "", base_name)

  tryCatch({
    cat(paste("Processing:", base_name, "...\n"))

    # --- A. Define Output Paths ---
    summary_out_path <- fs::path(summary_output_folder, glue("{file_sample_name}_summary.tsv"))
    log_out_path <- fs::path(log_output_folder, glue("{file_sample_name}_log.tsv"))

    # --- B. Read File Efficiently ---

    all_lines <- readLines(file_path, warn = FALSE)
    is_summary_line <- str_starts(all_lines, "UNMAPPED_SUMMARY")

    unmapped_count <- 0
    if (any(is_summary_line)) {
      unmapped_line <- all_lines[which(is_summary_line)[1]]
      unmapped_count_str <- str_extract(unmapped_line, "Total unmapped reads: \\d+") |> str_extract("\\d+")
      if (!is.na(unmapped_count_str)) unmapped_count <- as.numeric(unmapped_count_str)
    }

    # Keep only data lines (non-summary)
    data_lines <- all_lines[!is_summary_line]
    rm(all_lines)  # free memory

    # --- B.2. Check if any data lines exist ---
    if (length(data_lines) == 0) {
      
      cat(paste("  > WARNING:", base_name, "contained 0 data lines. Proceeding with 0 counts.\n"))
      total_primary_reads <- 0
      total_passed_reads <- 0
      num_filtered_out <- 0
      protein_summary <- data.table(rname = character(), read_count = integer())
      rm(data_lines)
      
    } else {
      
      # --- B.3. Write/Read Data Lines ---
      # Write temporary subset to disk and read directly via fread (fast + low RAM)
      tmp_file <- tempfile(fileext = ".tsv")
      writeLines(data_lines, tmp_file)
      rm(data_lines)

      # fread reads directly into memory-efficient data.table
      dt <- fread(tmp_file, sep = "\t", quote = "", showProgress = FALSE)
      unlink(tmp_file)

      # --- C. Check for valid data and columns ---
      required_cols <- c("is_primary", "mapq", "AS", "XS", "NM", "rname")
      
      if (nrow(dt) > 0 && all(required_cols %in% names(dt))) {
        
        # --- C.1. Filter for primary reads ---
        dt <- dt[is_primary == TRUE]
        total_primary_reads <- nrow(dt)

        if (total_primary_reads > 0) {
          # Coerce numeric columns (fast in-place conversion)
          numeric_cols <- c("mapq", "AS", "XS", "NM")
          for (col in numeric_cols) set(dt, j = col, value = as.numeric(dt[[col]]))

          # --- D. Apply thresholds efficiently ---
          dt_filt <- dt[
            mapq >= thresholds_list$MAPQ_min &
            AS   >= thresholds_list$AS_min &
            XS   <= thresholds_list$XS_max &
            NM   <= thresholds_list$NM_max
          ]

          total_passed_reads <- nrow(dt_filt)
          num_filtered_out <- total_primary_reads - total_passed_reads

          # --- E. Summarize counts per protein ---
          protein_summary <- dt_filt[, .(read_count = .N), by = rname][order(-read_count)]
          
          rm(dt_filt) # Cleanup

        } else {
          # Case: File had data, but 0 primary reads
          cat(paste("  > NOTE:", base_name, "contained 0 primary reads.\n"))
          total_passed_reads <- 0
          num_filtered_out <- 0
          protein_summary <- data.table(rname = character(), read_count = integer())
        }
        
        rm(dt) # Cleanup

      } else {
        # --- C/D/E. Handle empty or malformed file (e.g., "x", "no_alignments") ---
        cat(paste("  > WARNING:", base_name, "contained 0 valid data rows or was malformed. Proceeding with 0 counts.\n"))
        total_primary_reads <- 0
        total_passed_reads <- 0
        num_filtered_out <- 0
        protein_summary <- data.table(rname = character(), read_count = integer()) # Empty summary
        
        if (exists("dt")) rm(dt) # Cleanup
      }
    } # End of data processing logic

    # --- F. Write outputs ---
    fwrite(protein_summary, summary_out_path, sep = "\t", quote = FALSE)
    
    log_data <- data.table(
      file_name = base_name,
      processed_time = as.character(Sys.time()),
      total_primary_reads = total_primary_reads,
      filtered_out_thresholds = num_filtered_out,
      passed_summarized = total_passed_reads,
      unmapped_reads = unmapped_count
    )
    fwrite(log_data, log_out_path, sep = "\t", quote = FALSE)

    # Delete input after success
    fs::file_delete(file_path)

    cat(paste("  > Finished:", base_name, ". Summary and log saved.\n"))

    # Cleanup
    if (exists("protein_summary")) rm(protein_summary)
    if (exists("log_data")) rm(log_data)
    gc()

  }, error = function(e) {
    cat(paste("!!! ERROR processing", base_name, ":", e$message, "\n"))
    cat("  > File was NOT deleted.\n")
    
    # Clean up any partial objects in case of error
    if (exists("dt")) rm(dt)
    if (exists("dt_filt")) rm(dt_filt)
    if (exists("protein_summary")) rm(protein_summary)
    if (exists("log_data")) rm(log_data)
    if (exists("data_lines")) rm(data_lines)
    gc()
  })
}

summarize_with_threshold(tsv_file)