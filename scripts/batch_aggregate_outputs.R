library(glue)
library(batchtools)
library(dplyr)
library(tibble)
source("notebooks/R_scripts/aggregate_outputs.R")
source("notebooks/R_scripts/helpers_general.R")
wkdir <- "/central/groups/MazmanianLab/clara/microbiome_sample_search"
metadata_file <- "/central/scratch/clarayu/sample_meta_full_2025-05-20.tsv"
tsv_dir <- "/resnick/groups/MazmanianLab/clarayu/full_run/fastqs/tsv_outputs"


# ---- Get TSV files ----
tsv_files <- base::list.files(tsv_dir, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)

length(tsv_files)/40

tsv_chunks <- chunk_func(tsv_files, 48)

duckdb_path <- glue("{wkdir}/data/interim/duckdb")
dir.create(duckdb_path, recursive = TRUE, showWarnings = FALSE)

batchtools_params <- tibble(
    "tsv_files" = tsv_chunks
    ) %>% 
    tibble::rownames_to_column(var = "index") %>% 
    mutate(output_csv = glue("{duckdb_path}/chunk_{index}.csv")) %>%
    dplyr::select(-index) %>% 
    glimpse()

batchtools_params %>% View



# ---- Configure registry ----
cluster_run <- glue("{get_time()}_aggregate_outputs")
message("\n\nRUNNING:  ", cluster_run, "\n")
breg <- makeRegistry(
  file.dir = glue(
    "{wkdir}/.cluster_runs/",
    cluster_run
  ),
  seed = 42
)
breg$cluster.functions <- batchtools::makeClusterFunctionsSlurm(
  template = glue("{wkdir}/batch_templates/batchtools.slurm.tmpl"),
  scheduler.latency = 0.1,
  fs.latency = 1
)
# Submit Jobs ----
jobs <- batchMap(
  fun = create_matrix,
  args = batchtools_params,
  reg = breg
)
# jobs[, chunk := chunk(job.id, chunk.size = 20)]
# print(jobs[, .N, by = chunk])

submitJobs(jobs,
  resources = list(
    walltime = 100,
    memory = "100GB",
    ncpus = 2,
    max.concurrent.jobs = 9999
  )
)





