library(targets)
library(tidySingleCellExperiment)
library(purrr)
library(dplyr)
library(stringr)
library(Seurat)

# Read and save each sample individually to disk
x = readRDS("/vast/projects/thymic_survive/shen.m_M000637_Gray_batch1/extdata/SCEs/M000637_Gray_batch1.demultiplexed.SCE.rds")
colData = as.data.frame(colData(x))
rna_counts = counts(x)
x = CreateSeuratObject(rna_counts, meta.data=colData)
x |> dplyr::count(Sample)
seu_list <- x |> group_split(Sample)
outdir <- "/vast/projects/thymic_survive/results"
if (!dir.exists(outdir)) dir.create(outdir)
map(seu_list, ~ {
  sample = .x |> distinct(Sample) |> pull()
  saveRDS(.x, file.path(outdir, paste0(sample, "_seurat.rds")))
}, .progress = T)

# QC by HPCell
samples = list.files("/vast/projects/thymic_survive/results", full.names = T, pattern = "_seurat.rds")
samples <- samples |>
  setNames(
    samples |>
      basename() |>
      sub("_seurat\\.rds$", "", x = _)
  )
target_store = "/vast/scratch/users/shen.m/thymic_survive/"
job::job({
  library(HPCell)
  library(tarchetypes)
  samples |>
    
    initialise_hpc(
      store = target_store,
      gene_nomenclature = "ensembl",
      data_container_type = "seurat_rds",
      species = "mouse",
      computing_resources =  crew.cluster::crew_controller_slurm(
        name = "elastic",
        workers = 300,
        tasks_max = 20,
        seconds_idle = 30,
        crashes_error = 10,
        options_cluster = crew.cluster::crew_options_slurm(
          memory_gigabytes_required = c(10, 20, 40),
          cpus_per_task = c(2),
          time_minutes = c(60*4),
          verbose = TRUE,
          script_lines = c(
            "#!/bin/bash",
            "module purge",
            "module load quarto"
          )
        )
        
      ),
      verbosity = "summary",
      update = "thorough", 
      error = "continue",
      garbage_collection = 100, 
      workspace_on_error = TRUE
      
    ) |> 
    # Empty droplet
    remove_empty_DropletUtils() |> 
    
    # Annotation
    annotate_cell_type() |> 
    
    # Alive
    remove_dead_scuttle(target_annotation = "annotation_tbl", 
                        group_by = "mouserna.labels.coarse") |>
    
    # Doublet
    remove_doublets_scDblFinder() |>
    
    
    print()
})

# debugonce(alive_identification)
# alive_identification(data_object, empty_tbl, annotation_tbl, cell_type_column = "mouserna.labels.coarse",
#                      feature_nomenclature = gene_nomenclature, species_db = species)
tar_workspace(empty_report, store = target_store)
data_object <- map(data_object, ~{
  assay_name = .x@assays |> names() |> extract2(1)
  counts = assay(.x, assay_name)
  coldata <- as.data.frame(colData(.x))
  .y = CreateSeuratObject(counts, meta.data = coldata, assay = assay_name)
})


# Load quality control results from HPC pipeline
sample_names <- tar_read(sample_names, store = target_store) 
data_object <- tar_read(data_object, store = target_store) 
empty_tbl <- tar_read(empty_tbl, store = target_store) 
alive_tbl <- tar_read(alive_tbl,store = target_store) 
annotation_tbl <- tar_read(annotation_tbl, store = target_store) 
doublet_tbl <- tar_read(doublet_tbl, store = target_store) 


# ## Test Empty Droplets report
# 
# rmarkdown::render(
#   input =  file.path(rprojroot::find_package_root_file(), "inst/rmd/Empty_droplet_report.qmd"),
#   output_file = file.path(rprojroot::find_package_root_file(), "inst/rmd/Empty_droplet_report.html"),
#   params = list(x1 = tar_read(empty_tbl, store = target_store),
#                 x2 = tar_read(data_object, store = target_store),
#                 x3 = tar_read(alive_tbl, store = target_store),
#                 x4 = sample_names
# ))

# Apply QC filters to the data
data = readRDS("/vast/projects/thymic_survive/shen.m_M000637_Gray_batch1/extdata/SCEs/M000637_Gray_batch1.demultiplexed.SCE.rds")
data <- data |> left_join(empty_tbl |>bind_rows())
data <- data |> left_join(alive_tbl |> select(.cell, alive))
data <- data |> left_join(doublet_tbl)

# Filter for high-quality cells
qc_data <- data |> 
  filter(empty_droplet == FALSE,
         alive == TRUE,
         scDblFinder.class != "doublet")

cat("Original cell count:", ncol(data), "\n")
cat("Cells after QC filtering:", ncol(qc_data), "\n")
cat("QC retention rate:", round(ncol(qc_data)/ncol(data)*100, 2), "%\n")


# Check non-empty droplet rates per sample
data = readRDS("/vast/projects/thymic_survive/shen.m_M000637_Gray_batch1/extdata/SCEs/M000637_Gray_batch1.demultiplexed.SCE.rds")
data <- data |> left_join(empty_tbl |>bind_rows())
qc_tibble <- data |>
  summarise(
    total_cells = n(),
    empty_cells = sum(empty_droplet),
    empty_pct = empty_cells / total_cells * 100,
    .by = Sample
  ) |> arrange(Sample)

