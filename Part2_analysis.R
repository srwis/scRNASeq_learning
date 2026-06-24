# Single-Cell RNA-Seq Analysis - Part 2: Remaining 5 Compartments
# Loops through all remaining compartments, runs patient-specific filtering, V5 normalization,
# scaling, PCA, Harmony integration, SNN clustering, and UMAP. Saves RDS objects and UMAP plots.

message("Loading packages...")
library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(purrr)
library(patchwork)

# 1. Configuration of compartments and paths
expression_dir <- "data/SCP1884/expression"

compartments <- list(
  "CO_EPI" = list(path = file.path(expression_dir, "62a79393d8bced7ddefbf0d1"), name = "Colon_Epithelial"),
  "CO_IMM" = list(path = file.path(expression_dir, "62a7a911a54b79c09baa336a"), name = "Colon_Immune"),
  "TI_STR" = list(path = file.path(expression_dir, "62a7b26f0e85c7e6d0bb03f0"), name = "Ileum_Stromal"),
  "TI_EPI" = list(path = file.path(expression_dir, "62a7b6bf52be218275f15f43"), name = "Ileum_Epithelial"),
  "TI_IMM" = list(path = file.path(expression_dir, "62a7c0fd3c8dbbf858d3ac47"), name = "Ileum_Immune")
)

# Robust baseline QC thresholds per compartment
compartment_qc <- list(
  "CO_EPI" = list(min_feat = 200, max_feat = 4000, max_count = 15000, max_mt = 5),
  "CO_IMM" = list(min_feat = 200, max_feat = 2500, max_count = 10000, max_mt = 5),
  "TI_STR" = list(min_feat = 200, max_feat = 3000, max_count = 12000, max_mt = 5),
  "TI_EPI" = list(min_feat = 200, max_feat = 4000, max_count = 15000, max_mt = 5),
  "TI_IMM" = list(min_feat = 200, max_feat = 2500, max_count = 10000, max_mt = 5)
)

# Patient overrides (empty by default, structured for custom tuning)
patient_overrides <- list()

# Helper function to read raw matrices into Seurat objects (exactly as in Part 1)
create_seurat_obj <- function(path, comp_name) {
  matrix_path <- list.files(path, ".*raw.mtx", full.names = TRUE)
  barcodes_path <- list.files(path, ".*.barcodes.tsv", full.names = TRUE)
  features_path <-  list.files(path, ".*.features.tsv", full.names = TRUE)
  
  if (length(matrix_path) == 0 || length(barcodes_path) == 0 || length(features_path) == 0) {
    stop(paste("Required raw data files not found in", path))
  }
  
  sample_name <- matrix_path |> basename() |> substr(1, 6)
  
  message("Reading matrix files...")
  counts_matrix <- ReadMtx(
    mtx = matrix_path,
    cells = barcodes_path,
    features = features_path,
    cell.column = 1,
    feature.column = 2
  )
  
  obj <- CreateSeuratObject(counts = counts_matrix, project = sample_name)
  obj$Tissue <- sample_name
  obj$Compartment <- comp_name
  return(obj)
}

# 2. Main Sequential Loop
for (comp_id in names(compartments)) {
  comp_info <- compartments[[comp_id]]
  comp_path <- comp_info$path
  comp_name <- comp_info$name
  
  message("\n=========================================")
  message("PROCESSING COMPARTMENT: ", comp_name, " (", comp_id, ")")
  message("=========================================")
  
  # Step A: Load Raw Data
  message("Loading raw Seurat object...")
  se_raw <- create_seurat_obj(comp_path, comp_name)
  message("Loaded cells: ", ncol(se_raw), ", features: ", nrow(se_raw))
  
  # Step B: Quality Control by Patient
  message("Splitting object by patient identity...")
  patient_list <- SplitObject(se_raw, split.by = "orig.ident")
  
  message("Applying quality control filters per patient...")
  filtered_list <- imap(patient_list, function(obj, patient_id) {
    # Default thresholds for this compartment
    thresh <- compartment_qc[[comp_id]]
    
    # Check for patient-specific overrides
    if (!is.null(patient_overrides[[patient_id]])) {
      thresh <- patient_overrides[[patient_id]]
    }
    
    # Calculate mitochondrial percentage
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    
    # Identify passing cells
    keep_cells <- colnames(obj)[
      obj$nFeature_RNA > thresh$min_feat & 
      obj$nFeature_RNA < thresh$max_feat & 
      obj$nCount_RNA < thresh$max_count &
      obj$percent.mt < thresh$max_mt
    ]
    
    if (length(keep_cells) == 0) {
      warning(paste("Patient", patient_id, "has zero cells passing filters, skipping."))
      return(NULL)
    }
    
    # Subset
    obj <- subset(obj, cells = keep_cells)
    return(obj)
  })
  
  # Remove NULL elements
  filtered_list <- compact(filtered_list)
  
  if (length(filtered_list) == 0) {
    stop("No cells passed QC thresholds across all patient samples!")
  }
  
  message("Merging filtered patient samples back together...")
  if (length(filtered_list) == 1) {
    se_clean <- filtered_list[[1]]
  } else {
    se_clean <- merge(
      x = filtered_list[[1]],
      y = filtered_list[-1],
      project = paste0(comp_id, "_Clean")
    )
  }
  
  message("Original cell count: ", ncol(se_raw))
  message("Filtered cell count: ", ncol(se_clean))
  
  # Step C: Seurat normalization, feature selection, scaling
  message("Normalizing and finding variable features...")
  se_clean <- NormalizeData(se_clean, normalization.method = "LogNormalize", scale.factor = 10000)
  se_clean <- FindVariableFeatures(se_clean, selection.method = "vst", nfeatures = 2000)
  
  message("Scaling data...")
  se_clean <- ScaleData(se_clean)
  
  message("Running PCA...")
  se_clean <- RunPCA(se_clean, npcs = 30, verbose = FALSE)
  
  # Step D: Harmony batch integration
  message("Running Harmony integration across patient samples...")
  se_clean <- IntegrateLayers(
    object = se_clean,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    verbose = FALSE
  )
  
  # Step E: Clustering and UMAP
  message("Finding neighbors and clusters (resolution = 0.5)...")
  se_clean <- FindNeighbors(se_clean, reduction = "harmony", dims = 1:20)
  se_clean <- FindClusters(se_clean, resolution = 0.5, verbose = FALSE)
  
  message("Running UMAP embedding...")
  se_clean <- RunUMAP(se_clean, reduction = "harmony", dims = 1:20)
  
  # Step F: Generate Plots
  message("Generating integrated UMAP plots...")
  p1 <- DimPlot(se_clean, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
    ggtitle(paste(comp_name, "- Clusters")) +
    theme_minimal()
    
  p2 <- DimPlot(se_clean, reduction = "umap", group.by = "orig.ident", label = FALSE) +
    ggtitle(paste(comp_name, "- Patient Mixing")) +
    theme_minimal() +
    theme(legend.position = "none") # Hide large patient legend
    
  combined_plot <- p1 + p2
  
  # Save plots
  plot_file <- file.path("data", paste0("umap_", tolower(comp_id), ".png"))
  ggsave(plot_file, plot = combined_plot, width = 12, height = 5, dpi = 150)
  message("Saved UMAP plot to: ", plot_file)
  
  # Step G: Save Seurat Object
  rds_file <- file.path("data", paste0("integrated_", tolower(comp_id), ".rds"))
  message("Saving Seurat object to RDS...")
  saveRDS(se_clean, file = rds_file)
  message("Saved RDS object to: ", rds_file)
  
  # Step H: Clean up memory
  message("Cleaning up memory for next compartment...")
  rm(se_raw, patient_list, filtered_list, se_clean, p1, p2, combined_plot)
  gc()
}

message("\nAll compartments processed successfully!")
