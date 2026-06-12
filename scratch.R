library(dplyr)
library(Seurat)
library(patchwork)

matrix_path <- "data/SCP1884/expression/62a0c9f523a28233e1c9b06b/CO_STR.scp.raw.mtx"
barcodes_path <- "data/SCP1884/expression/62a0c9f523a28233e1c9b06b/CO_STR.scp.barcodes.tsv"
features_path <-  "data/SCP1884/expression/62a0c9f523a28233e1c9b06b/CO_STR.scp.features.tsv"

counts_matrix <- ReadMtx(
  mtx = matrix_path,
  cells = barcodes_path,
  features = features_path,
  cell.column = 1,     # Usually barcodes are in the first column
  feature.column = 2   # Change to 1 if your gene file only has 1 column (IDs instead of Symbols)
)


#then create Seurat object
co_STR <- CreateSeuratObject(counts = counts_matrix, project = "CO_STR")
co_STR$tissue <- "Colon_Stromal"

matrix_path <- "data/SCP1884/expression/62a79393d8bced7ddefbf0d1/CO_EPI.scp.raw.mtx"
barcodes_path <- "data/SCP1884/expression/62a79393d8bced7ddefbf0d1/CO_EPI.scp.barcodes.tsv"
features_path <-  "data/SCP1884/expression/62a79393d8bced7ddefbf0d1/CO_EPI.scp.features.tsv"

counts_matrix <- ReadMtx(
  mtx = matrix_path,
  cells = barcodes_path,
  features = features_path,
  cell.column = 1,     # Usually barcodes are in the first column
  feature.column = 2   # Change to 1 if your gene file only has 1 column (IDs instead of Symbols)
)


#then create Seurat object
co_STR <- CreateSeuratObject(counts = counts_matrix, project = "CO_STR")
co_STR$tissue <- "Colon_Stromal"

#now let's do QC
#going to just follow the recommendations for now
co_STR[["percent.mt"]] <- PercentageFeatureSet(co_STR, pattern = "^MT-")

#graph some stuff
# In the example below, we visualize QC metrics, and use these to filter cells.

#     We filter cells that have unique feature counts over 2,500 or less than 200
#     We filter cells that have >5% mitochondrial counts

# Visualize QC metrics as a violin plot
VlnPlot(co_STR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(co_STR, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(co_STR, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

co_STR <- subset(co_STR, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


# let write a purrr statement to read all these in

library("purrr")
expression_dir <- "data/SCP1884/expression/"

sample_dirs <- list.dirs(expression_dir)[-1]

this <- purrr::map(sample_dirs, function(path) {
   matrix_path <- list.files(path, ".*raw.mtx", full.names = TRUE)
   barcodes_path <- list.files(path, ".*.barcodes.tsv", full.names = TRUE)
   features_path <-  list.files(path, ".*.features.tsv", full.names = TRUE)
  
   sample_name <- matrix_path |> basename()  |> substr(1,6)

   counts_matrix <- ReadMtx(
     mtx = matrix_path,
     cells = barcodes_path,
     features = features_path,
     cell.column = 1,     # Usually barcodes are in the first column
     feature.column = 2   # Change to 1 if your gene file only has 1 column (IDs instead of Symbols)
   )
  # Create the Seurat object
   obj <- CreateSeuratObject(counts = counts_matrix, project = sample_name)
  
  # Tag with custom metadata
   obj$Tissue <- sample_name
  
   return(obj)
  
})

#and then let's go for the QC step

filtered_seurat_list <- map(this, function(obj) {
  
  # Calculate mitochondrial percentage
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # Apply strict filtering thresholds independently to each object
  # (Note: In a real workflow, you would visualize VlnPlots before setting these hard numbers)
  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  
  return(obj)
})

#and then merge everything

combined_atlas <- merge(
  x = filtered_seurat_list[[1]], 
  y = filtered_seurat_list[-1], 
  add.cell.ids = names(filtered_seurat_list),
  project = "Organ_Specific_CD_Atlas"
)

normed_atlas <- NormalizeData(combined_atlas)


normed_atlas <- FindVariableFeatures(normed_atlas, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(normed_atlas), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(normed_atlas)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
plot2


normed_atlas <- ScaleData(normed_atlas)

# 4. Dimensionality Reduction (PCA)
normed_atlas <- RunPCA(normed_atlas)


VizDimLoadings(normed_atlas, dims = 1:2, reduction = "pca")


DimHeatmap(normed_atlas, dims = 1, cells = 500, balanced = TRUE)
#thinking I want to go back and keep all compartments separate. 
#TODO: redo but separate compartments
#TODO: generate quarto report with parameters for each compartment 


processed_list <- map(filtered_seurat_list, function(obj) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  return(obj)
})
library(ggplot2)
plot_list <- map(processed_list, function(obj) {
  
  sample_name <- obj$Tissue[[1]] 

  # Identify the top 10 most variable genes in this specific compartment
  top10 <- head(VariableFeatures(obj), 10)
  
  # Plot variable features with and without labels
  plot1 <- VariableFeaturePlot(obj) + 
           ggtitle(paste("Variable Features:", sample_name))
  
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  
  return(plot2)
})


# Apply Scaling, PCA, Clustering, and UMAP to each compartment independently
fully_analyzed_list <- map(processed_list, function(obj) {
  
  # --- PHASE 3: SCALING AND PCA ---
  # The ScaleData() function shifts the expression of each gene, so that the mean expression across cells is 0[cite: 198].
  # It also scales the expression of each gene, so that the variance across cells is 1[cite: 198].
  obj <- ScaleData(obj)
  
  # We perform linear dimensional reduction by running RunPCA()[cite: 206]. 
  # By default, only the previously determined variable features are used as input[cite: 206].
  obj <- RunPCA(obj)
  
  # --- PHASE 4: CLUSTERING & VISUALIZATION ---
  # We embed cells in a graph structure, building a K-nearest neighbor (KNN) graph[cite: 233, 234].
  # This step is performed using the FindNeighbors() function, and takes as input the previously defined dimensionality of the dataset[cite: 235]. 
  # (Note: We use 1:10 as a standard baseline, but ideally you would check an ElbowPlot for each object first)
  obj <- FindNeighbors(obj, dims = 1:10)
  
  # The FindClusters() function implements this procedure, and contains a resolution parameter that sets the 'granularity' of the downstream clustering[cite: 237].
  # Increased values lead to a greater number of clusters[cite: 237].
  obj <- FindClusters(obj, resolution = 0.5)
  
  # We run non-linear dimensional reduction techniques like UMAP to visualize and explore these datasets[cite: 240, 245].
  obj <- RunUMAP(obj, dims = 1:10)
  
  return(obj)
})


DimPlot(fully_analyzed_list[[1]], reduction = "umap", label = TRUE) + 
  ggtitle("Colon Stroma Compartment")

library(Seurat)
library(tidyverse)

# 1. Split the Colon Stromal object into a list of objects based on patient identity
# (Assuming your active identities are set to the patient IDs, or use split.by = "orig.ident")
stromal_list <- SplitObject(colon_stromal, split.by = "ident")

# 2. Define custom thresholds ONLY for the tricky patients.
# You don't need to type out all 30 patients if most of them are normal!
patient_thresholds <- list(
  "I124246" = list(min_feat = 500, max_feat = 6000, max_count = 25000, max_mt = 8),
  "N104689" = list(min_feat = 500, max_feat = 5500, max_count = 20000, max_mt = 10),
  "H106265" = list(min_feat = 200, max_feat = 2000, max_count = 8000,  max_mt = 5)
)

# 3. Define a safe "default" baseline for all the other patients who behave normally
default_thresh <- list(min_feat = 200, max_feat = 3000, max_count = 12000, max_mt = 5)

# 4. Apply the tailored filters using imap
filtered_stromal_list <- imap(stromal_list, function(obj, patient_id) {
  
  # Fetch the custom threshold for this patient, OR use the default if not listed
  thresh <- patient_thresholds[[patient_id]]
  if (is.null(thresh)) {
    thresh <- default_thresh
  }
  
  # Calculate MT% if you haven't already
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # Apply the strict, tailored filtering
  obj <- subset(obj, 
    subset = nFeature_RNA > thresh$min_feat & 
             nFeature_RNA < thresh$max_feat & 
             nCount_RNA < thresh$max_count &
             percent.mt < thresh$max_mt
  )
  
  return(obj)
})

# 5. Merge the pristine patient data back into a single Colon Stromal object
colon_stromal_clean <- merge(
  x = filtered_stromal_list[[1]], 
  y = filtered_stromal_list[-1]
)
