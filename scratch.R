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

