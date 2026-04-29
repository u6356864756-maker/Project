install.packages("Seurat")
install.packages("BiocManager")
BiocManager::install("BUSpaRse") # optional, for utilities


library(europepmc)
library(ggupset)
library(wordcloud)
library(ggplot2)
library(Seurat)
library(patchwork)
library(Seurat)
library(dplyr)
library(Matrix)

library(Seurat)
library(Matrix)

# ── Set your paths ──────────────────────────────────────────────
raw_dir  <- "C:/Users/hunte/OneDrive/Desktop/R/Project/SCP2019/expression/raw"
norm_dir <- "C:/Users/hunte/OneDrive/Desktop/R/Project/SCP2019/expression/normalized"

# ── Load RAW counts ─────────────────────────────────────────────
raw_counts <- ReadMtx(
  mtx      = file.path(raw_dir, "Carotid_Expression_Matrix_raw_counts_V1.mtx"),
  cells    = file.path(raw_dir, "Carotid_Expression_Matrix_barcodes_V1.tsv"),
  features = file.path(raw_dir, "Carotid_Expression_Matrix_genes_V1.tsv")
)

# ── Load NORMALIZED counts ───────────────────────────────────────
norm_counts <- ReadMtx(
  mtx      = file.path(norm_dir, "Carotid_Expression_Matrix_norm_V1.mtx"),
  cells    = file.path(norm_dir, "Carotid_Expression_Matrix_barcodes_norm_V1.tsv"),
  features = file.path(norm_dir, "Carotid_Expression_Matrix_genes_norm_V1.tsv")
)

# ── Create Seurat object from RAW counts (correct approach) ──────
df <- CreateSeuratObject(
  counts       = raw_counts,
  project      = "SCP2019_VSMC_Carotid",
  min.cells    = 3,
  min.features = 200
)

# ── Store normalized matrix as a separate assay ──────────────────
# First align cells — normalized may have slightly different barcodes
shared_cells <- intersect(colnames(df), colnames(norm_counts))
df   <- df[, shared_cells]

df[["normalized"]] <- CreateAssayObject(
  data = norm_counts[, shared_cells]  # 'data' slot = normalized
)

# ── Verify ───────────────────────────────────────────────────────
df          # should show ~6049 cells, ~33694 genes
Assays(df)  # should show "RNA" and "normalized"


raw_genes  <- rownames(df[["RNA"]])
norm_genes <- rownames(df[["normalized"]])

length(raw_genes)               # 22,993
length(norm_genes)              # should be ~33,694
length(intersect(raw_genes, norm_genes))  # how many overlap
length(setdiff(norm_genes, raw_genes))    # genes only in normalized
length(setdiff(raw_genes, norm_genes))    # genes only in raw


##### PipeLine Start 
DefaultAssay(df) <- "RNA"


#Making PErcent MT col 
df[['percent.mt']]<-PercentageFeatureSet(df, pattern = "^MT-")
VlnPlot(df, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "percent.mt")
#The amount of mitochondrial DNA is essentially 0. Because of this no samples will be removed. 

#normalize data 
#Normalizing the DATA 
df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 1e4)

#Variable Features
df<- FindVariableFeatures(df, selection.method = "vst", nfeatures = 2000)
top10<- head(VariableFeatures(df))


#Step 4 Scaling (standarization) method 1 
all.genes <- rownames(df)
#list of all genes present 
df<-ScaleData(df, features = all.genes)


df <-RunPCA(df, features = VariableFeatures(df))
#print x # of PCA vectors, and names of # of genes/PCA vector 
print(df[["pca"]], dims = 1:5, nfeatures = 6)

ElbowPlot(df, ndims = 30)
#PC 1–4: steep drop, clearly informative
#PC 5–10: moderate drop, still capturing real biology
#PC 10–15: curve starts flattening noticeably
#PC 15+: mostly noise, marginal gains

#Use PCA 15 

#Cluster
df <- FindNeighbors(df, dims = 1:15)
df <- FindClusters(df, resolution = 0.5)
Idents(df)

df <-RunUMAP(df, dims = 1:15)

DimPlot(df, reduction = "umap", label = TRUE)
DimPlot(df, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)


markers <- FindAllMarkers(df, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)
library(tidyverse)
#markers is going to contain distinct markers for each group indicidually 
markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)


cluster1.markers <- FindMarkers(df, ident.1 = 0, min.pct = 0.25)
dim(cluster1.markers) # 1.3 k genes
head(cluster1.markers, n = 5)
