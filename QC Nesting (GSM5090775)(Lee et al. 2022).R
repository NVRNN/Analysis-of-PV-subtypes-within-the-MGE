library(tidyverse)
library(Seurat)
library(SeuratDisk)


#BiocManager::install("scDblFinder")
#install.packages("scDblFinder")
library(scDblFinder)
# install.packages("tidyverse")
#install.packages("Seurat")

# fast map had to be installed as linux didnt come with it https://blog.zenggyu.com/posts/en/2018-01-29-installing-r-r-packages-e-g-tidyverse-and-rstudio-on-ubuntu-linux/index.html
# instal this in linux terminal https://blog.zenggyu.com/posts/en/2018-01-29-installing-r-r-packages-e-g-tidyverse-and-rstudio-on-ubuntu-linux/index.html
#installing *source* package ‘fastmap’
#install.packages("fastmap")

# install.packages("tidyverse")

# for instilation with current version of R
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("SeuratDisk")

#if (!requireNamespace("remotes", quietly = TRUE)) {
#  install.packages("remotes")
#}
#remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
# project_name <- "230622NathanielNapoli"
project_name <- "Nestin"
user <- Sys.info()["login"]

base_dir <- paste("/media", user, "TOSHIBA EXT/Bioinformatics_Nestin", sep = "/")
data_dir <- paste(base_dir, "data", sep = "/")
out_dir <- paste(base_dir, "analysis_output", sep = "/")
dir.create(out_dir, showWarnings = FALSE)

setwd(out_dir)

# QC 
library(dplyr)
library(patchwork)

#install.packages("Matrix", type = "source")
#install.packages("irlba", type = "source")

Nestin <-LoadH5Seurat("Seurat_Nestin.h5seurat")

# Doublet detection -----------------------------------------------------------------------------------------------

# Plan will be:
# - identify and remove doublets
# - perform test clustering and investigate other QC metrics
# - set criteria for filtering; filter

# Calculate mitochondrial content
grep(pattern = "^mt-", x = rownames(Nestin), value = TRUE)
Nestin[["percent.mt"]] <- PercentageFeatureSet(Nestin, pattern = "^mt-")

# Quick looks at QC metrics
VlnPlot(object = Nestin, features = "nCount_RNA")
ggsave(filename = "01-QC-nCountRNA.png")

VlnPlot(object = Nestin, features = "nFeature_RNA")
ggsave(filename = "01-QC-nFeatureRNA.png")

VlnPlot(object = Nestin, features = "percent.mt")
ggsave(filename = "01-QC-percent_mt.png")

#QC based on features and mitochondrial expression
Nestin <- subset(Nestin, subset = nFeature_RNA > 1500 & nFeature_RNA < 4500 & percent.mt < 10)

#Normalize the data
Nestin <- NormalizeData(Nestin, normalization.method = "LogNormalize", scale.factor = 10000)

Nestin <- NormalizeData(Nestin)

#Identification of highly variable features (feature selection)
Nestin <- FindVariableFeatures(Nestin, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Nestin), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Nestin)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
ggsave(filename = "01-QC-variable genes.png")

#Scaling the data
all.genes <- rownames(Nestin)
Nestin <- ScaleData(Nestin, features = all.genes)

#Perform linear dimensional reduction PCA
Nestin <- RunPCA(Nestin, features = VariableFeatures(object = Nestin))

# Examine and visualize PCA results a few different ways
print(Nestin[["pca"]], dims = 1:5, nfeatures = 5)
#graphing
VizDimLoadings(Nestin, dims = 1:2, reduction = "pca")

DimPlot(Nestin, reduction = "pca")
ggsave(filename = "PC1.png")

DimHeatmap(Nestin, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(Nestin, dims = 1:15, cells = 500, balanced = TRUE)
ggsave(filename = "PC1-15.png")

#Dimensionality of the data
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
#Nestin <- JackStraw(Nestin, num.replicate = 100)
#Nestin <- ScoreJackStraw(Nestin, dims = 1:20)

#JackStrawPlot(Nestin, dims = 1:15)

ElbowPlot(Nestin, ndims = 50)

#Cluster
Nestin <- FindNeighbors(Nestin, dims = 1:25)
Nestin <- FindClusters(Nestin, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(Nestin), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')

Nestin <- RunUMAP(Nestin, dims = 1:25)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Nestin, reduction = "umap")
ggsave(filename = "UMAP_Cluster.png")

#save umap
saveRDS(Nestin, file = "Nestin umap")

cluster1.markers <- FindMarkers(Nestin, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 10)

cluster5.markers <- FindMarkers(Nestin, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

Nestin.markers <- FindAllMarkers(Nestin, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Nestin.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

 

Nestin.markers <- FindMarkers(Nestin, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(Nestin, features = c("Shh", "ttf1"))
ggsave(filename = "Shh.png")

FeaturePlot(Nestin, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))

# Make UMAP and visualise genes from paper Progenitors 
FeaturePlot(Nestin,
            features = c("Gfap", "Slc1a3", "Hes1", "Hes5", "Lhx2", "Notch1", "Notch2", "Notch3", "Vim", "Sox2", "Sall3", "Sparc", "Pdpn"), min.cutoff = "q9",
            reduction = "umap")
ggsave("progenitor.png")

# Neurons 
FeaturePlot(Nestin,
            features = c("Tubb3", "Mapt", "Dcx", "Nrxn3", "Stmn2", "Stmn3", "Gng3","Scn3a", "Map2", "Gap43", "Tmem130", "L1cam", "Gria2"), min.cutoff = "q9",
            reduction = "umap")
ggsave("Neurons.png")

# VZ
FeaturePlot(Nestin,
            features = c("Fabp7", "Slc1a3", "1190002H23rik", "Snap23", "Fgfbp3", "Rhpn1"), min.cutoff = "q9",
            reduction = "umap")
ggsave(filename = "VZ.png")

# SVZ
FeaturePlot(Nestin,
            features = c("Tubb3", "Lhx6", "St18", "Nrxn3", "Dlx6os1", "Nudt4", "Gad2", "Dcx", "Runx1t1"), min.cutoff = "q9",
            reduction = "umap")
ggsave("SVZ.png")

#dMGE E14.5
FeaturePlot(Nestin,
            features = c("Gm10837", "Gm10717", "Gm6984", "Gm10282", "Atp5g1", "Gm10801", "Eid2", "Uqcr11", "Serf2", "Sfi1", "Gm8730", "Gm9846", "Gm8759", "Atp5g2", "Gm10106", "Gria3", "Cwc22", "Sepw1", "Gm10718" ), min.cutoff = "q9",
            reduction = "umap")
#dmge E12.5
FeaturePlot(Nestin,
            features = c("Nrxn1", "Jag1", "Nkx6-2","Gm5069", "Npy", "Nek7", "8430410K20rik", "Mks1", "Nrp1","Ebf1","Hspa12a","Wnk3-PS""2310014H01rik","Mavs","Lgals1","Zfp568","Ptdss2","Got2","C130036l24rik" ), min.cutoff = "q9",
            reduction = "umap")
ggsave("dMGE E12.png")

# assorted markers
#vMGE
FeaturePlot(Nestin,
            features = c("Shh", "Ywhab", "Lhx8", "Nkx2-1", "Zic1","Mbip", "Asb4", "Sulf2", "Sez6", "Crabp2", "Lhx6","Dach1", "Serpine2", "Th", "Cdkn1c", "Pde1c", "Etv1", "Spp1", "Lipg", "Ephb3", "Zic4", ), min.cutoff = "q9",
            reduction = "umap")
ggsave(filename = "vMGE.png")

FeaturePlot(Nestin,
            features = c("Shh", "Ywhab", "Lhx8", "Nkx2-1", "Zic1","Mbip", "Asb4", "Sulf2", "Sez6", "Crabp2", "Lhx6","Dach1"), min.cutoff = "q9",
            reduction = "umap")
ggsave(filename = "vMGE.png")

#assorted markers
FeaturePlot(Nestin,
            features = c("Nkx6-2", "Ywhab", "Ywhag", "Ywhah", "Ccnd1", "Ccnd2", "Fabp7", "Gpsm2", "Dcx", "Sox2"), min.cutoff = "q9",
            reduction = "umap")
ggsave(filename = "progenitor.png")

Nestin.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
DoHeatmap(Nestin, features = top10$gene) + NoLegend()












# Doublet detection
dbl <- scDblFinder(sce = GetAssayData(Nestin), samples = seurat@meta.data$sample)
# 168 doublets detected (6.3%)

# Summarise
dbl@colData %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "barcode") %>%
  select(barcode, scDblFinder.sample, scDblFinder.class) %>%
  group_by(scDblFinder.sample, scDblFinder.class) %>%
  summarise(n = n()) %>%
  mutate(total = sum(n),
         proportion = n / total) %>%
  filter(scDblFinder.class == "doublet")

# Add doublet info to seurat object
Nestin@meta.data <-
  Nestin@meta.data %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "barcode") %>%
  left_join(y = dbl@colData %>% as_tibble(rownames = NA) %>% rownames_to_column(var = "barcode") %>% select(barcode, scDblFinder.class),
            by = "barcode") %>%
  column_to_rownames(var = "barcode")
dim(Nestin)
# [1] 18788  2669

Nestin@meta.data %>% head()

# Remove doublets
singlets <- subset(Nestin, subset = scDblFinder.class == "singlet")
dim(singlets)
# [1] 18788  2558