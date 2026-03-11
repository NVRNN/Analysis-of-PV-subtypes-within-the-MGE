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
project_name <- "athena"
user <- Sys.info()["login"]

base_dir <- paste("/media", user, "TOSHIBA EXT/bioinformatics", sep = "/")
data_dir <- paste(base_dir, "data", sep = "/")
out_dir <- paste(base_dir, "analysis_output", sep = "/")
dir.create(out_dir, showWarnings = FALSE)

setwd(out_dir)


# Setup -----------------------------------------------------------------------------------------------------------

marker_genes <- c("FGDBP3", "ALPL", "BTBD7", "WNT7B", "GLI1", "CCND1")



# Load counts data ------------------------------------------------------------------------------------------------

# Metadata
meta <- read_tsv(file = paste(data_dir, "metadata.tsv", sep = "/"))
meta

# Check composition
meta %>%
  select(tissue, age) %>%
  table()


# Counts
counts <- read_tsv(paste(data_dir, "GSE109796_Oscar.GEO.singleCell.gene.count.txt.gz", sep = "/"))
counts
dim(counts)
# [1] 37310  2670


# Make the column names of the counts data match the metadata
meta
colnames(counts) %>% head()

counts %>%
  rename_with(~ str_split_fixed(string = ., pattern = "_", n = 2)[,1])
# Error about duplicate names

# Are there duplicates in the metadata?
meta %>%
  group_by(sample) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
# All good there, check the counts

counts %>%
  colnames() %>%
  gsub(pattern = "-", replacement = "_", x = .) %>%
  str_split_fixed(string = ., pattern = "_", n = 5)[,1:3]
  
counts <- 
  counts %>%
  # select(-geneList) %>%
  rename_with(~ sub(pattern = "-2-", replacement = "_", x = .)) %>%
  rename_with(~ gsub(pattern = "-", replacement = "_", x = .)) %>%
  rename_with(~ ifelse(. == "geneList", "geneList", str_match(string = ., pattern = "^([A-Z]\\d_\\d+_[A-Z]\\d+)_")[,2]))
counts


# Convert hyphens to underscores and remove -2- from meta
meta <- 
  meta %>%
  mutate(sample = sub(pattern = "-2-", replacement = "_", x = sample),
         sample = gsub(pattern = "-", replacement = "_", x = sample))
meta


# Check that all counts samples have corresponding metadata
summary(meta$sample %in% colnames(counts))
# Good



# Ensure gene names are unique ------------------------------------------------------------------------------------

test <- CreateSeuratObject(counts = 
                             counts %>% 
                             mutate(geneList = str_split_fixed(string = geneList, pattern = "\\|", n = 2)[,2]) %>%
                             column_to_rownames(var = "geneList"),
                           meta.data = meta %>% column_to_rownames(var = "sample"))


counts %>% pull(geneList) %>% length()
counts %>% pull(geneList) %>% unique() %>% length()
# Gene names currently unique due to Ensembl prefix


# Look into non-unique genes
counts %>% 
  mutate(geneList = str_split_fixed(string = geneList, pattern = "\\|", n = 2)[,2]) %>%
  group_by(geneList) %>%
  mutate(n = n()) %>%
  select(geneList, n, everything()) %>%
  ungroup() %>%
  arrange(desc(n), geneList) %>%
  filter(n > 1) %>% pull(geneList) %>% unique()
# 357 genes
# Need to figure out how to deal with counts values


# Which is the most highly expressed?
counts %>% 
  mutate(geneList = str_split_fixed(string = geneList, pattern = "\\|", n = 2)[,2]) %>%
  group_by(geneList) %>%
  mutate(n = n()) %>%
  select(geneList, n, everything()) %>%
  ungroup() %>%
  filter(n > 1) %>%
  select(-n) %>%
  pivot_longer(-geneList, names_to = "cell", values_to = "count") %>%
  group_by(geneList) %>%
  summarise(total_counts = sum(count)) %>%
  arrange(desc(total_counts)) %>% View()
# QK is the most abundantly counted gene (at least in raw counts)

# Look at this gene
counts %>%
  filter(grepl(pattern = "\\|QK", x = geneList)) %>%
  pivot_longer(-geneList) %>%
  group_by(geneList) %>%
  summarise(sum = sum(value))
# Counts are not equally split, so perhaps reasonable just to add the values from the different isoforms

counts %>%
  filter(grepl(pattern = "\\|U6$", x = geneList)) %>%
  pivot_longer(-geneList) %>%
  group_by(geneList) %>%
  summarise(sum = sum(value)) %>%
  arrange(desc(sum)) %>%
  ggplot(aes(x = geneList, y = sum)) +
  geom_point() +
  theme(axis.text.x = element_blank())
# Let's just add counts for gene isoforms together


# Operate only on duplicated genes
collapsed_genes_long <- counts %>%
  mutate(geneList = str_split_fixed(string = geneList, pattern = "\\|", n = 2)[,2]) %>%
  group_by(geneList) %>%
  mutate(n = n()) %>%
  select(geneList, n, everything()) %>%
  ungroup() %>%
  filter(n > 1) %>%
  select(-n) %>%
  pivot_longer(-geneList, names_to = "cell", values_to = "count") %>%
  group_by(cell, geneList) %>%
  summarise(total_counts = sum(count)) #%>%
  # pivot_wider(names_from = cell, values_from = total_counts)
collapsed_genes_long

# Remove duplicate genes from the main counts data
unique_gene_counts_long <- 
  counts %>%
  mutate(geneList = str_split_fixed(string = geneList, pattern = "\\|", n = 2)[,2]) %>%
  group_by(geneList) %>%
  mutate(n = n()) %>%
  select(geneList, n, everything()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n) %>%
  pivot_longer(-geneList, names_to = "cell", values_to = "total_counts") %>%
  select(cell, geneList, total_counts)
unique_gene_counts_long

# Combine the two gene sets
clean_counts <- 
  rbind(unique_gene_counts_long, collapsed_genes_long) %>%
  pivot_wider(names_from = "cell", values_from = "total_counts")
clean_counts

# cleaning out not needed packages
rm(collapsed_genes_long, counts, unique_gene_counts_long)
#garbage collection to remove old stuff
gc()

# Quick check - compare to previous look at this gene
clean_counts %>%
  select(geneList, "C1_101_A1") %>%
  filter(geneList == "2410002O22RIK")
# 861


# See which cells have <50k counts (a filter applied in the paper)
clean_counts %>%
  pivot_longer(-geneList, names_to = "barcode", values_to = "count") %>%
  group_by(barcode) %>%
  summarise(n = sum(count)) %>%
  mutate(low_counts = n < 50000) %>%
  filter(low_counts)
# No cells with < 50k counts; must have been removed in the data already
# Can't do the next part of QC filtering; from the paper:
# "exonic read distribution, read distribution across different chromosomes, GC content distribution and gene expression distribution"
# - removed 666 cells



# Create a Seurat object ------------------------------------------------------------------------------------------

seurat <- CreateSeuratObject(counts = clean_counts %>% column_to_rownames(var = "geneList"),
                             meta.data = meta %>% column_to_rownames(var = "sample"),
                             min.cells = 10)
dim(seurat)
# [1] 18788  2669

# Check out the metadata
seurat@meta.data %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "barcode")


# Add sample name to metadata
seurat@meta.data <- 
  seurat@meta.data %>%
  mutate(sample = paste(tissue, age, sep = "_"))

seurat@meta.data %>% head()


# Check that cells numbers are correct
seurat@meta.data %>%
  select(sample) %>%
  table()
# TODOx: Check these match the paper: 1000% confirmed accurate
# Yeah pretty close

# cleaning out not needed packages
rm(clean_counts, collapsed_genes_long, counts, meta, unique_gene_counts_long)
#garbage collection to remove old stuff
gc()

# Doublet detection -----------------------------------------------------------------------------------------------

# Plan will be:
# - identify and remove doublets
# - perform test clustering and investigate other QC metrics
# - set criteria for filtering; filter

# Calculate mitochondrial content
grep(pattern = "^MT-", x = rownames(seurat), value = TRUE)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

# Quick looks at QC metrics
VlnPlot(object = seurat, features = "nCount_RNA", group.by = "sample")
ggsave(filename = "01-QC-nCountRNA.png")

VlnPlot(object = seurat, features = "nFeature_RNA", group.by = "sample")
ggsave(filename = "01-QC-nFeatureRNA.png")

VlnPlot(object = seurat, features = "percent.mt", group.by = "sample") + geom_hline(yintercept = 10)
ggsave(filename = "01-QC-percent_mt.png")


# Doublet detection
dbl <- scDblFinder(sce = GetAssayData(seurat), samples = seurat@meta.data$sample)
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
seurat@meta.data <-
  seurat@meta.data %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "barcode") %>%
  left_join(y = dbl@colData %>% as_tibble(rownames = NA) %>% rownames_to_column(var = "barcode") %>% select(barcode, scDblFinder.class),
            by = "barcode") %>%
  column_to_rownames(var = "barcode")
dim(seurat)
# [1] 18788  2669

seurat@meta.data %>% head()

# Remove doublets
singlets <- subset(seurat, subset = scDblFinder.class == "singlet")
dim(singlets)
# [1] 18788  2558

rm(seurat, clean_counts, collapsed_genes_long, counts, dbl, meta, unique_gene_counts_long)
gc()


# Diagnostic clustering -------------------------------------------------------------------------------------------

# See how QC metrics look in each cluster to help guide filtering thresholds
# Don't need to keep this after filtering
# Get the process from the Seurat tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

scaled_singlets <-
  singlets %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(singlets))

scaled_singlets <-
  scaled_singlets %>%
  RunPCA(features = VariableFeatures(object = scaled_singlets)) %>%
  FindNeighbors(dims = 1:30) %>%
  RunTSNE(dims = 1:30) %>%
  RunUMAP(dims = 1:30)

scaled_singlets %>% FindClusters(resolution = 0.8) %>% DimPlot() # 9 clusters
scaled_singlets %>% FindClusters(resolution = 0.4) %>% DimPlot() # 6 clusters
scaled_singlets %>% FindClusters(resolution = 1.5) %>% DimPlot() # 11 clusters
scaled_singlets %>% FindClusters(resolution = 2) %>% DimPlot() # 13 clusters
# Let's choose this resolution

scaled_singlets <- 
  scaled_singlets %>% 
  FindClusters(resolution = 2)

scaled_singlets %>% 
  DimPlot(reduction = "tsne", split.by = "age") 
# Doesn't really look like the paper's images, but we haven't done QC yet, and I think they processed the ages separately

# Look at QC metrics by cluster
scaled_singlets %>% VlnPlot(features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "seurat_clusters")
# Filtering percent.mt at 10% wouldn't be removing any entire cell-types, so we'll choose that threshold
# Cells were already required to have > 50k transcripts
# No cells are expressing way more genes than the other, and same for transcripts, so no other filters required
# Can also plot these metrics on the UMAP

scaled_singlets %>%
  FeaturePlot(features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))



# QC filtering ----------------------------------------------------------------------------------------------------

filtered <-
  scaled_singlets %>% subset(subset = percent.mt < 10)

dim(filtered)

# Save filtered data
SaveH5Seurat(object = filtered, filename = "filtered_data.h5Seurat")
rm(singlets, scaled_singlets)
gc()


# Cell cycle regression -----------------------------------------------------------------------------------------------------
# https://satijalab.org/seurat/articles/cell_cycle_vignette
# Load scRNA data if required
filtered <- LoadH5Seurat("filtered_data.h5Seurat")

# TODO: Do we want to integrate? All samples? Same tissue only?
# Integrate all samples
# - Interested in gene expression programs during differentiation
# - Paper keeps different tissues together
# perform cell cycle regression before Intergration https://github.com/satijalab/seurat/issues/5879

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

filtered <- NormalizeData(filtered)
filtered <- FindVariableFeatures(filtered, selection.method = "vst")
filtered <- ScaleData(filtered, features = rownames(filtered))

filtered <- RunPCA(filtered, features = VariableFeatures(filtered), ndims.print = 6:10, nfeatures.print = 10)

DimHeatmap(filtered, dims = c(8, 10))

filtered <- CellCycleScoring(filtered, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(filtered[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(filtered, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
filtered <- RunPCA(filtered, features = c(s.genes, g2m.genes))
DimPlot(filtered)

filtered <- ScaleData(filtered, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(filtered))

# Now, a PCA on the variable genes no longer returns components associated with cell cycle
filtered <- RunPCA(filtered, features = VariableFeatures(filtered), nfeatures.print = 10)

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
filtered <- RunPCA(filtered, features = c(s.genes, g2m.genes))
DimPlot(filtered, reduction = "pca")

# cell cycle differences removed yay

# Save filtered data after regressing cell cycle
SaveH5Seurat(object = filtered, 
             filename = "filtered.cell_cycle.h5Seurat",
             overwrite = TRUE)



# Integration -----------------------------------------------------------------------------------------------------

# Change this to TRUE if you want to reload the data at this point
reload <- false

if(reload){
  filtered <- LoadH5Seurat(file = "filtered.cell_cycle.h5Seurat")
}

library(patchwork)
# way to look at meta data of object
# filtered@meta.data %>% head() 


#separate into identifying lable (Sample in this case)
filtered.list <- SplitObject(filtered, split.by = "sample")

# normalize and identify variable features for each dataset independently 
filtered.list <- lapply(X = filtered.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = filtered.list)

anchors<- FindIntegrationAnchors(object.list = filtered.list, anchor.features = features)

# this command creates an 'integrated' data assay
# TODO: Change anchors.combined to "integrated" to make things clearer
integrated <- IntegrateData(anchorset = anchors)
integrated

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
# TODO: Copy the line above but change RunUMAP to RunTSNE
integrated <- RunTSNE(integrated, reduction = "pca", dims = 1:30)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.8) # Found 6 when 0.5 clusters7 when 0.8

# Visualization
p1 <- DimPlot(integrated, reduction = "umap", group.by = "sample")
p2 <- DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(integrated, reduction = "umap", split.by = "sample")

#need to install following packages 
# install.packages('BiocManager')
# BiocManager::install('multtest')
# install.packages('metap')
library(BiocManager)
library(multtest)
library(metap)

# install.packages("qqconf")

# Identify conserved cell type markers change cluster investigated 0-5 to investigate
# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(integrated) <- "RNA"
cluster_0_markers <- FindConservedMarkers(integrated, ident.1 = 0, grouping.var = "sample", verbose = FALSE)
head(cluster_0_markers, n = 9)

# Make UMAP and input features to investigate identified in last step
FeaturePlot(integrated,
            features = rownames(cluster_0_markers)[1:9],
            min.cutoff = "q9",
            reduction = "umap")

# Make UMAP and visualise genes from paper Progenitors 
FeaturePlot(integrated,
            features = c("GFAP", "SLC1A3", "HES1", "HES5", "LHX2", "NOTCH1", "NOTCH2", "NOTCH3", "VIM", "SOX2", "SAL3", "SPARC", "PDPN"), min.cutoff = "q9",
            reduction = "umap")
# Neurons 
FeaturePlot(integrated,
            features = c("TUBB3", "MAPT", "DCX", "NRXN3", "STMN2", "STMN3", "GNG3","SCN3A", "MAP2", "GAP43", "TMEM130", "L1CAM", "GRIA2"), min.cutoff = "q9",
            reduction = "umap")
# VZ
FeaturePlot(integrated,
            features = c("FABP7", "SLC1A3", "1190002H23RIK", "SNAP23", "FGFBP3", "RHPN1"), min.cutoff = "q9",
            reduction = "umap")
# SVZ
FeaturePlot(integrated,
            features = c("RSRC2", "NRXN3", "ST18", "CD24A", "DLX6OS1"), min.cutoff = "q9",
            reduction = "umap")
#dMGE
FeaturePlot(integrated,
            features = c("NRXN1", "JAG1", "NKX6-2","GM5069", "NPY", "NEK7", "8430410K20RIK", "MKS1", "NRP1","EBF1","HSPA12A","WNK3-PS""2310014H01RIK","MAVS","LGALS1","ZFP568","PTDSS2","GOT2","C130036L24RIK"), min.cutoff = "q9",
            reduction = "umap")
#vMGE
FeaturePlot(integrated,
            features = c("Lhx8", "Nkx2-1", "Zic1","Mbip", "Asb4", "Sulf2", "Sez6", "Crabp2", "Lhx6","Dach1"), min.cutoff = "q9",
            reduction = "umap")
#PV1 subtype markers 
#for nestin data set FeaturePlot(integrated,
#            features = c("Lima1", "Ccnd2", "Ywhaz","Clip2", "Zfp191", "Gm5805", "U4atac", "Ndufa13", "Scaper","Lhx6", "Fam110a", "Clasrp", "S100a16", "1810009a15Rik", "Fam100a", "Tomm20", "Gm17306", "Pdcd11", "Slc29a4"), "Usf1", "Fkbp3", min.cutoff = "q9",
 #           reduction = "umap")
FeaturePlot(integrated,
            features = c("LIMA1", "CCND2", "YWHAZ","CLIP2", "ZFP191", "GM5805", "U4ATAC", "NDUFA13", "SCAPER","LHX6", "FAM110A", "CLASRP", "S100A16", "FAM100A", "TOMM20", "GM17306", "PDCD11", "SLC29A4", "USF1", "FKBP3"), min.cutoff = "q9",
            reduction = "umap")
#PV2 subtype markers
FeaturePlot(integrated,
            features = c("TTI1", "SEPT4", "PTPRG","FAM190A", "ZFP398", "TMEM5", "4921531C22RIK", "NR2C1", "THOC5","LCMT1", "RSPH3A", "SPATA13", "DYRK1B", "ELP3", "KIF26A", "38330406C13RIK", "EXOC6", "ATP1A2", "RABL3", "D11WSU47E", "POFUT1"), min.cutoff = "q9",
            reduction = "umap")
#PV3 subtype markers 
FeaturePlot(integrated,
            features = c("TRPS1", "CXCR7", "4930412F15RIK","PI4KB", "AMPH", "ACP6", "ERBB4", "LETM1", "CEP135","CPLX1", "PRKD2", "SYNJ2BP", "RNF139", "GUCY1A3", "RBBP5", "DNMBP", "RNF38", "NMNAT2", "1810014B01RIK", "4933421E11RIK"), min.cutoff = "q9",
            reduction = "umap")
#PV4 subtype markers 
FeaturePlot(integrated,
            features = c("ITPR1", "TIAM2", "7SK","GM17659", "B430010I23RIK", "DECR2", "ZFYVE9", "HIST1H4D" "FOXN3", "ABT1", "PAM", "LYST", "ABCD2", "GPBP1L1", "LPMK", "ATG2B", "DOM3Z", "TEF", "SECISBP2L", "MDN1"), min.cutoff = "q9",
            reduction = "umap")

FeaturePlot(integrated,
            features = c("SERPINE2", "TH", "CDKN1C","PDE1C", "ETV1", "SPP1", "LIPG", "EPHB3", "ZIC4","2610035D17RIK"), min.cutoff = "q9",
            reduction = "umap")
integrated <- RenameIdents(integrated, `0` = "CCND2", `1` = "ALPL", `2` = "BTBD7",
                                `3` = "WNT7B", `4` = "GLI1", `5` = "CCND1", `6` = "YWHAZ")
DimPlot(integrated, label = TRUE)

TSNEPlot(integrated, features = c("CCND2", "CCND1", "SOX2", "GLI1", "SST", "YWHAZ", "YWHAE", "SOX10", "GSX2", "SHH", "TTF-1", "GPSM2", "Fabp7", "Lhx6", "Gad1", "Dcx", "Etv1", "Sst", "Olig2"), min.cutoff = "q9")
FeaturePlot(integrated,
            features = c("CCND2", "CCND1", "SOX2", "GLI1", "SST", "YWHAZ", "YWHAE", "SOX10", "GSX2", "SHH", "TTF1", "GPSM2", "FABP7", "LHX6", "GAD1", "DCX", "ETV1", "SST", "OLIG2"), min.cutoff = "q9",
            reduction = "tsne")


# Try a different approach to finding markers; this will find cluster-specific markers
cluster_markers <- FindAllMarkers(object = integrated, min.pct = 0.25, logfc.threshold = 0.25)

# Look at the output
cluster_markers %>% head()

# Make a dot plot to show markers genes across clusters
top5_markers <-
  cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

top5_markers

DotPlot(object = integrated, features = unique(top5_markers$gene)) +
  coord_flip()



# To explore different clustering resolutions -----------------------------

DefaultAssay(integrated) <- "integrated"

# Play around with clustering resolution without saving the output
FindClusters(integrated, resolution = 1)

integrated <- FindClusters(integrated, resolution = 0.5) # Found 6 clusters

# Visualization
p1 <- DimPlot(integrated, reduction = "umap", group.by = "sample")
p2 <- DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(integrated, reduction = "umap", split.by = "sample")
DefaultAssay(integrated) <- "RNA"

# Find markers for each cluster
cluster_markers <- FindAllMarkers(object = integrated, min.pct = 0.25, logfc.threshold = 0.25)

# Look at the output
cluster_markers %>% head()

# Make a dot plot to show markers genes across clusters
top5_markers <-
  cluster_markers %>% 
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

top5_markers

DotPlot(object = integrated, features = unique(top5_markers$gene)) +
  coord_flip()



#########################################################
integrated@meta.data %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "barcode") %>%
  group_by(sample, seurat_clusters) %>%
  summarise(n = n())


VlnPlot(integrated, features = c("CCND1", "TTF1"))



#loading in the filtered cell cycle data 
filtered <- LoadH5Seurat("~/Downloads/filtered.cell_cycle.h5Seurat")
