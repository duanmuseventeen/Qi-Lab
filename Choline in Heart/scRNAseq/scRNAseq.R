rm(list = ls());gc()

# https://blog.csdn.net/zfyyzhys/article/details/147021127
# Load pkgs---------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(Startrac)
library(SeuratWrappers)
library(clustree) # use for determine the optimal resolution
# library(ROGUE) # use for determine the optimal resolution
library(harmony)
library(stringr)
library(decontX)
# library(scDblFinder)
library(DoubletFinder)
library(Augur)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)
library(ggalluvial)
library(patchwork)
library(tidyr)
library(DESeq2)
library(RColorBrewer)
library(GSVA)
# library(slingshot) # Trajectory
# library(tradeSeq) # Trajectory
# library(monocle) # Trajectory
library(monocle3) # Trajectory
library("CellChat") # Communication
# library("iTALK") # Communication
# library(infercnv) # CNV
# library(copykat) # CNV
# library(ArchR)
# Load data---------------------------------------------------------------------
setwd("Choline/")
# options(future.globals.maxSize = 1000 * 1024^2) # 设置为1GB
set.seed(1011)

# https://github.com/friedpine/scRNASeq_ESCC/blob/main/3c.R
colors_list = c('#E76254','#EF8A47','#f4a494','#FFE6B7','#AADCE0','#528FAD',
                '#a4549c','#1E466E','#C7B8BD','#8C4834','#C17E65','#645cac',
                '#EFD061','#547857','#c49c94','#f7b6d2','#dbdb8d')

setwd("Choline/Result/scRNAseq/")
hub <- c("choline","sham")
names(hub) <- hub

Heart <- hub %>%
  Read10X %>%
  CreateSeuratObject(
    min.cells = 0,
    min.features = 0,
    project = "Heart",
    assay = "RNA")

setwd("Choline/")

dim(Heart@assays$RNA$counts)
# [1] 33696 26164

str_detect(rownames(Heart), "^ENSMUSG") %>% sum
# [1] 248
str_detect(rownames(Heart), "^Gm[0-9]+") %>% sum
# [1] 9433

scobj <- Heart[!str_detect(rownames(Heart), "^ENSMUSG"),]
scobj <- scobj[!str_detect(rownames(scobj), "^Gm[0-9]+"),]

dim(scobj@assays$RNA$counts)
# [1] 24015 26164

table(scobj@meta.data$orig.ident)
# choline    sham 
# 9673   16491
# QC----------------------------------------------------------------------------
# sum(rownames(scobj) %>% str_detect("Rp[sl]"))
# sum(rownames(scobj) %>% str_detect("Rp[0-9sl]"))
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^mt-")
scobj[["percent.rp"]] <- PercentageFeatureSet(scobj, pattern = "^Rp[sl]")
scobj[["percent.hb"]] <- PercentageFeatureSet(scobj, pattern = "^Hb[^(p)]")

pdf("Result/scRNAseq/QC1.pdf")
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
VlnPlot(scobj, features = c("nFeature_RNA"), ncol = 1)
VlnPlot(scobj, features = c("nCount_RNA"), ncol = 1)
VlnPlot(scobj, features = c("percent.mt"), ncol = 1) + scale_y_continuous(breaks = c(10,20))
VlnPlot(scobj, features = c("percent.rp"), ncol = 1)
VlnPlot(scobj, features = c("percent.hb"), ncol = 1)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
data.frame(
  nFeature = scobj@meta.data$nFeature_RNA,
  group = scobj@meta.data$orig.ident
) %>%
  ggplot(aes(x = nFeature, color = group)) +
  geom_density() +
  geom_vline(xintercept = c(200,300,400,500,1000,5000,7000,8000), color = "gray50", linetype = 2) +
  scale_x_continuous(breaks = c(200,300,400,500,750,1000,3000,4000,5000)) +
  theme_classic()

data.frame(
  nCount = scobj@meta.data$nCount_RNA,
  group = scobj@meta.data$orig.ident
) %>%
  ggplot(aes(x = nCount, color = group)) +
  geom_density() +
  scale_x_continuous(limits = c(0,5000), breaks = c(100,200,300,400,500)) +
  theme_classic()
dev.off()

scobj <- subset(
  scobj, 
  subset = 
    nFeature_RNA > 200 &
    nFeature_RNA < 3000 & 
    nCount_RNA > 200 &
    nCount_RNA < 7500 &
    percent.mt < 25 &
    percent.rp < 20 &
    percent.hb < 1)  

pdf("Result/scRNAseq/QC2.pdf")
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rp","percent.hb"), ncol = 5)
VlnPlot(scobj, features = c("nFeature_RNA"), ncol = 1)
VlnPlot(scobj, features = c("nCount_RNA"), ncol = 1)
VlnPlot(scobj, features = c("percent.mt"), ncol = 1) + scale_y_continuous(breaks = c(10,20))
VlnPlot(scobj, features = c("percent.rp"), ncol = 1)
VlnPlot(scobj, features = c("percent.hb"), ncol = 1)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
data.frame(
  nFeature = scobj@meta.data$nFeature_RNA,
  group = scobj@meta.data$orig.ident
) %>%
  ggplot(aes(x = nFeature, color = group)) +
  geom_density() +
  geom_vline(xintercept = c(200,300,400,500,1000,5000,7000,8000), color = "gray50", linetype = 2) +
  # scale_x_continuous(breaks = c(200,300,400,500,750,1000,3000,4000,5000)) +
  theme_classic()

data.frame(
  nCount = scobj@meta.data$nCount_RNA,
  group = scobj@meta.data$orig.ident
) %>%
  ggplot(aes(x = nCount, color = group)) +
  geom_density() +
  scale_x_continuous(limits = c(0,20000),
                     breaks = c(200,300,400,500,750,1000,3000,4000,5000)) +
  theme_classic()
dev.off()

table(scobj@meta.data$orig.ident)
# choline    sham 
# 7956   12938

# dim(scobj)
# [1] 24015 20894
save(scobj, file = "Result/scRNAseq/scobj.qc.Rdata")
# blacklist---------------------------------------------------------------------
blacklist <- readxl::read_excel("C:/D/R project/Multi-Omics-Routine/scRNAseq/blacklist.xlsx")
sum(blacklist$Mouse %in% rownames(scobj))
# # [1] 24

scobj <- scobj[!(rownames(scobj) %in% c(blacklist$Mouse, "Malat1")),]
scobj <- scobj[!str_detect(rownames(scobj), "^mt-"),]
scobj <- scobj[!str_detect(rownames(scobj), "^Rp[sl]"),]
scobj <- scobj[!str_detect(rownames(scobj), "^Hb[^(p)]"),]
scobj <- scobj[!str_detect(rownames(scobj), "^Hsp"),]
scobj <- scobj[!str_detect(rownames(scobj), "^Dnaj"),]
scobj <- scobj[!str_detect(rownames(scobj), "Rik$"),]

dim(scobj)
# [1] 21320 20894

scobj@meta.data %>%
  dplyr::select(orig.ident, nCount_RNA) %>%
  group_by(orig.ident) %>%
  mutate(min = min(nCount_RNA)) %>%
  distinct(min, .keep_all = TRUE)
# orig.ident nCount_RNA   min
# <fct>           <dbl> <dbl>
#   1 choline          2055   500
# 2 sham              937   606

scobj[["RNA"]] <- split(scobj[["RNA"]], f = scobj$orig.ident)
# Nomalization & harmony & cluster----------------------------------------------
# Processes including cell cycle and tissue dissociation may influence the expression of not only the associated signature genes, but also the whole transcriptome (51).
# To minimize the impact of those processes, the S phase score and G2M phase score were calculated with function CellCycleScoring, and the DIG score was calculated with function AddModuleScore. 
# Then those scores were regressed out in the Seurat pipeline.

if(1000){
  scobj.1000 <- scobj %>% 
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony") %>% 
    JoinLayers(assay = "RNA")
  
  scobj.1000.20 <- scobj.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.1000.25 <- scobj.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.1000.30 <- scobj.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Result/scRNAseq/vst.1000.h clustree 20 25 30.pdf")
  ElbowPlot(scobj.1000, reduction = "pca", ndims = 50)
  clustree(scobj.1000.20, prefix = "RNA_snn_res.")
  clustree(scobj.1000.25, prefix = "RNA_snn_res.")
  clustree(scobj.1000.30, prefix = "RNA_snn_res.")
  dev.off()
}
if(1500){
  scobj.1500 <- scobj %>% 
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony") %>% 
    JoinLayers(assay = "RNA")
  
  scobj.1500.20 <- scobj.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.1500.25 <- scobj.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.1500.30 <- scobj.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Result/scRNAseq/vst.1500.h clustree 20 25 30.pdf")
  ElbowPlot(scobj.1500, reduction = "pca", ndims = 50)
  clustree(scobj.1500.20, prefix = "RNA_snn_res.")
  clustree(scobj.1500.25, prefix = "RNA_snn_res.")
  clustree(scobj.1500.30, prefix = "RNA_snn_res.")
  dev.off()
}
if(2000){
  scobj.2000 <- scobj %>% 
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony") %>% 
    JoinLayers(assay = "RNA")
  
  scobj.2000.20 <- scobj.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.2000.25 <- scobj.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.2000.30 <- scobj.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Result/scRNAseq/vst.2000.h clustree 20 25 30.pdf")
  ElbowPlot(scobj.2000, reduction = "pca", ndims = 50)
  clustree(scobj.2000.20, prefix = "RNA_snn_res.")
  clustree(scobj.2000.25, prefix = "RNA_snn_res.")
  clustree(scobj.2000.30, prefix = "RNA_snn_res.")
  dev.off()
}

scobj.harmony <- scobj.2000.25 %>% 
  RunUMAP(reduction = "harmony", dims = 1:25)
# Visualization-----------------------------------------------------------------
scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.1")
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.1", reduction = "umap", label = T)
# doublets----------------------------------------------------------------------
save(scobj.harmony, file = "Result/scRNAseq/scobj.harmony.Rdata")
if(TRUE){
  scobj.harmony.split <- SplitObject(scobj.harmony, split.by = "orig.ident") # into list
  
  for (i in 1:length(scobj.harmony.split)) {
    # pK Identification (ground-truth) -------------------------------------------
    sweep.list <- paramSweep(scobj.harmony.split[[i]], PCs = 1:25)
    sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    
    pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)])) 
    ## Homotypic Doublet Proportion Estimate -------------------------------------
    homotypic.prop <- modelHomotypic(scobj.harmony.split[[i]]@meta.data$RNA_snn_res.0.1)  
    ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi <- round(0.1 * nrow(scobj.harmony.split[[i]]@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    ## Run DoubletFinder with varying classification stringencies ----------------
    # scobj.harmony.split[[i]] <- doubletFinder(scobj.harmony.split[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    
    scobj.harmony.split[[i]] <- doubletFinder(
      scobj.harmony.split[[i]], PCs = 1:25, pN = 0.25, pK = pK, 
      nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  }
  
  save(scobj.harmony.split, file = "Result/scRNAseq/scobj.harmony.split.Rdata")
  
  Singlet <- c(
    rownames(scobj.harmony.split[[1]]@meta.data) [scobj.harmony.split[[1]]@meta.data$DF.classifications_0.25 == "Singlet"],
    rownames(scobj.harmony.split[[2]]@meta.data) [scobj.harmony.split[[2]]@meta.data$DF.classifications_0.25 == "Singlet"])
  
  scobj.harmony@meta.data$id <- rownames(scobj.harmony@meta.data)
  
  # scobj.h.dblfinder <- subset(scobj.harmony, subset = id %in% Singlet)
  scobj.harmony@meta.data$dblfinder <- NA
  scobj.harmony@meta.data$dblfinder <- "doublet"
  scobj.harmony@meta.data$dblfinder[scobj.harmony@meta.data$id %in% Singlet] <- "singlet"
  
  dim(scobj.harmony)
  # [1] 21320 20894
  table(scobj.harmony@meta.data$dblfinder)
  # doublet singlet 
  # 1110   19784
}

pdf("/data/FYM/64.ESCC_SA_12/contamination.pdf")
CM <- c('Pln', 'Ttn', 'Tnnc2') # 1 9
ACM <- c('Tnnt2', 'Ttn', 'Myh6', 'Myh7') # 1 9
VCM <- c('Myh7', 'Myl2', 'Fhl2') # 1 9
EC <- c('Vwf', 'Cdh5', 'Thbd', 'Tek', 'Eng', 'Pecam1', 'Npr3') # 3 4 ?6 ?8
FB <- c('Fn1', 'Fbln1', 'Vim', 'Col1a1', 'Col1a2', 'Dcn') # 5
M <- c('Cd163', 'S100a8', 'Csf1r', 'C5ar1', 'Cd74') # 2
Myeloid	<- c('Cd163', 'Cd14', 'C1qa', 'C1qb', 'Adgre1', 'Ms4a7') # 2
Monocyte <- c('Vcan', 'Ccr2', 'Cd14') # 2
neutrophil <- c('Csf3r', 'S100a8', 'Retnlg') # 12
P <- c('Hcn1', 'Hcn4', 'Cacna1d') # pacemaker cells (P cells)
Peri <- c('Rgs5', 'Abcc9', 'Kcnj8') # 
SMC <- c('Myh11', 'Acta2','Cnn1', 'Tagln', 'Myocd', 'Cald1', 'Mylk', 'Notch3', 
         'Pdgfrb', 'Cspg4') #
Meso <- c('Msln', 'Wt1', 'Bnc1') # 
Neural <- c('Plp1', 'Nrxn1', 'Nrxn3') 
Adi <- c('Gpam', 'Fasn', 'Lep') # 4
Lymph <- c('Cd3e', 'Il7r', 'Cd40lg') # 7
T_cell <- c('Cd3e', 'Cd3d', 'Cd2', 'Gzmk') # 7
B_cell <- c('Cd79a', 'Ms4a1', 'Cd19') # 11
NK <- c('Gzmb', 'Rpf1') # 7
Mast <- c('Kit', 'Cpa3') 
RBC <- c('Hba1', 'Hba2', 'Hbb', 'Hbd', 'Hbe1', 'Hbg1', 'Hbg2', 'Hbm', 'Hbq1', 'Hbz') # 
CellCycle	<- c('Mki67', 'Ccna2', 'Ccnb2', 'Pcna', 'Stmn1')

DotPlot(scobj.harmony, features = c(
  CM, ACM, VCM, EC, FB, M, Myeloid, Monocyte, neutrophil, P, Peri, SMC,
  Meso, Neural, Adi, Lymph, T_cell, B_cell, NK, Mast, RBC, CellCycle
) %>% unique, 
group.by = "RNA_snn_res.0.4") + 
  scale_color_viridis() +
  labs(x = "", y = "") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
# Re-run========================================================================
scobj[['dblfinder']] <- scobj.harmony$dblfinder

save(scobj, file ="Result/scRNAseq/scobj_for_2nd.Rdata")

# 1000 without regress cell cycle
if(TRUE){
  scobj.harmony.1000 <- scobj %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony") %>% 
    JoinLayers(assay = "RNA")
  
  scobj.harmony.1000.20 <- scobj.harmony.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.harmony.1000.25 <- scobj.harmony.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.harmony.1000.30 <- scobj.harmony.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Result/scRNAseq/vst.1000.h_clustree 20 25 30(singlet) without regression.pdf")
  ElbowPlot(scobj.harmony.1000, reduction = "pca", ndims = 50)
  clustree(scobj.harmony.1000.20, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.1000.25, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.1000.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.harmony.1000,
       scobj.harmony.1000.20, scobj.harmony.1000.25, scobj.harmony.1000.30,
       file = "Result/scRNAseq/scobj.harmony.1000(without regression).Rdata")
}
# 1500 without regress cell cycle
if(TRUE){
  scobj.harmony.1500 <- scobj %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony") %>% 
    JoinLayers(assay = "RNA")
  
  scobj.harmony.1500.20 <- scobj.harmony.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.harmony.1500.25 <- scobj.harmony.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.harmony.1500.30 <- scobj.harmony.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Result/scRNAseq/vst.1500.h_clustree 20 25 30(singlet) without regression.pdf")
  ElbowPlot(scobj.harmony.1500, reduction = "pca", ndims = 50)
  clustree(scobj.harmony.1500.20, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.1500.25, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.1500.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.harmony.1500,
       scobj.harmony.1500.20, scobj.harmony.1500.25, scobj.harmony.1500.30,
       file = "Result/scRNAseq/scobj.harmony.1500(without regression).Rdata")
}
# 2000 without regress cell cycle
if(TRUE){
  scobj.harmony.2000 <- scobj %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony") %>% 
    JoinLayers(assay = "RNA")
  
  scobj.harmony.2000.20 <- scobj.harmony.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.harmony.2000.25 <- scobj.harmony.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.harmony.2000.30 <- scobj.harmony.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Result/scRNAseq/vst.2000.h_clustree 20 25 30(singlet) without regression.pdf")
  ElbowPlot(scobj.harmony.2000, reduction = "pca", ndims = 50)
  clustree(scobj.harmony.2000.20, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.2000.25, prefix = "RNA_snn_res.")
  clustree(scobj.harmony.2000.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.harmony.2000,
       scobj.harmony.2000.20, scobj.harmony.2000.25, scobj.harmony.2000.30,
       file = "Result/scRNAseq/scobj.harmony.2000(without regression).Rdata")
}

scobj.harmony <- scobj.harmony.1500.30 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)
# Visualization=================================================================
DimPlot(scobj.harmony, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.1", reduction = "umap", label = T) +
DimPlot(scobj.harmony, group.by = "RNA_snn_res.0.2", reduction = "umap", label = T)
# Annotation====================================================================
scobj.harmony <- RegroupIdents(scobj.harmony, metadata = "RNA_snn_res.0.1")
## Findmarkers------------------------------------------------------------------
Clist <- list()
n <- 1
for (i in 0:13) {
  Clist[[n]] <- FindMarkers(scobj.harmony, group.by = "RNA_snn_res.0.2", ident.1 = as.character(i))
  write.csv(Clist[[n]], paste0("D:/C",n-1,".csv"))
  
  n <- n + 1
}

res.cosg <- COSG::cosg(scobj.harmony, groups = "all")
data.frame(
  symbol = res.cosg$names$`4`,
  symbol = res.cosg$scores$`4`
)

CM <- c('Pln', 'Actc1', 'Tnni3', 'Tnnt2', 'Myl2', 'Myl3') 
ACM <- c('Tnnt2', 'Ttn', 'Myh6', 'Myh7') 
VCM <- c('Myl2', 'Fhl2')
EC <- c('Vwf', 'Cdh5', 'Thbd', 'Tek', 'Eng', 'Pecam1') 
FB <- c('Fbln1', 'Col1a1', 'Col1a2', 'Dcn') # 5
M <- c('Cd74', 'Csf1r', 'C5ar1') # 'Cd163', 
Myeloid	<- c('Cd68', 'Cd14', 'C1qa', 'C1qb', 'Adgre1', 'Ms4a7') # 2
Monocyte <- c('Ccr2', 'Cd14') # 2
neutrophil <- c('S100a8', 'S100a9', 'Csf3r', 'S100a8', 'Retnlg') # 12
P <- c('Hcn1', 'Hcn4', 'Cacna1d') # pacemaker cells (P cells)
Peri <- c('Rgs5', 'Abcc9', 'Kcnj8') # 
SMC <- c('Myh11', 'Acta2','Tagln', 'Mylk', 'Notch3', 
         'Pdgfrb') #
Meso <- c('Msln', 'Wt1', 'Bnc1') # 
Neural <- c('Plp1') 
Adi <- c('Gpam', 'Fasn', 'Lep') # 4
Lymph <- c('Cd3e', 'Il7r') # 7
T_cell <- c('Cd2', 'Cd3e', 'Cd3d') # 7
B_cell <- c('Cd19', 'Cd79a', 'Ms4a1') # 11
NK <- c('Gzmb', 'Rpf1') 
Mast <- c('Kit', 'Cpa3') 
# RBC <- c('Hba1', 'Hba2', 'Hbb', 'Hbd', 'Hbe1', 'Hbg1', 'Hbg2', 'Hbm', 'Hbq1', 'Hbz') # 
# CellCycle	<- c('Mki67', 'Ccna2', 'Ccnb2', 'Pcna', 'Stmn1')

table(scobj.harmony$RNA_snn_res.0.1)
# 0     1     2     3     4     5     6     7     8     9    10    11 
# 13371  2649  1361   808   736   303   160   127   114    63    47    45 

DotPlot(scobj.harmony, features = c(
  CM, EC, FB, M, Myeloid, Monocyte, neutrophil, Peri, SMC,
  # P, Meso, Adi, Mast
  Neural, T_cell, B_cell
  # , NK, RBC, CellCycle
) %>% unique, 
group.by = "RNA_snn_res.0.2") + 
  scale_color_viridis() +
  labs(x = "", y = "") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

VlnPlot(scobj.harmony, features = c("Dcn","Col1a1",'Fn1', 'Igfbp6', 'Mfap5','Vim', 'Pecam1',"Vwf","Eng",'Cdh5'), 
        group.by = "RNA_snn_res.0.2")
VlnPlot(scobj.harmony, features = c("Vwf","Eng",'Cdh5','Pecam1',CM), 
        group.by = "RNA_snn_res.0.2")

scobj.harmony$cell_type <- "Unknown"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(0)] <- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(1)] <- "Mesenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(2)] <- "Mesenchyme"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(3)] <- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(4)] <- "Cardiomyocyte"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(5)] <- "Mesenchyme" # 内皮间充质转化
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(6)] <- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(7)] <- "Lymphoid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(8)] <- "Myeloid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(9)] <- "Cardiomyocyte"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(10)]<- "Lymphoid"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(11)]<- "Neureonal"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(12)]<- "Endothelial"
scobj.harmony$cell_type[scobj.harmony$RNA_snn_res.0.2 %in% c(13)]<- "Endothelial"

# save(scobj.harmony, file = "Result/scRNAseq/method4/scobj.harmony(singlet-withoutreg.cc-2000-25) anno.Rdata")
load("Result/scRNAseq/method4/scobj.harmony(singlet-withoutreg.cc-2000-25) anno.Rdata")

table(scobj.harmony$RNA_snn_res.0.2)
# 0     1     2     3     4     5     6     7     8     9    10    11    12    13 
# 11137  2656  1365  1289   920   808   748   303   160   127   114    63    47    47
table(scobj.harmony$cell_type)
# Cardiomyocyte   Endothelial      Lymphoid    Mesenchyme       Myeloid     Neureonal 
# 1047         12520           417          4829           908            63
# Figure2=======================================================================
fig <- subset(scobj.harmony, subset = cell_type != "doublet")
fig$sample <- fig$orig.ident
fig$RNA_snn_res.0.2 <- factor(fig$RNA_snn_res.0.4, levels = c(0:13))
fig$cell_type <- factor(fig$cell_type, 
                        levels = c("Cardiomyocyte", "Endothelial", "Mesenchyme", 
                                   "Lymphoid", "Myeloid", "Neureonal"))
fig$orig.ident <- factor(fig$orig.ident, levels = c("sham", "choline"))

table(fig$cell_type)
# Cardiomyocyte   Endothelial    Mesenchyme      Lymphoid       Myeloid     Neureonal 
# 1097         12246          4660           401           858            63
patient_col <- c("#c79494", "#3d7d35")

p1 <- DimPlot(fig, reduction = "umap", group.by = "cell_type", cols = colors_list,
              label = TRUE, pt.size = 0.1, raster = FALSE)
p2 <- DimPlot(fig, reduction = "umap",
              group.by = "orig.ident", cols = patient_col, 
              label = TRUE, pt.size = 0.1, raster = FALSE)

p <- p1 + p2 + plot_layout(ncol = 2, nrow = 1)
ggsave(p, filename = "Figure1 UMAP.pdf", width=12, height=5, units="in")

p3 <- DotPlot(fig, features = c(
  "Pln",    "Actc1",  "Tnni3",  "Tnnt2",  "Myl2",   "Myl3", "Cdh5",
  "Thbd", "Tek",   "Eng",    "Pecam1", "Fbln1",  "Col1a1", "Col1a2", "Dcn","Rgs5",
  "Abcc9", "Kcnj8",  "Myh11",  "Acta2",  "Tagln",  "Mylk",   "Notch3", "Pdgfrb",
  "Cd2", "Cd19",   "Cd79a",  "Ms4a1",  "Csf1r",  "C5ar1",
  "Cd68","Cd14","C1qa","C1qb","C1qc","Adgre1","Ms4a7",  "Ccr2","Cd14",
  "S100a8","S100a9", "Csf3r","S100a8","Retnlg", "Plp1" ) %>% unique, group.by = "cell_type") + 
  scale_color_viridis() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p3, filename = "Figure1 MARKER.pdf", width=6, height=8, units="in")

dat <- fig@meta.data %>% 
  dplyr::select(orig.ident, cell_type)
p <- dat %>% 
  mutate(value = 1) %>% 
  group_by(orig.ident, cell_type) %>%
  summarise(n = n()) %>% 
  group_by(orig.ident) %>%
  mutate(sum = sum(n)) %>% 
  mutate(n = n / sum) %>% 
  ggplot(aes(x = orig.ident, y = n, fill = cell_type, stratum = cell_type, alluvium = cell_type))+
  geom_stratum(width = 0.5, color='white') +
  geom_alluvium(alpha = 0.5,
                width = 0.5,
                curve_type = "linear") +
  scale_color_manual(values = colors_list) +
  scale_fill_manual(values = colors_list) +
  labs(x = "", y = "Percent", fill = "Cell type") +
  theme_bw()+
  theme(text = element_text(size = 20))
ggsave(p, filename = "Figure1 percent.pdf", width=6, height=6, units="in")
# Enrichment Analysis-----------------------------------------------------------
myKEGGplot <- function(result, sizerange = c(2,8), filter = "Human Diseases", topn = 10){
  result <- result %>% 
    arrange(pvalue) %>% 
    filter(pvalue < 0.05) %>% 
    filter(category != filter)
  if(nrow(result) >= topn){
    result <- result[,1:topn]
  }
  
  p <- result %>% 
    # filter(Description %in% sele) %>% 
    mutate(n1 = GeneRatio %>% str_remove("/[0-9]*$") %>% as.numeric,
           n2 = GeneRatio %>% str_remove("^[0-9]*/") %>% as.numeric,
           GeneRatio = n1/n2) %>% 
    arrange(desc(GeneRatio)) %>% 
    mutate(Description = factor(Description, level = Description)) %>% 
    ggplot(aes(x = n1/n2, y = Description, size = Count, fill = p.adjust)) +
    geom_point(shape = 21) +
    scale_size_continuous(range = sizerange) +
    scale_fill_gradient(low = "#347eba", high = "#de6764", transform = "reverse") +
    labs(size = "counts", colour = "p.adjust", y = "", x = "GeneRatio") +
    theme_bw() +
    theme(text = element_text(size = 16))
  
  return(p)
}
## ORA==========================================================================
if(epithelial){
  load("1/1.DEA/bulk 拆分SMC/drinking/pseudobulk_Epithelial.Rdata")
  load("1/1.DEA/bulk 拆分SMC/drinking/Epithelial_DESeq2-DEGs.Rdata")
  
  sum(DEG_choline_vs_noncholine$pvalue < 0.05) # [1] 1419
  sum(DEG_choline_vs_noncholine$padj < 0.05) # [1] 67
  sum(DEG_choline_vs_noncholine$padj < 0.05 & abs(DEG_choline_vs_noncholine$log2FoldChange) > 1) # [1] 67
  
  degs <- rownames(DEG_choline_vs_noncholine)[DEG_choline_vs_noncholine$padj < 0.05 & abs(DEG_choline_vs_noncholine$log2FoldChange) > 1]
  geneList <- clusterProfiler::bitr(degs, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  
  ERA_KEGG <- enrichKEGG(
    gene = geneList$ENTREZID,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 1, 
    qvalueCutoff = 1)
  ERA_KEGG <- setReadable(ERA_KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  
  p <- myKEGGplot(ERA_KEGG@result)
  
  save(ERA_KEGG, file = "1/2.ORA/epithelial_ERA_KEGG.Rdata")
  ggsave(p, filename = "epithelial_ORA_KEGG.pdf", 
         path = "1/2.ORA/", width=6, height=3, units="in")
  
  ERA_GO <- enrichGO(
    gene = geneList$ENTREZID,
    OrgDb = "org.Hs.eg.db",
    keyType = "ENTREZID", 
    ont = "ALL",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    readable = T)
  ERA_GO <- setReadable(ERA_GO, 'org.Hs.eg.db', 'ENTREZID')
  save(ERA_GO, file = "1/2.ORA/epithelial_ERA_GO.Rdata")
  
  sum(ERA_GO@result$p.adjust < 0.05) # [1] 4
  ## KEGG=========================================================================
  load("../ESCC_rmENSG_del27/pseudobulk_Epithelial.Rdata")
  load("../ESCC_rmENSG_del27/Epithelial_DESeq2-DEGs.Rdata")
  
  DEG_nonMRPvsMRP <- DEG_nonMRPvsMRP %>% 
    tibble::rownames_to_column("SYMBOL") %>% 
    arrange(desc(log2FoldChange)) %>% 
    left_join(geneList, by = "SYMBOL") %>% 
    filter(complete.cases(ENTREZID))
  
  geneList2 <- DEG_nonMRPvsMRP$log2FoldChange
  names(geneList2) <- DEG_nonMRPvsMRP$ENTREZID
  GSEA_KEGG <- gseKEGG(geneList = geneList2, 
                       seed = 1011,
                       organism = "hsa",
                       minGSSize    = 10,
                       pvalueCutoff = 1,
                       verbose      = FALSE, 
                       eps = 0)
  GSEA_KEGG <- setReadable(GSEA_KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  save(GSEA_KEGG, file = "../ESCC_rmENSG_del27/GSEA_KEGG.Rdata")
  sum(GSEA_KEGG@result$p.adjust < 0.05) # [1] 2
  
  p <- dotplot(GSEA_KEGG)
  ggsave(p, filename = "Figure1-Enrichment Analysis-GSEA KEGG.pdf", 
         path = "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/", width=10, height=10, units="in")
  export::table2excel(GSEA_KEGG@result, "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/GSEA_KEGG.xlsx")
  
  GSEA_GO <- gseGO(geneList = geneList2,
                   ont = "ALL",
                   OrgDb = "org.Hs.eg.db",
                   minGSSize    = 10,
                   pvalueCutoff = 1,
                   verbose      = FALSE, 
                   eps = 0)
  GSEA_GO <- setReadable(GSEA_GO, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  save(GSEA_GO, file = "../ESCC_rmENSG_del27/GSEA_GO.Rdata")
  sum(GSEA_GO@result$p.adjust < 0.05) # [1] 136
  export::table2excel(GSEA_GO@result, "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/GSEA_GO.xlsx")
  
  p <- dotplot(GSEA_GO)
  ggsave(p, filename = "Figure1-Enrichment Analysis-GSEA GO.pdf", 
         path = "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/", width=10, height=10, units="in")
  
  # mytheme <- theme_bw() +
  #   theme(text = element_text(size = 12), 
  #         plot.margin = ggplot2::margin(30,30,30,30), 
  #         panel.grid = element_blank(), 
  #         legend.title = element_blank(),
  #         legend.background = element_rect(linetype = 1, colour = "#555555"),
  #         axis.title.x = element_text(vjust = -4),
  #         axis.title.y = element_text(vjust = 4),
  #         legend.position = "right")
  # gg_GOE_bar <- function(NES, Topn = 10){
  #   Top <- NES@result
  #   
  #   TopBPUp <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "BP" & NES > 0) %>%
  #     arrange(desc(NES))
  #   
  #   TopBPDown <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "BP" & NES < 0) %>%
  #     arrange(NES)
  #   
  #   TopCCUp <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "CC" & NES > 0) %>%
  #     arrange(desc(NES))
  #   
  #   TopCCDown <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "CC" & NES < 0) %>%
  #     arrange(NES)
  #   
  #   TopMFUp <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "MF" & NES > 0) %>%
  #     arrange(desc(NES))
  #   
  #   TopMFDown <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "MF" & NES < 0) %>%
  #     arrange(NES)
  #   
  #   GO_EA2 <- NES
  #   GO_EA2@result <- GO_EA2@result %>% 
  #     filter(Description %in% c(TopBPUp$Description[1:Topn], TopBPDown$Description[1:Topn],
  #                               TopCCUp$Description[1:Topn], TopCCDown$Description[1:Topn],
  #                               TopMFUp$Description[1:Topn], TopMFDown$Description[1:Topn]))
  #   
  #   tmp <- GO_EA2@result %>% arrange(NES) %>%
  #     mutate(ID = factor(ID, levels = ID),
  #            color = ifelse(NES < 0, 0, 1) %>% factor(labels = c("down","up")))
  #   ggplot(tmp, aes(x = NES, y = ID, fill = color)) +
  #     geom_vline(xintercept = 0, color = "gray80") +
  #     geom_col(width = 0.8) +
  #     scale_fill_manual(values = c("#448844","#AB3A29")) +
  #     labs(y = "") +
  #     guides(fill = "none") +
  #     facet_grid(ONTOLOGY~., scale='free') + 
  #     mytheme
  # }
  # gg_GOE_bar(GSEA_GO)
  ## Hallmark=====================================================================
  require(GSVA)
  
  load("1/1.DEA/bulk 拆分SMC/drinking/pseudobulk_Epithelial.Rdata")
  
  counts <- counts %>% dplyr::select(
    all_of(c(c("sc-04","sc-06","sc-07","sc-10"),
             c("sc-01","sc-02","sc-03","sc-05","sc-08","sc-09","sc-11","sc-12")))
  )
  
  hallmarks <- qusage::read.gmt("C:/D/Project/63.ESCC/h.all.v2023.1.Hs.symbols.gmt") #返回的是表格
  gsvaPar <- ssgseaParam(exprData = as.matrix(counts), 
                         geneSets = hallmarks,
                         normalize = TRUE)
  gsva_data <- gsva(gsvaPar, verbose = FALSE)
  
  annotation_cols <- data.frame(
    Group = c(rep(0,4), rep(1,8)) %>% 
      factor(levels = c(0,1), labels = c("noncholine","choline")),
    row.names = colnames(gsva_data)
  )
  pheatmap::pheatmap(gsva_data, 
                     clustering_method = "ward.D",
                     annotation_col = annotation_cols,
                     annotation_colors = list(
                       Group = c(choline = "#C00000", noncholine = "black")),
                     cluster_rows = TRUE,
                     cluster_cols = FALSE,
                     scale = "row")
}
if(myeloid){
  load("../ESCC_rmENSG_del27/Figures/Figure1 DEA/Myeloid_DESeq2-DEGs.Rdata")
  load("../ESCC_rmENSG_del27/Figures/Figure1 DEA/pseudobulk_Myeloid.Rdata")
  
  sum(DEG_nonMRPvsMRP$pvalue < 0.05) # [1] 1819
  sum(DEG_nonMRPvsMRP$padj < 0.05) # [1] 34
  sum(DEG_nonMRPvsMRP$padj < 0.05 & abs(DEG_nonMRPvsMRP$log2FoldChange) > 1) # [1] 29
  
  degs <- rownames(DEG_nonMRPvsMRP)[DEG_nonMRPvsMRP$padj < 0.05 & abs(DEG_nonMRPvsMRP$log2FoldChange) > 1]
  geneList <- clusterProfiler::bitr(degs, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  
  ERA_KEGG <- enrichKEGG(
    gene = geneList$ENTREZID,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 1, 
    qvalueCutoff = 1)
  ERA_KEGG <- setReadable(ERA_KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  save(ERA_KEGG, file = "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/Myeloid_ERA_KEGG.Rdata")
  sum(ERA_KEGG@result$p.adjust < 0.05) # 10
  p <- dotplot(ERA_KEGG, showCategory = sum(ERA_KEGG@result$p.adjust < 0.05))+ 
    scale_y_discrete(labels=function(x) str_wrap(x, width=70))
  
  ggsave(p, filename = "Myeloid Enrichment Analysis-ORA KEGG.pdf", 
         path = "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/", width=10, height=10, units="in")
  export::table2excel(ERA_KEGG@result, "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/Myeloid_ERA_KEGG.xlsx")
  
  ERA_GO <- enrichGO(
    gene = geneList$ENTREZID,
    OrgDb = "org.Hs.eg.db",
    keyType = "ENTREZID", 
    ont = "ALL",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    readable = T)
  ERA_GO <- setReadable(ERA_GO, 'org.Hs.eg.db', 'ENTREZID')
  save(ERA_GO, file = "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/Myeloid_ERA_GO.Rdata")
  export::table2excel(ERA_GO@result, "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/Myeloid_ERA_GO.xlsx")
  
  sum(ERA_GO@result$p.adjust < 0.05) # [1] 31
  
  BPtop10 <- ERA_GO@result %>% arrange(p.adjust) %>% filter(ONTOLOGY == "BP");BPtop10 <- BPtop10$ID[1:10]
  CCtop10 <- ERA_GO@result %>% arrange(p.adjust) %>% filter(ONTOLOGY == "CC");CCtop10 <- CCtop10$ID[1:10]
  MFtop10 <- ERA_GO@result %>% arrange(p.adjust) %>% filter(ONTOLOGY == "MF");MFtop10 <- MFtop10$ID[1:10]
  
  ERA_GO@result <- ERA_GO@result %>% 
    filter(ID %in% c(BPtop10, CCtop10, MFtop10)) 
  p <- dotplot(ERA_GO, showCategory = 31) + 
    # facet_grid(ONTOLOGY~., scale='free') + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50))
  ggsave(p, filename = "Myeloid Enrichment Analysis-ORA GO.pdf", 
         path = "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/", width=10, height=10, units="in")
  ## KEGG=========================================================================
  # load("../ESCC_rmENSG_del27/pseudobulk_Epithelial.Rdata")
  # load("../ESCC_rmENSG_del27/Epithelial_DESeq2-DEGs.Rdata")
  
  DEG_nonMRPvsMRP <- DEG_nonMRPvsMRP %>% 
    tibble::rownames_to_column("SYMBOL") %>% 
    arrange(desc(log2FoldChange)) %>% 
    left_join(geneList, by = "SYMBOL") %>% 
    filter(complete.cases(ENTREZID))
  
  geneList2 <- DEG_nonMRPvsMRP$log2FoldChange
  names(geneList2) <- DEG_nonMRPvsMRP$ENTREZID
  GSEA_KEGG <- gseKEGG(geneList = geneList2, 
                       seed = 1011,
                       organism = "hsa",
                       minGSSize    = 10,
                       pvalueCutoff = 1,
                       verbose      = FALSE, 
                       eps = 0)
  GSEA_KEGG <- setReadable(GSEA_KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  save(GSEA_KEGG, file = "../ESCC_rmENSG_del27/GSEA_KEGG.Rdata")
  sum(GSEA_KEGG@result$p.adjust < 0.05) # [1] 0
  
  p <- dotplot(GSEA_KEGG)
  ggsave(p, filename = "Figure1-Enrichment Analysis-GSEA KEGG.pdf", 
         path = "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/", width=10, height=10, units="in")
  export::table2excel(GSEA_KEGG@result, "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/GSEA_KEGG.xlsx")
  
  GSEA_GO <- gseGO(geneList = geneList2,
                   ont = "ALL",
                   OrgDb = "org.Hs.eg.db",
                   minGSSize    = 10,
                   pvalueCutoff = 1,
                   verbose      = FALSE, 
                   eps = 0)
  GSEA_GO <- setReadable(GSEA_GO, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  save(GSEA_GO, file = "../ESCC_rmENSG_del27/GSEA_GO.Rdata")
  sum(GSEA_GO@result$p.adjust < 0.05) # [1] 0
  export::table2excel(GSEA_GO@result, "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/GSEA_GO.xlsx")
  
  p <- dotplot(GSEA_GO)
  ggsave(p, filename = "Figure1-Enrichment Analysis-GSEA GO.pdf", 
         path = "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/", width=10, height=10, units="in")
  
  # mytheme <- theme_bw() +
  #   theme(text = element_text(size = 12), 
  #         plot.margin = ggplot2::margin(30,30,30,30), 
  #         panel.grid = element_blank(), 
  #         legend.title = element_blank(),
  #         legend.background = element_rect(linetype = 1, colour = "#555555"),
  #         axis.title.x = element_text(vjust = -4),
  #         axis.title.y = element_text(vjust = 4),
  #         legend.position = "right")
  # gg_GOE_bar <- function(NES, Topn = 10){
  #   Top <- NES@result
  #   
  #   TopBPUp <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "BP" & NES > 0) %>%
  #     arrange(desc(NES))
  #   
  #   TopBPDown <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "BP" & NES < 0) %>%
  #     arrange(NES)
  #   
  #   TopCCUp <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "CC" & NES > 0) %>%
  #     arrange(desc(NES))
  #   
  #   TopCCDown <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "CC" & NES < 0) %>%
  #     arrange(NES)
  #   
  #   TopMFUp <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "MF" & NES > 0) %>%
  #     arrange(desc(NES))
  #   
  #   TopMFDown <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "MF" & NES < 0) %>%
  #     arrange(NES)
  #   
  #   GO_EA2 <- NES
  #   GO_EA2@result <- GO_EA2@result %>% 
  #     filter(Description %in% c(TopBPUp$Description[1:Topn], TopBPDown$Description[1:Topn],
  #                               TopCCUp$Description[1:Topn], TopCCDown$Description[1:Topn],
  #                               TopMFUp$Description[1:Topn], TopMFDown$Description[1:Topn]))
  #   
  #   tmp <- GO_EA2@result %>% arrange(NES) %>%
  #     mutate(ID = factor(ID, levels = ID),
  #            color = ifelse(NES < 0, 0, 1) %>% factor(labels = c("down","up")))
  #   ggplot(tmp, aes(x = NES, y = ID, fill = color)) +
  #     geom_vline(xintercept = 0, color = "gray80") +
  #     geom_col(width = 0.8) +
  #     scale_fill_manual(values = c("#448844","#AB3A29")) +
  #     labs(y = "") +
  #     guides(fill = "none") +
  #     facet_grid(ONTOLOGY~., scale='free') + 
  #     mytheme
  # }
  # gg_GOE_bar(GSEA_GO)
  ## Hallmark=====================================================================
  require(GSVA)
  
  load("1/1.DEA/bulk 拆分SMC/drinking/pseudobulk_Myeloid.Rdata")
  
  counts <- counts %>% dplyr::select(
    all_of(c(c("sc-04","sc-06","sc-07","sc-10"),
             c("sc-01","sc-02","sc-03","sc-05","sc-08","sc-09","sc-11","sc-12")))
  )
  
  hallmarks <- qusage::read.gmt("C:/D/Project/63.ESCC/h.all.v2023.1.Hs.symbols.gmt") #返回的是表格
  gsvaPar <- ssgseaParam(exprData = as.matrix(counts), 
                         geneSets = hallmarks,
                         normalize = TRUE)
  gsva_data <- gsva(gsvaPar, verbose = FALSE)
  
  annotation_cols <- data.frame(
    Group = c(rep(0,4), rep(1,8)) %>% 
      factor(levels = c(0,1), labels = c("noncholine","choline")),
    row.names = colnames(gsva_data)
  )
  pheatmap::pheatmap(gsva_data, 
                     clustering_method = "ward.D",
                     annotation_col = annotation_cols,
                     annotation_colors = list(
                       Group = c(choline = "#C00000", noncholine = "black")),
                     cluster_rows = TRUE,
                     cluster_cols = FALSE,
                     scale = "row")
}
if(mesenchymal){
  load("../ESCC_rmENSG_del27/Figures/Figure1 DEA/T cell_DESeq2-DEGs.Rdata")
  load("../ESCC_rmENSG_del27/Figures/Figure1 DEA/pseudobulk_T cell.Rdata")
  
  sum(DEG_nonMRPvsMRP$pvalue < 0.05) # [1] 1582
  sum(DEG_nonMRPvsMRP$padj < 0.05) # [1] 38
  sum(DEG_nonMRPvsMRP$padj < 0.05 & abs(DEG_nonMRPvsMRP$log2FoldChange) > 1) # [1] 35
  
  degs <- rownames(DEG_nonMRPvsMRP)[DEG_nonMRPvsMRP$padj < 0.05 & abs(DEG_nonMRPvsMRP$log2FoldChange) > 1]
  geneList <- clusterProfiler::bitr(degs, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  
  ERA_KEGG <- enrichKEGG(
    gene = geneList$ENTREZID,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 1, 
    qvalueCutoff = 1)
  ERA_KEGG <- setReadable(ERA_KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  save(ERA_KEGG, file = "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/T cell_ERA_KEGG.Rdata")
  sum(ERA_KEGG@result$p.adjust < 0.05) # 1
  p <- dotplot(ERA_KEGG, showCategory = 10)+ 
    scale_y_discrete(labels=function(x) str_wrap(x, width=70))
  
  ggsave(p, filename = "T cell_Enrichment Analysis-ORA KEGG.pdf", 
         path = "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/", width=10, height=10, units="in")
  export::table2excel(ERA_KEGG@result, "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/T cell_ERA_KEGG.xlsx")
  
  ERA_GO <- enrichGO(
    gene = geneList$ENTREZID,
    OrgDb = "org.Hs.eg.db",
    keyType = "ENTREZID", 
    ont = "ALL",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    readable = T)
  ERA_GO <- setReadable(ERA_GO, 'org.Hs.eg.db', 'ENTREZID')
  save(ERA_GO, file = "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/T cell_ERA_GO.Rdata")
  export::table2excel(ERA_GO@result, "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/T cell_ERA_GO.xlsx")
  
  sum(ERA_GO@result$p.adjust < 0.05) # [1] 67
  
  BPtop10 <- ERA_GO@result %>% arrange(p.adjust) %>% filter(ONTOLOGY == "BP");BPtop10 <- BPtop10$ID[1:10]
  CCtop10 <- ERA_GO@result %>% arrange(p.adjust) %>% filter(ONTOLOGY == "CC");CCtop10 <- CCtop10$ID[1:10]
  MFtop10 <- ERA_GO@result %>% arrange(p.adjust) %>% filter(ONTOLOGY == "MF");MFtop10 <- MFtop10$ID[1:10]
  
  ERA_GO@result <- ERA_GO@result %>% 
    filter(ID %in% c(BPtop10, CCtop10, MFtop10)) 
  p <- dotplot(ERA_GO, showCategory = 67) + 
    # facet_grid(ONTOLOGY~., scale='free') + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=80))
  ggsave(p, filename = "T cell_Enrichment Analysis-ORA GO.pdf", 
         path = "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/", width=10, height=10, units="in")
  ## KEGG=========================================================================
  # load("../ESCC_rmENSG_del27/pseudobulk_Epithelial.Rdata")
  # load("../ESCC_rmENSG_del27/Epithelial_DESeq2-DEGs.Rdata")
  
  DEG_nonMRPvsMRP <- DEG_nonMRPvsMRP %>% 
    tibble::rownames_to_column("SYMBOL") %>% 
    arrange(desc(log2FoldChange)) %>% 
    left_join(geneList, by = "SYMBOL") %>% 
    filter(complete.cases(ENTREZID))
  
  geneList2 <- DEG_nonMRPvsMRP$log2FoldChange
  names(geneList2) <- DEG_nonMRPvsMRP$ENTREZID
  GSEA_KEGG <- gseKEGG(geneList = geneList2, 
                       seed = 1011,
                       organism = "hsa",
                       minGSSize    = 10,
                       pvalueCutoff = 1,
                       verbose      = FALSE, 
                       eps = 0)
  GSEA_KEGG <- setReadable(GSEA_KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  save(GSEA_KEGG, file = "../ESCC_rmENSG_del27/GSEA_KEGG.Rdata")
  sum(GSEA_KEGG@result$p.adjust < 0.05) # [1] 0
  
  p <- dotplot(GSEA_KEGG)
  ggsave(p, filename = "Figure1-Enrichment Analysis-GSEA KEGG.pdf", 
         path = "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/", width=10, height=10, units="in")
  export::table2excel(GSEA_KEGG@result, "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/GSEA_KEGG.xlsx")
  
  GSEA_GO <- gseGO(geneList = geneList2,
                   ont = "ALL",
                   OrgDb = "org.Hs.eg.db",
                   minGSSize    = 10,
                   pvalueCutoff = 1,
                   verbose      = FALSE, 
                   eps = 0)
  GSEA_GO <- setReadable(GSEA_GO, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  save(GSEA_GO, file = "../ESCC_rmENSG_del27/GSEA_GO.Rdata")
  sum(GSEA_GO@result$p.adjust < 0.05) # [1] 0
  export::table2excel(GSEA_GO@result, "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/GSEA_GO.xlsx")
  
  p <- dotplot(GSEA_GO)
  ggsave(p, filename = "Figure1-Enrichment Analysis-GSEA GO.pdf", 
         path = "C:/D/Project/63.ESCC/ESCC_rmENSG_del27/", width=10, height=10, units="in")
  
  # mytheme <- theme_bw() +
  #   theme(text = element_text(size = 12), 
  #         plot.margin = ggplot2::margin(30,30,30,30), 
  #         panel.grid = element_blank(), 
  #         legend.title = element_blank(),
  #         legend.background = element_rect(linetype = 1, colour = "#555555"),
  #         axis.title.x = element_text(vjust = -4),
  #         axis.title.y = element_text(vjust = 4),
  #         legend.position = "right")
  # gg_GOE_bar <- function(NES, Topn = 10){
  #   Top <- NES@result
  #   
  #   TopBPUp <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "BP" & NES > 0) %>%
  #     arrange(desc(NES))
  #   
  #   TopBPDown <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "BP" & NES < 0) %>%
  #     arrange(NES)
  #   
  #   TopCCUp <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "CC" & NES > 0) %>%
  #     arrange(desc(NES))
  #   
  #   TopCCDown <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "CC" & NES < 0) %>%
  #     arrange(NES)
  #   
  #   TopMFUp <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "MF" & NES > 0) %>%
  #     arrange(desc(NES))
  #   
  #   TopMFDown <- Top %>% 
  #     filter(p.adjust < .05 & ONTOLOGY == "MF" & NES < 0) %>%
  #     arrange(NES)
  #   
  #   GO_EA2 <- NES
  #   GO_EA2@result <- GO_EA2@result %>% 
  #     filter(Description %in% c(TopBPUp$Description[1:Topn], TopBPDown$Description[1:Topn],
  #                               TopCCUp$Description[1:Topn], TopCCDown$Description[1:Topn],
  #                               TopMFUp$Description[1:Topn], TopMFDown$Description[1:Topn]))
  #   
  #   tmp <- GO_EA2@result %>% arrange(NES) %>%
  #     mutate(ID = factor(ID, levels = ID),
  #            color = ifelse(NES < 0, 0, 1) %>% factor(labels = c("down","up")))
  #   ggplot(tmp, aes(x = NES, y = ID, fill = color)) +
  #     geom_vline(xintercept = 0, color = "gray80") +
  #     geom_col(width = 0.8) +
  #     scale_fill_manual(values = c("#448844","#AB3A29")) +
  #     labs(y = "") +
  #     guides(fill = "none") +
  #     facet_grid(ONTOLOGY~., scale='free') + 
  #     mytheme
  # }
  # gg_GOE_bar(GSEA_GO)
  ## Hallmark=====================================================================
  require(GSVA)
  
  load("1/1.DEA/bulk 拆分SMC/drinking/pseudobulk_Mesenchyme.Rdata")
  
  counts <- counts %>% dplyr::select(
    all_of(c(c("sc-04","sc-06","sc-07","sc-10"),
             c("sc-01","sc-02","sc-03","sc-05","sc-08","sc-09","sc-11","sc-12")))
  )
  
  hallmarks <- qusage::read.gmt("C:/D/Project/63.ESCC/h.all.v2023.1.Hs.symbols.gmt") #返回的是表格
  gsvaPar <- ssgseaParam(exprData = as.matrix(counts), 
                         geneSets = hallmarks,
                         normalize = TRUE)
  gsva_data <- gsva(gsvaPar, verbose = FALSE)
  
  annotation_cols <- data.frame(
    Group = c(rep(0,4), rep(1,8)) %>% 
      factor(levels = c(0,1), labels = c("noncholine","choline")),
    row.names = colnames(gsva_data)
  )
  pheatmap::pheatmap(gsva_data, 
                     clustering_method = "ward.D",
                     annotation_col = annotation_cols,
                     annotation_colors = list(
                       Group = c(choline = "#C00000", noncholine = "black")),
                     cluster_rows = TRUE,
                     cluster_cols = FALSE,
                     scale = "row")
  
}
# Roe & Augur & MiloR-----------------------------------------------------------
# Roe----
data <- fig@meta.data
data$majorCluster = data$cell_type
data$patient = data$sample
data$loc = data$orig.ident
Roe <- calTissueDist(data,
                     byPatient = F,
                     colname.cluster = "majorCluster", # 不同细胞亚群
                     colname.patient = "patient", # 不同样本
                     colname.tissue = "loc", # 不同组织
                     method = "chisq", # "chisq", "fisher", and "freq" 
                     min.rowSum = 0) 
# sham   choline
# Cardiomyocyte 1.0340044 0.9439089
# Endothelial   1.0671458 0.8892414
# Mesenchyme    0.8525121 1.2432849
# Lymphoid      0.8551181 1.2389863
# Myeloid       0.8915667 1.1788634
# Neureonal     0.9178487 1.1355105

p_roe <- Roe %>% as.data.frame %>% 
  mutate(
    Var1 = factor(Var1, levels = c("Cardiomyocyte", "Endothelial", "Mesenchyme", 
                                   "Lymphoid", "Myeloid", "Neureonal")),
    Var2 = factor(Var2, levels = c("sham","choline"))) %>%
  mutate(label = format(Freq, digits = 2)) %>%
  ggplot(aes(x = Var2, y = Var1, fill = Freq, label = label)) +
  geom_tile() +
  geom_text() +
  labs(x = "", y = "") +
  scale_fill_gradient2(low = "#30bbea", mid = "white", high = "#e73a36", midpoint = 1) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_classic() +
  theme(axis.line = element_blank(), axis.ticks = element_blank())

ggsave(p_roe, filename = "Figure1 Roe.pdf", width=4, height=4.5, units="in")
# Augur----
# augur <- calculate_auc(fig, cell_type_col = "cell_type", label_col = "orig.ident")
# augur$AUC
# 
# save(augur, file = "Result/scRNAseq/method3/augur.Rdata")
#===============================================================================
# SECTION II - Cell Cell Communication Comparison===============================
#===============================================================================
# Part I: Data input & processing and initialization of CellChat object---------
# Load data---------------------------------------------------------------------
fig <- subset(scobj.harmony, subset = cell_type != "doublet")
fig$sample <- fig$orig.ident
fig$RNA_snn_res.0.2 <- factor(fig$RNA_snn_res.0.4, levels = c(0:12))
fig$cell_type <- factor(fig$cell_type, 
                        levels = c("Cardiomyocyte", "Endothelial", "Mesenchyme", 
                                   "Lymphoid", "Myeloid", "Neureonal"))
fig$orig.ident <- factor(fig$orig.ident, levels = c("sham", "choline"))

table(fig@meta.data$cell_type)
# Cardiomyocyte   Endothelial    Mesenchyme      Lymphoid       Myeloid     Neureonal 
# 1097         12246          4660           401           858            63 

fig.sham <- subset(fig, subset = sample == "sham")
fig.chol <- subset(fig, subset = sample == "choline")
#### 1.Create cellchat obj------------------------------------------------------
library(CellChat)
cc.sham <-createCellChat(object = fig.sham, group.by = "cell_type", assay = "RNA")
cc.chol <-createCellChat(object = fig.chol, group.by = "cell_type", assay = "RNA")
#### 2.Load LR DB---------------------------------------------------------------
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- subsetDB(CellChatDB)

cc.sham@DB <- CellChatDB.use
cc.chol@DB <- CellChatDB.use
#### 3.Preparement--------------------------------------------------------------
cc.sham <- subsetData(cc.sham) # This step is necessary even if using the whole database
cc.chol <- subsetData(cc.chol)

cc.sham <- identifyOverExpressedGenes(cc.sham)
cc.sham <- identifyOverExpressedInteractions(cc.sham)

cc.chol<- identifyOverExpressedGenes(cc.chol)
cc.chol<- identifyOverExpressedInteractions(cc.chol)

cc.sham <- smoothData(cc.sham, adj = PPI.mouse) # 2.1.2 adj = is forced
cc.chol<- smoothData(cc.chol, adj = PPI.mouse) # 2.1.2 adj = is forced
# Part II: Inference of cell-cell communication network-------------------------
cc.sham <- computeCommunProb(cc.sham, 
                             raw.use = FALSE, 
                             population.size = TRUE,
                             type =  "truncatedMean", 
                             trim = 0.2) 

cc.chol<- computeCommunProb(cc.chol, 
                            raw.use = FALSE, 
                            population.size = TRUE,
                            type =  "truncatedMean", 
                            trim = 0.2) 

# save(cc.sham, cc.chol, file = "Result/scRNAseq/method3/2. Cellchat/bulk-比较-细胞通讯.Rdata")
load("Result/scRNAseq/method3/2. Cellchat/bulk-比较-细胞通讯.Rdata")

cc.sham <- filterCommunication(cc.sham, min.cells = 10)
cc.sham <- computeCommunProbPathway(cc.sham)
cc.chol<- filterCommunication(cc.chol, min.cells = 10)
cc.chol<- computeCommunProbPathway(cc.chol)

cc.sham <- aggregateNet(cc.sham)
cc.chol<- aggregateNet(cc.chol)

cc.sham <- netAnalysis_computeCentrality(cc.sham, slot.name = "netP")
cc.chol<- netAnalysis_computeCentrality(cc.chol, slot.name = "netP")

cc.sham <- computeNetSimilarity(cc.sham, type = "functional")
cc.sham <- netEmbedding(cc.sham, type = "functional")
cc.sham <- netClustering(cc.sham, type = "functional")
cc.chol<- computeNetSimilarity(cc.chol, type = "functional")
cc.chol<- netEmbedding(cc.chol, type = "functional")
cc.chol<- netClustering(cc.chol, type = "functional")

object.list <- list(sham = cc.sham, choline = cc.chol)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat
# An object of class CellChat created from a merged object with multiple datasets 
# 1237 signaling genes.
# 19325 cells. 
# CellChat analysis of single cell RNA-seq data! 
# Part I: Identify altered interactions and cell populations
# CellChat employs a top-down approach, i.e., starting with the big picture and then refining it in a greater detail on the signaling mechanisms, to identify signaling changes at different levels, including altered interactions, cell populations, signaling pathways and ligand-receptor pairs. First, CellChat starts with a big picture to answer the following questions:
# Whether the cell-cell communication is enhanced or not
# The interaction between which cell types is significantly changed
# How the major sources and targets change from one condition to another
# Compare the total number of interactions and interaction strength
# To answer the question on whether the cell-cell communication is enhanced or not, CellChat compares the total number of interactions and interaction strength of the inferred cell-cell communication networks from different biological conditions.

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

# Compare the number of interactions and interaction strength among different cell populations
# To identify the interaction between which cell populations showing significant changes, CellChat compares the number of interactions and interaction strength among different cell populations using a circle plot with differential interactions (option A), a heatmap with differential interactions (option B) and two circle plots with the number of interactions or interaction strength per dataset (option C). Alternatively, users can examine the differential number of interactions or interaction strength among coarse cell types by aggregating the cell-cell communication based on the defined cell groups (option D).
# 
# (A) Circle plot showing differential number of interactions or interaction strength among different cell populations across two datasets
# The differential number of interactions or interaction strength in the cell-cell communication network between two datasets can be visualized using circle plot, where red
# (or blue
# ) colored edges represent increased
# (or decreased
# ) signaling in the second dataset compared to the first one.

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# (B) Heatmap showing differential number of interactions or interaction strength among different cell populations across two datasets
# CellChat can also show differential number of interactions or interaction strength in greater details using a heatmap. The top colored bar plot represents the sum of each column of the absolute values displayed in the heatmap (incoming signaling). The right colored bar plot represents the sum of each row of the absolute values (outgoing signaling). Therefore, the bar height indicates the degree of change in terms of the number of interactions or interaction strength between the two conditions. In the colorbar, red
# (or blue
# ) represents increased
# (or decreased
# ) signaling in the second dataset compared to the first one.

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

# Compare the major sources and targets in a 2D space
# Comparing the outgoing and incoming interaction strength in a 2D space allows ready identification of the cell populations with significant changes in sending or receiving signals between different datasets.
# 
# Identify cell populations with significant changes in sending or receiving signals between different datasets by following option A, and the signaling changes of specific cell populations by following option B.
# 
# (A) Identify cell populations with significant changes in sending or receiving signals
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
# gg <- list()
# for (i in 1:length(object.list)) {
#   gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
# }
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
p1 <- netAnalysis_signalingRole_scatter(object.list[[1]], title = names(object.list)[1], weight.MinMax = weight.MinMax)
p2 <- netAnalysis_signalingRole_scatter(object.list[[2]], title = names(object.list)[2], weight.MinMax = weight.MinMax)
p1 <- p1 + scale_x_continuous(limits = c(0,0.4)) + scale_y_continuous(limits = c(0,0.25))
p2 <- p2 + scale_x_continuous(limits = c(0,0.4)) + scale_y_continuous(limits = c(0,0.25))
p1 + p2

# (B) Identify the signaling changes of specific cell populations
# Furthermore, we can identify the specific signaling changes of Inflam.DC and cDC1 between NL and LS.
# 
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Endothelial")
# #> Visualizing differential outgoing and incoming signaling changes from NL to LS
# #> The following `from` values were not present in `x`: 0
# #> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Endothelial", signaling.exclude = c("MIF"))
# #> Visualizing differential outgoing and incoming signaling changes from NL to LS
# #> The following `from` values were not present in `x`: 0, 2
# #> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2))

# Part II: Identify altered signaling with distinct network architecture and interaction strength
# Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

# Compute and visualize the pathway distance in the learned joint manifold
# CellChat can identify the signaling networks with larger (or smaller) difference based on their Euclidean distance in the shared two-dimensions space. Larger distance implies larger difference of the communication networks between two datasets in terms of either functional or structure similarity. It should be noted that we only compute the distance of overlapped signaling pathways between two datasets. Those signaling pathways that are only identified in one dataset are not considered here. If there are more than three datasets, one can do pairwise comparisons by modifying the parameter comparison in the function rankSimilarity.

rankSimilarity(cellchat, type = "functional") # the result are disturbed by seed
#> Compute the distance of signaling networks between datasets 1 2

# Identify altered signaling with distinct interaction strength=================
# By comparing the information flow/interaction strength of each signaling pathway, CellChat identifies signaling pathways that: (i) turn off, (ii) decrease, (iii) turn on, or (iv) increase, by changing their information flow at one condition as compared to another condition. Identify the altered signaling pathways or ligand-receptor pairs based on the overall information flow by following option A, and based on the outgoing (or incoming) signaling patterns by following option B.

# (A) Compare the overall information flow of each signaling pathway or ligand-receptor pair
# CellChat can identify the conserved and context-specific signaling pathways by simply comparing the information flow for each signaling pathway, which is defined by the sum of communication probability among all pairs of cell groups in the inferred network (i.e., the total weights in the network).

# This bar chart can be plotted in a stacked mode or not. Significant signaling pathways were ranked based on differences in the overall information flow within the inferred networks between NL and LS skin. When setting do.stat = TRUE, a paired Wilcoxon test is performed to determine whether there is a significant difference of the signaling information flow between two conditions. The top signaling pathways colored red are enriched in NL skin, and these colored greens are enriched in the LS skin.
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2

# (B) Compare outgoing (or incoming) signaling patterns associated with each cell population
# The above analysis summarize the information from the outgoing and incoming signaling together. CellChat can also compare the outgoing (or incoming) signaling pattern between two datasets, allowing to identify signaling pathways/ligand-receptors that exhibit different signaling patterns. We can combine all the identified signaling pathways from different datasets and thus compare them side by side, including outgoing signaling, incoming signaling and overall signaling by aggregating outgoing and incoming signaling together.

# CellChat uses heatmap plot to show the contribution of signals (signaling pathways or ligand-receptor pairs) to cell groups in terms of outgoing or incoming signaling. In this heatmap, colobar represents the relative signaling strength of a signaling pathway across cell groups (Note that values are row-scaled). The top colored bar plot shows the total signaling strength of a cell group by summarizing all signaling pathways displayed in the heatmap. The right grey bar plot shows the total signaling strength of a signaling pathway by summarizing all cell groups displayed in the heatmap.
library(ComplexHeatmap)
if(incoming){
  if(TRUE){
    object <- object.list[[1]]
    pattern = c("outgoing", "incoming", "all")[2]
    slot.name = "netP"
    
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
      outgoing[, i] <- centr[[i]]$outdeg
      incoming[, i] <- centr[[i]]$indeg
    }
    if (pattern == "outgoing") {
      mat <- t(outgoing)
    } else if (pattern == "incoming") {
      mat <- t(incoming)
    } else if (pattern == "all") {
      mat <- t(outgoing + incoming)
    }
    # if (!is.null(signaling)) {
    #   mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    #   mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    #   idx <- match(rownames(mat1), signaling)
    #   mat[idx[!is.na(idx)], ] <- mat1
    #   dimnames(mat) <- list(signaling, colnames(mat1))
    # }
    mat.ori <- mat
    mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
    # mat[mat == 0] <- NA
    sham <- mat
  }
  if(TRUE){
    object <- object.list[[2]]
    pattern = c("outgoing", "incoming", "all")[2]
    slot.name = "netP"
    
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
      outgoing[, i] <- centr[[i]]$outdeg
      incoming[, i] <- centr[[i]]$indeg
    }
    if (pattern == "outgoing") {
      mat <- t(outgoing)
    } else if (pattern == "incoming") {
      mat <- t(incoming)
    } else if (pattern == "all") {
      mat <- t(outgoing + incoming)
    }
    # if (!is.null(signaling)) {
    #   mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    #   mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    #   idx <- match(rownames(mat1), signaling)
    #   mat[idx[!is.na(idx)], ] <- mat1
    #   dimnames(mat) <- list(signaling, colnames(mat1))
    # }
    mat.ori <- mat
    mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
    # mat[mat == 0] <- NA
    choline <- mat
  }  
  
  all(rownames(sham) == rownames(choline))
  
  diff.in <- sham %>% as.data.frame %>% dplyr::select(`Endothelial`) %>%
    tibble::rownames_to_column("pathway") %>% 
    full_join(choline %>% as.data.frame %>% dplyr::select(`Endothelial`) %>% 
                tibble::rownames_to_column("pathway"), by = "pathway") 
  diff.in[is.na(diff.in)] <- 0
  diff.in <- diff.in %>% 
    mutate(diff.in = `Endothelial.y` - `Endothelial.x`) %>% 
    arrange(desc(diff.in))
  
  write.csv(diff.in, "2.CellChat/compare/Mesenchyme差异通路(incoming).csv")
  
  diffpw.in <- diff.in$pathway[abs(diff.in$diff) > 0.1]
}
if(outcoming){
  if(TRUE){
    object <- object.list[[1]]
    pattern = c("outgoing", "incoming", "all")[1]
    slot.name = "netP"
    
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
      outgoing[, i] <- centr[[i]]$outdeg
      incoming[, i] <- centr[[i]]$indeg
    }
    if (pattern == "outgoing") {
      mat <- t(outgoing)
    } else if (pattern == "incoming") {
      mat <- t(incoming)
    } else if (pattern == "all") {
      mat <- t(outgoing + incoming)
    }
    # if (!is.null(signaling)) {
    #   mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    #   mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    #   idx <- match(rownames(mat1), signaling)
    #   mat[idx[!is.na(idx)], ] <- mat1
    #   dimnames(mat) <- list(signaling, colnames(mat1))
    # }
    mat.ori <- mat
    mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
    # mat[mat == 0] <- NA
    sham <- mat
  }
  if(TRUE){
    object <- object.list[[2]]
    pattern = c("outgoing", "incoming", "all")[1]
    slot.name = "netP"
    
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
      outgoing[, i] <- centr[[i]]$outdeg
      incoming[, i] <- centr[[i]]$indeg
    }
    if (pattern == "outgoing") {
      mat <- t(outgoing)
    } else if (pattern == "incoming") {
      mat <- t(incoming)
    } else if (pattern == "all") {
      mat <- t(outgoing + incoming)
    }
    # if (!is.null(signaling)) {
    #   mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    #   mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    #   idx <- match(rownames(mat1), signaling)
    #   mat[idx[!is.na(idx)], ] <- mat1
    #   dimnames(mat) <- list(signaling, colnames(mat1))
    # }
    mat.ori <- mat
    mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
    # mat[mat == 0] <- NA
    choline <- mat
  }  
  
  all(rownames(sham) == rownames(choline))
  
  diff.out <- sham %>% as.data.frame %>% dplyr::select(`Endothelial`) %>%
    tibble::rownames_to_column("pathway") %>% 
    full_join(choline %>% as.data.frame %>% dplyr::select(`Endothelial`) %>% 
                tibble::rownames_to_column("pathway"), by = "pathway") 
  diff.out[is.na(diff.out)] <- 0
  diff.out <- diff.out %>% 
    mutate(diff.out = `Endothelial.y` - `Endothelial.x`) %>% 
    arrange(desc(diff.out))
  
  write.csv(diff.out, "2.CellChat/compare/Mesenchyme差异通路(out).csv")
  
  diffpw.out <- diff.out$pathway[abs(diff.out$diff) > 0.1]
}
if(all){
  if(TRUE){
    object <- object.list[[1]]
    pattern = c("outgoing", "incoming", "all")[3]
    slot.name = "netP"
    
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
      outgoing[, i] <- centr[[i]]$outdeg
      incoming[, i] <- centr[[i]]$indeg
    }
    if (pattern == "outgoing") {
      mat <- t(outgoing)
    } else if (pattern == "incoming") {
      mat <- t(incoming)
    } else if (pattern == "all") {
      mat <- t(outgoing + incoming)
    }
    # if (!is.null(signaling)) {
    #   mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    #   mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    #   idx <- match(rownames(mat1), signaling)
    #   mat[idx[!is.na(idx)], ] <- mat1
    #   dimnames(mat) <- list(signaling, colnames(mat1))
    # }
    mat.ori <- mat
    mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
    # mat[mat == 0] <- NA
    sham <- mat
  }
  if(TRUE){
    object <- object.list[[2]]
    pattern = c("outgoing", "incoming", "all")[3]
    slot.name = "netP"
    
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
      outgoing[, i] <- centr[[i]]$outdeg
      incoming[, i] <- centr[[i]]$indeg
    }
    if (pattern == "outgoing") {
      mat <- t(outgoing)
    } else if (pattern == "incoming") {
      mat <- t(incoming)
    } else if (pattern == "all") {
      mat <- t(outgoing + incoming)
    }
    # if (!is.null(signaling)) {
    #   mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    #   mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    #   idx <- match(rownames(mat1), signaling)
    #   mat[idx[!is.na(idx)], ] <- mat1
    #   dimnames(mat) <- list(signaling, colnames(mat1))
    # }
    mat.ori <- mat
    mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
    # mat[mat == 0] <- NA
    choline <- mat
  }  
  
  all(rownames(sham) == rownames(choline))
  
  diff <- sham %>% as.data.frame %>% dplyr::select(`Endothelial`) %>%
    tibble::rownames_to_column("pathway") %>% 
    full_join(choline %>% as.data.frame %>% dplyr::select(`Endothelial`) %>% 
                tibble::rownames_to_column("pathway"), by = "pathway") 
  diff[is.na(diff)] <- 0
  diff <- diff %>% 
    mutate(diff = `Endothelial.y` - `Endothelial.x`) %>% 
    arrange(desc(diff))
  
  write.csv(diff, "2.CellChat/compare/Mesenchyme差异通路(all).csv")
  
  diffpw.all <- diff$pathway[abs(diff$diff) > 0.1]
}

ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "all", signaling = pathway.union, title = names(object.list)[1], width = 5, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "all", signaling = pathway.union, title = names(object.list)[2], width = 5, height = 20)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# ------------------------------------------------------------------------------
ht1 = netAnalysis_signalingRole_heatmap(
  object.list[[1]], pattern = "incoming", signaling = diffpw.in, 
  title = names(object.list)[1], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(
  object.list[[2]], pattern = "incoming", signaling = diffpw.in, 
  title = names(object.list)[2], width = 5, height = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(
  object.list[[1]], pattern = "outgoing", signaling = diffpw.out, 
  title = names(object.list)[1], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(
  object.list[[2]], pattern = "outgoing", signaling = diffpw.out, 
  title = names(object.list)[2], width = 5, height = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(
  object.list[[1]], pattern = "all", signaling = diffpw.all, 
  title = names(object.list)[1], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(
  object.list[[2]], pattern = "all", signaling = diffpw.all, 
  title = names(object.list)[2], width = 5, height = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# Part III: Identify the up-gulated and down-regulated signaling ligand-receptor pairs====
# Identify dysfunctional signaling by comparing the communication probabities
# CellChat can compare the communication probabilities mediated by L-R pairs from certain cell groups to other cell groups. This can be done by setting comparison in the function netVisual_bubble.
pathways <- c(
  unique(subsetCommunication(cellchat, slot.name = "netP")$nondrink$pathway_name),
  unique(subsetCommunication(cellchat, slot.name = "netP")$drink$pathway_name)) %>% 
  unique
for(i in downp){
  p <- netVisual_bubble(cellchat, 
                        sources.use = levels(cc.chol@idents), 
                        targets.use = levels(cc.chol@idents),  
                        comparison = c(1, 2), 
                        signaling = i,
                        remove.isolate = TRUE,
                        angle.x = 45)
  ggsave(p, filename = paste0("bubble_pathway-",i,".pdf"), width=15, height=8, units="in")
}

p <- netVisual_bubble(cellchat, 
                      sources.use = c(1:10), 
                      targets.use = 7,  
                      comparison = c(1, 2), 
                      remove.isolate = TRUE,
                      angle.x = 45)
write.csv(p$data, "所有靶向Treg的通路.csv")
# Identify dysfunctional signaling by using differential expression analysis
# The above method for identifying the upgulated and down-regulated signaling is perfomed by comparing the communication probability between two datasets for each L-R pair and each pair of cell groups. Alternative, we can identify the upgulated and down-regulated signaling ligand-receptor pairs based on the differential expression analysis (DEA). Specifically, we perform differential expression analysis between two biological conditions (i.e., NL and LS) for each cell group, and then obtain the upgulated and down-regulated signaling based on the fold change of ligands in the sender cells and receptors in the receiver cells.
# Of note, users may observe the same LR pairs appearing in both the up-regulated and down-regulated results due to the fact that DEA between conditions is performed for each cell group. To perform DEA between conditions by ignoring cell group information, users can set group.DE.combined = TRUE in identifyOverExpressedGenes for CellChat v2.1.1.
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "drink"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 
cellchat <- identifyOverExpressedGenes(
  cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, 
  features.name = features.name, only.pos = FALSE, thresh.pc = 0.2, 
  thresh.fc = log2(1.2), thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, ligand.logFC = NULL, receptor.logFC = log2(1.5), receptor.pvalues = 0.05)
# net.up <- subsetCommunication(cellchat, net = net, datasets = "LS",ligand.logFC = 0.05, receptor.logFC = NULL)
# net.up <- subsetCommunication(cellchat, net = net, datasets = "NL",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, ligand.logFC = NULL, receptor.logFC = -log2(1.5), receptor.pvalues = 0.05)
# net.down <- subsetCommunication(cellchat, net = net, datasets = "NL",ligand.logFC = -0.05, receptor.logFC = NULL)
# Since the signaling genes in the net.up and net.down might be complex with multi-subunits, we can do further deconvolution to obtain the individual signaling genes.

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
# Users can also find all the significant outgoing/incoming/both signaling according to the customized features and cell groups of interest

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1:10), 
                        targets.use = "CD4 Treg_FOXP3", comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,title.name = "upregulated LR")
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1:10), 
                        targets.use = "CD4 Treg_FOXP3", comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,title.name = "downregulated LR")
#> Comparing communications on a merged object
gg1 + gg2

df <- findEnrichedSignaling(object.list[[2]], features = c("CCL19", "CXCL12"), idents = c("Inflam. FIB", "COL11A1+ FIB"), pattern ="outgoing")
# Visualize the identified up-regulated and down-regulated signaling ligand-receptor pairs
# CellChat can visualize the identified up-regulated and down-regulated signaling ligand-receptor pairs using bubble plot (option A), chord diagram (option B) or wordcloud (option C).
# (C) Wordcloud plot
# Visualize the enriched ligands, signaling,or ligand-receptor pairs in one condition compared to another condition using wordcloud
# visualize the enriched ligands in the first condition
computeEnrichmentScore(net.down, species = 'mouse', variable.both = TRUE)
# visualize the enriched ligands in the second condition
computeEnrichmentScore(net.up, species = 'mouse', variable.both = TRUE)

# Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram=====
# Similar to the CellChat analysis of individual dataset, CellChat can visually compare cell-cell communication networks using hierarchy plot, circle plot, chord diagram, or heatmap. More details on the visualization can be found in the CellChat analysis of individual dataset.

# Edge color/weight, node color/size/shape: 
# In all visualization plots, edge colors are consistent with the sources as sender, 
# and edge weights are proportional to the interaction strength. 
# Thicker edge line indicates a stronger signal. In the Hierarchy plot and Circle plot, 
# circle sizes are proportional to the number of cells in each cell group. 
# In the hierarchy plot, solid and open circles represent source and target, respectively. 
# In the Chord diagram, the inner thinner bar colors represent the targets that receive signal from the corresponding outer bar. The inner bar size is proportional to the signal strength received by the targets. Such inner bar is helpful for interpreting the complex chord diagram. Note that there exist some inner bars without any chord for some cell groups, please just igore it because this is an issue that has not been addressed by circlize package.
pathways <- intersect(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
pdf("两组间网络.pdf")
for (i in pathways) {
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = i) # control the edge weights across different datasets
  par(mfrow = c(1,2), xpd=TRUE)
  for (j in 1:length(object.list)) {
    netVisual_aggregate(object.list[[j]], signaling = i, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(i, names(object.list)[j]))
  }
}
dev.off()


pathways.show <- c("IL4", "CD48", "CDH1","SIRP")[2] 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}
# 如果遇到
# 错误于slot(x, slot.name[i])$prob[, , attribute[i]]: 下标出界
# 原因是一个组中没有改信号通路的信号，请根据实际情况，分组绘图

par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Chord diagram
pathways.show <- c("COLLAGEN", "LAMININ", "FN1")[1] 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

# For the chord diagram, CellChat has an independent function netVisual_chord_cell to flexibly visualize the signaling network by adjusting different parameters in the circlize package. For example, we can define a named char vector group to create multiple-group chord diagram, e.g., grouping cell clusters into different cell types.
# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(object.list[[1]]@idents)
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Plot the aggregated cell-cell communication network at the signaling pathway level
# Using chord diagram, CellChat provides two functions netVisual_chord_cell and netVisual_chord_gene for visualizing cell-cell communication with different purposes and different levels. netVisual_chord_cell is used for visualizing the cell-cell communication between different cell groups (where each sector in the chord diagram is a cell group), and netVisual_chord_gene is used for visualizing the cell-cell communication mediated by mutiple ligand-receptors or signaling pathways (where each sector in the chord diagram is a ligand, receptor or signaling pathway.)

par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 4, targets.use = c(5:8), lab.cex = 0.5, title.name = paste0("Signaling from Inflam.FIB - ", names(object.list)[i]))
}

# compare all the interactions sending from fibroblast to inflamatory immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1,2, 3, 4), targets.use = c(8,10),  title.name = paste0("Signaling received by Inflam.DC and .TC - ", names(object.list)[i]), legend.pos.x = 10)
}
# show all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1:10), targets.use = 7,
                       big.gap = 3, small.gap = 0.25,
                       slot.name = "netP", title.name = paste0("Signaling pathways sending to Treg - ", names(object.list)[i]), 
                       legend.pos.x = 10)
}
# NB: Please ignore the note when generating the plot such as “Note: The first link end is drawn out of sector ‘MIF’.”. If the gene names are overlapped, you can adjust the argument small.gap by decreasing the value.

# Part V: Compare the signaling gene expression distribution between different datasets====
# We can plot the gene expression distribution of signaling genes related to L-R pairs or signaling pathway using a Seurat wrapper function plotGeneExpression.

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("nondrink", "drink")) # set factor level
p <- plotGeneExpression(cellchat, 
                        signaling = "CD34", 
                        split.by = "datasets", 
                        colors.ggplot = T,
                        color.use = c("#548235","#C00000"),
                        type = "violin")

for (i in c("IL4","VISTA","CD34","SPP1",#"GH",
            "SELL","CDH1","SIRP","CD70","CD48")) {
  p <- plotGeneExpression(cellchat, 
                          signaling = i, 
                          split.by = "datasets", 
                          colors.ggplot = T,
                          color.use = c("#548235","#C00000"),
                          type = "violin")
  ggsave(p, filename = paste0(i, ".pdf"), width=8, height=8, units="in")
}
VlnPlot(scobj.Mesen.h, group.by = "drinking", features = c("CD48","CD244"))
VlnPlot(scobj.Mesen.h, group.by = "drinking", features = c("GH1","GHR"))
VlnPlot(scobj.Mesen.h, group.by = "drinking", features = c("GP1BA","ITGAM","ITGB2"))
# https://github.com/sqjin/CellChat/issues/317
# Cannot find  GH . Please input a correct name! 


mypathway <- "IL4"
a <- cellchat@LR$drink$LRsig[
  which(cellchat@LR$drink$LRsig$pathway_name == mypathway),]
b <- cellchat@LR$nondrink$LRsig[
  which(cellchat@LR$nondrink$LRsig$pathway_name == mypathway),]
unique(
  c(unique(c(a$ligand, b$ligand)),
    unique(c(a$receptor, b$receptor)))
)

#===============================================================================
# SECTION III - Cell Cell Communication=========================================
#===============================================================================
fig <- subset(scobj.harmony, subset = cell_type != "doublet")
fig$sample <- fig$orig.ident
fig$RNA_snn_res.0.2 <- factor(fig$RNA_snn_res.0.4, levels = c(0:12))
fig$cell_type <- factor(fig$cell_type, 
                        levels = c("Cardiomyocyte", "Endothelial", "Mesenchyme", 
                                   "Lymphoid", "Myeloid", "Neureonal"))
fig$orig.ident <- factor(fig$orig.ident, levels = c("sham", "choline"))

table(fig@meta.data$cell_type)
# Cardiomyocyte   Endothelial    Mesenchyme      Lymphoid       Myeloid     Neureonal 
# 1097         12246          4660           401           858            63 
#### 1.Create cellchat obj------------------------------------------------------
library(CellChat)
cellChat <- createCellChat(object = fig, group.by = "cell_type", assay = "RNA")

cellChat
# An object of class CellChat created from a single dataset 
# 21357 genes.
# 19325 cells. 
# CellChat analysis of single cell RNA-seq data! 
#### 2.Load LR DB---------------------------------------------------------------
# # use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- subsetDB(CellChatDB)
cellChat@DB <- CellChatDB.use
#### 3.Preparement--------------------------------------------------------------
cellChat <- subsetData(cellChat) # This step is necessary even if using the whole database

cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)

cellChat <- smoothData(cellChat, adj = PPI.mouse) # 2.1.2 adj = is forced
# Part II: Inference of cell-cell communication network-------------------------
cellChat <- computeCommunProb(cellChat, 
                              raw.use = FALSE, 
                              population.size = TRUE,
                              type =  "truncatedMean", 
                              trim = 0.2) 

cellChat <- filterCommunication(cellChat, min.cells = 10)

cellChat <- computeCommunProbPathway(cellChat)

# save(cellChat, file = "Result/scRNAseq/method3/2. Cellchat/bulk cellChat.Rdata")
load("Result/scRNAseq/method3/2. Cellchat/bulk cellChat.Rdata")

# subsetCommunication(cellChat, slot.name = "net") %>% write.csv("2.CellChat/bulk/bulk cellchat-LR.csv")
# subsetCommunication(cellChat, slot.name = "netP")%>% write.csv("2.CellChat/bulk/bulk cellchat-pathway.csv")
#### 5.Visualization------------------------------------------------------------
cellChat <- aggregateNet(cellChat)

# write.csv(cellChat@net$count %>% as.data.frame, "2.CellChat/bulk/bulk cellchat count.csv")
# write.csv(cellChat@net$weight %>% as.data.frame, "2.CellChat/bulk/bulk cellchat weight.csv")

pdf("2.CellChat/bulk cellchat 1.pdf", width = 8, height = 8)
pheatmap::pheatmap(cellChat@net$count, 
                   color = colorRampPalette(c("white","#007947"))(200),
                   cluster_rows = F, cluster_cols = F)
pheatmap::pheatmap(cellChat@net$weight, 
                   color = colorRampPalette(c("white","#007947"))(200),
                   cluster_rows = F, cluster_cols = F)
dev.off()

# CellChat can also visualize the aggregated cell-cell communication network. For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellChat@idents))
par(mfrow = c(1,2), xpd=TRUE)
# netVisual_circle(cellChat@net$count, vertex.weight = 1, weight.scale = T, label.edge= F, title.name = "Number of interactions")
# netVisual_circle(cellChat@net$weight, vertex.weight = 1, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_circle(cellChat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Visualize the inferred cell-cell communication network
# ..vertex.weight	The weight of vertex: either a scale value or a vector
#   Default is a scale value being 1, indicating all vertex is plotted in the same 
#   size;
#   Set 'vertex.weight' as a vector to plot vertex in different size; setting 
#   'vertex.weight = NULL' will have vertex with different size that are portional 
#   to the number of cells in each cell group.


# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. 
# Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.

mat <- cellChat@net$count
par(mfrow = c(2,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

mat <- cellChat@net$weight
par(mfrow = c(2,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
#===============================================================================
# SECTION IV - Mesenchyme =======================================================
#===============================================================================
load("Result/scRNAseq/method4/scobj.harmony(singlet-withoutreg.cc-2000-25) anno.Rdata")

scobj.harmony@active.assay <- "RNA"
scobj.Mesen <- scobj.harmony %>% subset(subset = cell_type %in% c("Mesenchyme"))

dim(scobj.Mesen)
# [1] 21320  4829

scobj.Mesen <- CreateSeuratObject(
  counts = scobj.Mesen@assays$RNA$counts, 
  assay = "RNA", 
  meta.data = scobj.Mesen@meta.data, 
  project = "Heart-mesenchyme cell")
# QC----------------------------------------------------------------------------
table(scobj.Mesen@meta.data$orig.ident)
# choline    sham 
# 2266    2563
# Nomalization & harmony & cluster----------------------------------------------
# 1000 without regress cell cycle
if(FALSE){
  scobj.Mesen.1000 <- scobj.Mesen %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony")
  
  scobj.Mesen.1000.10 <- scobj.Mesen.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:10) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Mesen.1000.15 <- scobj.Mesen.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:15) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Mesen.1000.20 <- scobj.Mesen.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Mesen.1000.25 <- scobj.Mesen.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Mesen.1000.30 <- scobj.Mesen.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Mesenchyme vst.1000.h_clustree 10 15 20 25 30(singlet) without regression.pdf")
  ElbowPlot(scobj.Mesen.1000, reduction = "pca", ndims = 50)
  clustree(scobj.Mesen.1000.10, prefix = "RNA_snn_res.")
  clustree(scobj.Mesen.1000.15, prefix = "RNA_snn_res.")
  clustree(scobj.Mesen.1000.20, prefix = "RNA_snn_res.")
  clustree(scobj.Mesen.1000.25, prefix = "RNA_snn_res.")
  clustree(scobj.Mesen.1000.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.Mesen.1000, scobj.Mesen.1000.10, scobj.Mesen.1000.15,
       scobj.Mesen.1000.20, scobj.Mesen.1000.25, scobj.Mesen.1000.30,
       file = "scobj.Mesen.1000.Rdata")
}
# 1500 without regress cell cycle
if(TRUE){
  scobj.Mesen.1500 <- scobj.Mesen %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony")
  
  scobj.Mesen.1500.10 <- scobj.Mesen.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:10) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Mesen.1500.15 <- scobj.Mesen.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:15) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Mesen.1500.20 <- scobj.Mesen.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Mesen.1500.25 <- scobj.Mesen.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Mesen.1500.30 <- scobj.Mesen.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Mesenchyme vst.1500.h_clustree 10 15 20 25 30(singlet) without regression.pdf")
  ElbowPlot(scobj.Mesen.1500, reduction = "pca", ndims = 50)
  clustree(scobj.Mesen.1500.10, prefix = "RNA_snn_res.")
  clustree(scobj.Mesen.1500.15, prefix = "RNA_snn_res.")
  clustree(scobj.Mesen.1500.20, prefix = "RNA_snn_res.")
  clustree(scobj.Mesen.1500.25, prefix = "RNA_snn_res.")
  clustree(scobj.Mesen.1500.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.Mesen.1500, scobj.Mesen.1500.10, scobj.Mesen.1500.15,
       scobj.Mesen.1500.20, scobj.Mesen.1500.25, scobj.Mesen.1500.30,
       file = "scobj.Mesen.1500.Rdata")
}
# 2000 without regress cell cycle
if(FALSE){
  scobj.Mesen.2000 <- scobj.Mesen %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony")
  
  scobj.Mesen.2000.10 <- scobj.Mesen.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:10) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Mesen.2000.15 <- scobj.Mesen.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:15) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Mesen.2000.20 <- scobj.Mesen.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Mesen.2000.25 <- scobj.Mesen.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Mesen.2000.30 <- scobj.Mesen.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Mesenchyme vst.2000.h_clustree 10 15 20 25 30(singlet) without regression.pdf")
  ElbowPlot(scobj.Mesen.2000, reduction = "pca", ndims = 50)
  clustree(scobj.Mesen.2000.10, prefix = "RNA_snn_res.")
  clustree(scobj.Mesen.2000.15, prefix = "RNA_snn_res.")
  clustree(scobj.Mesen.2000.20, prefix = "RNA_snn_res.")
  clustree(scobj.Mesen.2000.25, prefix = "RNA_snn_res.")
  clustree(scobj.Mesen.2000.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.Mesen.2000, scobj.Mesen.2000.10, scobj.Mesen.2000.15,
       scobj.Mesen.2000.20, scobj.Mesen.2000.25, scobj.Mesen.2000.30,
       file = "scobj.Mesen.2000.Rdata")
}

# 1000 30 0.6 0-8
# 1500 30 0.8 0-10
# 2000 10 0.7 0-10
# 2000 20 0.7 0-9

scobj.Mesen.h <- scobj.Mesen.1500.30 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)
# Visualization-----------------------------------------------------------------
DimPlot(scobj.Mesen.h, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.Mesen.h, group.by = "RNA_snn_res.0.1", reduction = "umap", label = T)+
DimPlot(scobj.Mesen.h, group.by = "RNA_snn_res.0.5", reduction = "umap", label = T)+
DimPlot(scobj.Mesen.h, group.by = "RNA_snn_res.0.7", reduction = "umap", label = T)+
DimPlot(scobj.Mesen.h, group.by = "RNA_snn_res.0.8", reduction = "umap", label = T)
# Annotation====================================================================
## Findmarkers------------------------------------------------------------------
Clist <- list()
n <- 1
for (i in 0:10) {
  Clist[[n]] <- FindMarkers(scobj.Mesen.h, group.by = "RNA_snn_res.0.8", ident.1 = as.character(i))
  write.csv(Clist[[n]], paste0("D:/C",n-1,".csv"))
  
  n <- n + 1
}
## Annotation-------------------------------------------------------------------
scobj.Mesen.h$cell_type <- "Unknown"
scobj.Mesen.h$cell_type[scobj.Mesen.h$RNA_snn_res.0.8 %in% c(0)]  <- "FB_Ogn"
scobj.Mesen.h$cell_type[scobj.Mesen.h$RNA_snn_res.0.8 %in% c(1)]  <- "FB_CD74"
scobj.Mesen.h$cell_type[scobj.Mesen.h$RNA_snn_res.0.8 %in% c(2)]  <- "FB_Egfr"
scobj.Mesen.h$cell_type[scobj.Mesen.h$RNA_snn_res.0.8 %in% c(3)]  <- "SMC_Acta2"
scobj.Mesen.h$cell_type[scobj.Mesen.h$RNA_snn_res.0.8 %in% c(4)]  <- "SMC_Rgs5"
scobj.Mesen.h$cell_type[scobj.Mesen.h$RNA_snn_res.0.8 %in% c(5)]  <- "Pericyte_Abcc9"
scobj.Mesen.h$cell_type[scobj.Mesen.h$RNA_snn_res.0.8 %in% c(6)]  <- "SMC_Cdh5"
scobj.Mesen.h$cell_type[scobj.Mesen.h$RNA_snn_res.0.8 %in% c(7)]  <- "FB_Mmp2"
scobj.Mesen.h$cell_type[scobj.Mesen.h$RNA_snn_res.0.8 %in% c(8)]  <- "FB_Postn"
scobj.Mesen.h$cell_type[scobj.Mesen.h$RNA_snn_res.0.8 %in% c(9)]  <- "FB_Mfap5"
scobj.Mesen.h$cell_type[scobj.Mesen.h$RNA_snn_res.0.8 %in% c(10)] <- "FB_Igf1"

table(scobj.Mesen.h$RNA_snn_res.0.8)
# 0    1    2    3    4    5    6    7    8    9   10 
# 1120  746  556  416  410  405  342  332  247  192   63 
table(scobj.Mesen.h$cell_type)
# FB_CD74        FB_Egfr        FB_Igf1       FB_Mfap5        FB_Mmp2 
# 746            556             63            192            332 
# FB_Ogn       FB_Postn Pericyte_Abcc9      SMC_Acta2       SMC_Cdh5 
# 1120            247            405            416            342 
# SMC_Rgs5 
# 410

DimPlot(scobj.Mesen.h, group.by = "cell_type", reduction = "umap", label = T)

FeaturePlot(scobj.Mesen.h, features = "nCount_RNA")
VlnPlot(scobj.Mesen.h, features = "nCount_RNA", group.by = "cell_type")
FeaturePlot(scobj.Mesen.h, features = "nFeature_RNA")
VlnPlot(scobj.Mesen.h, features = "nFeature_RNA", group.by = "cell_type")

tmp <- scobj.Mesen.h@meta.data %>% 
  filter(cell_type == 'FB_CD74')
summary(tmp$nCount_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 300.0   448.0   606.5   732.4   882.5  2609.0 
summary(tmp$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 263.0   393.0   514.5   598.6   723.0  1799.0 

save(scobj.Mesen.h, file = "Result/scRNAseq/method4/2. Mesenchyme/scobj.Mesen.h anno.Rdata")
# Function Annotation===========================================================
# NF_CFD------------------------------------------------------------------------
C3 <- FindMarkers(scobj.Mesen.h, group.by = "RNA_snn_res.0.4", ident.1 = "3")

geneList <- bitr(rownames(C3)[C3$p_val_adj < 0.05 & C3$pct.1 > 0.2 & abs(C3$avg_log2FC) >1], 
                 fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
ORA_KEGG <- enrichKEGG(
  gene = geneList$ENTREZID,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 1, 
  qvalueCutoff = 1)
ORA_KEGG <- setReadable(ORA_KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
p <- dotplot(ORA_KEGG, showCategory = 20) + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p, filename = "NF_CFD KEGG.pdf", width=7, height=6, units="in")
write.csv(ORA_KEGG@result, "NF_CFD KEGG.csv")

ORA_Reactome <- ReactomePA::enrichPathway(
  gene= geneList$ENTREZID, 
  pvalueCutoff = 1, 
  readable=TRUE, 
  organism = "human")
p <- clusterProfiler::dotplot(ORA_Reactome)
ggsave(p, filename = "NF_CFD REACTOME.pdf", width=7, height=6, units="in")
write.csv(ORA_Reactome@result, "NF_CFD REACTOME.csv")
save(ORA_Reactome, file = "NF_CFD REACTOME.Rdata")
## Hallmark=====================================================================
require(GSVA)

load("4.Mesenchyme/scobj.Mesen.h anno removing dbl.Rdata")
counts <- subset(scobj.Mesen.h, subset = cell_type %in% c("NF_CFD")) %>% 
  AggregateExpression(group.by = "orig.ident")
counts <- counts$RNA %>% as.data.frame
dim(counts)
# [1] 25414    12

hallmarks <- qusage::read.gmt("C:/D/Project/63.ESCC/h.all.v2023.1.Hs.symbols.gmt") #返回的是表格
# gsva_matrix <- GSVA::gsva(as.matrix(count), hallmarks, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
# https://blog.csdn.net/dolanjiang/article/details/144673628

# https://zhuanlan.zhihu.com/p/518145829
# version <= 1.50.5
# 使用GSVA需要输入基因表达矩阵和基因集。 基因集即为我们上一步所得list；
# 基因表达矩阵可以使用logCPM、logRPKM、logTPM（GSVA参数kcdf选择"Gaussian"，默认）或
# counts数据（参数kcdf选择"Poisson"）。 GSVA还支持BiocParallel，可设置参数
# parallel.sz进行多核计算。

gsvaPar <- ssgseaParam(exprData = as.matrix(counts), 
                       geneSets = hallmarks,
                       normalize = TRUE)
gsva_data <- gsva(gsvaPar, verbose = FALSE)

annotation_cols <- data.frame(
  Group = c(1,1,1,0,1,0,0,1,1,0,1,1) %>%
    factor(levels = c(0,1), labels = c("noncholine","choline")),
  row.names = colnames(gsva_data)
)
pheatmap::pheatmap(gsva_data, 
                   clustering_method = "ward.D",
                   annotation_col = annotation_cols,
                   annotation_colors = list(
                     Group = c(noncholine = "black", choline = "#C00000")),
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   scale = "row")
# Figure2=======================================================================
fig <- subset(scobj.Mesen.h, subset = cell_type != "doublets")
fig$sample <- fig$orig.ident
fig$RNA_snn_res.0.8 <- factor(fig$RNA_snn_res.0.4, levels = c(0:10))
fig$orig.ident <- factor(fig$orig.ident, levels = c("sham", "choline"))
table(fig$cell_type)
# FB_CD74        FB_Egfr        FB_Igf1       FB_Mfap5        FB_Mmp2 
# 746            556             63            192            332 
# FB_Ogn       FB_Postn Pericyte_Abcc9      SMC_Acta2       SMC_Cdh5 
# 1120            247            405            416            342 
# SMC_Rgs5 
# 410

patient_col <- c("#c79494", "#3d7d35")

p1 <- DimPlot(fig, reduction = "umap", group.by = "cell_type", cols = colors_list,
              label = TRUE, pt.size = 0.1, raster = FALSE)
p2 <- DimPlot(fig, reduction = "umap",
              group.by = "orig.ident", cols = patient_col, 
              label = TRUE, pt.size = 0.1, raster = FALSE)
p <- p1 + p2 + plot_layout(ncol = 2, nrow = 1)
ggsave(p, filename = "Mesenchyme cell UMAP.pdf", path = "Figure/", width=12, height=5, units="in")

p <- DotPlot(fig, features = c(
  "Cd74", "H2-Aa", "H2-Eb1", "H2-Ab1", # "Cxcl9", "Cxcl12",  "Smn1", "Pdgfra",
  # "Cthrc1",  "Aifm2",
  'Dcn', # 'C3', 'C7',
  "Col1a1", "Col1a2","Col1a3", "Fbln1", 
  "Egfr", 'Igf1', 'Mfap5', "Mmp2", "Ogn", 
  "Vcan", "Col8a1","Postn",
  "Rsg5", "Abcc9", "Kcnj8",  'Pdgfrb', 'Cspg4',
  'Myl9','Myh11', 'Acta2','Cnn1', 'Tagln', 'Myocd', 'Cald1', 'Mylk', 
  "Cdh5","Pecam1","Eng","Tek"
) %>% unique,
group.by = "cell_type") +
  scale_color_viridis() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, filename = "Mesenchyme MARKER.pdf", width=10, height=10, units="in")

dat <- fig@meta.data %>% 
  dplyr::select(orig.ident, cell_type)
p <- dat %>% 
  mutate(value = 1) %>% 
  group_by(orig.ident, cell_type) %>%
  summarise(n = n()) %>% 
  group_by(orig.ident) %>%
  mutate(sum = sum(n)) %>% 
  mutate(n = n / sum) %>% 
  ggplot(aes(x = orig.ident, y = n, fill = cell_type, stratum = cell_type, alluvium = cell_type))+
  geom_stratum(width = 0.5, color='white') +
  geom_alluvium(alpha = 0.5,
                width = 0.5,
                curve_type = "linear") +
  scale_color_manual(values = colors_list) +
  scale_fill_manual(values = colors_list) +
  labs(x = "", y = "Percent", fill = "Cell type") +
  theme_bw()+
  theme(text = element_text(size = 20))
ggsave(p, filename = "Mesenchyme percent.pdf", width=6, height=6, units="in")
# pseudobulk--------------------------------------------------------------------
geom_volcano <- function(dat,pos.num = 10,neg.num = -10){
  require(ggplot2)
  require(ggrepel)
  
  dat <- dat %>% 
    tibble::rownames_to_column("SYMBOL") %>% 
    mutate(color = case_when(log2FoldChange >  1 & padj < 0.05 ~ 2,
                             log2FoldChange < -1 & padj < 0.05 ~ 1,
                             TRUE ~ 0)) %>% 
    mutate(label = factor(color, levels = c(0,1,2), labels = c("Non-Sig","Down","Up")),
           color = factor(color))
  dat_up <- dat %>% 
    filter(label == "Up") %>% 
    arrange(desc(log2FoldChange))
  dat_up <- dat_up[1:ifelse(nrow(dat_up) >= 10, 10, nrow(dat_up)),]
  if(is.na(dat_up$SYMBOL[1]) & nrow(dat_up) == 1){
    dat_up[1,] <- c(NA,0,0,0,0,1,1,NA,NA)
    dat_up$log2FoldChange <- dat_up$log2FoldChange %>% as.numeric
    dat_up$padj <- dat_up$padj %>% as.numeric
  }
  dat_down <- dat %>% 
    filter(label == "Down") %>% 
    arrange(log2FoldChange)
  dat_down <- dat_down[1:ifelse(nrow(dat_down) >= 10, 10, nrow(dat_down)),]
  if(is.na(dat_down$SYMBOL[1]) & nrow(dat_down) == 1){
    dat_down[1,] <- c("",0,0,0,0,1,1,NA,NA)
    dat_down$log2FoldChange <- dat_down$log2FoldChange %>% as.numeric
    dat_down$padj <- dat_down$padj %>% as.numeric
  }
  
  ggplot(dat, aes(x = log2FoldChange, y = -log10(padj), col = log2FoldChange, label = SYMBOL)) +
    geom_point(
      # aes(size = !!rlang::sym(abundance))
    )+
    geom_vline(xintercept = c(-1,1), color = "gray80", linetype = 2) +
    geom_hline(yintercept = c(1.30103), color = "gray80", linetype = 2) +
    ylab(expression(-log[10]~(adj.~P~value))) +
    xlab("Log2(Fold Change)") +
    labs(color = "") +
    scale_color_gradient2(high = "red3", mid = "white", low = "blue3",
                          midpoint = 0, na.value = "grey80"
                          #  space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour"
    )+
    scale_size_continuous(range = c(0.1, 4)) +
    geom_text_repel(
      data = dat_up,
      color = "red3",
      size = 5,
      nudge_x = pos.num - as.numeric(dat_up$log2FoldChange),
      segment.size=0.3,
      segment.color="grey",
      direction="y",
      hjust= 0,
      max.overlaps = Inf) +
    geom_text_repel(
      data= dat_down,
      color="blue3",
      size=5,
      nudge_x = neg.num - as.numeric(dat_down$log2FoldChange),
      segment.size = 0.3,
      segment.color = "grey",
      direction="y",
      hjust= 1,
      max.overlaps = Inf) +
    # labs(size = expression("Abundance (log2)"),
    #      color = expression("Direction signed"),
    #      title = trait.names[i]) +
    theme_minimal() +
    theme(legend.position = "right",
          legend.title.align = 0, # left align
          legend.title = element_text(margin = margin(t = 15, unit = "pt")) # add more space on top of legend titles
          #legend.spacing.y = unit(1,"cm")
    ) +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(size=20),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18),
          aspect.ratio = 1/1.2, panel.grid.major = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    theme(plot.title = element_text(hjust = 0.5, face = "italic", colour="grey50", size=20))
}

# smoking
for (i in unique(fig$cell_type)) {
  tmp <- fig %>% subset(subset = cell_type == i)
  
  counts <- AggregateExpression(tmp, group.by = "orig.ident")
  counts <- counts$RNA %>% as.data.frame
  
  group <- ifelse(colnames(counts) %in% c("sc-01","sc-04","sc-06","sc-07","sc-10"), "nonsmoker", "smoker")
  condition<- factor(group, levels = c("nonsmoker","smoker"))
  coldata <- data.frame(row.names = colnames(counts), condition)
  dds <- DESeqDataSetFromMatrix(countData = counts %>% as.matrix,
                                colData = coldata,
                                design = ~condition)
  dds <- DESeq(dds)  
  DEG_smoker_vs_nonsmoker <- results(dds, cooksCutoff=FALSE, 
                                     name="condition_smoker_vs_nonsmoker", 
                                     independentFiltering = FALSE) %>% 
    as.data.frame %>% 
    arrange(padj)
  DEG_smoker_vs_nonsmoker <- na.omit(DEG_smoker_vs_nonsmoker)
  
  save(counts, file = paste0("4.Mesenchyme/smoking/pseudobulk_",i,".Rdata"))
  save(DEG_smoker_vs_nonsmoker, file = paste0("4.Mesenchyme/smoking/",i,"_DESeq2-DEGs.Rdata"))
  export::table2excel(DEG_smoker_vs_nonsmoker, add.rownames = TRUE, sheetName = i,
                      paste0("4.Mesenchyme/smoking/",i,"_DESeq2-DEGs.xlsx"))
  
  print(i)
}

# drinking
for (i in unique(fig$cell_type)) {
  tmp <- fig %>% subset(subset = cell_type == i)
  
  counts <- AggregateExpression(tmp, group.by = "orig.ident")
  counts <- counts$RNA %>% as.data.frame
  
  group <- ifelse(colnames(counts) %in% c("sc-04","sc-06","sc-07","sc-10"), "noncholine", "choline")
  condition<- factor(group, levels = c("noncholine","choline"))
  coldata <- data.frame(row.names = colnames(counts), condition)
  dds <- DESeqDataSetFromMatrix(countData = counts %>% as.matrix,
                                colData = coldata,
                                design = ~condition)
  dds <- DESeq(dds)  
  # vsd <- assay(vst(dds, blind = FALSE))
  DEG_choline_vs_noncholine <- results(dds, cooksCutoff=FALSE, 
                                       name="condition_choline_vs_noncholine", 
                                       independentFiltering = FALSE) %>% 
    as.data.frame %>% 
    arrange(padj)
  DEG_choline_vs_noncholine <- na.omit(DEG_choline_vs_noncholine)
  
  save(counts, file = paste0("4.Mesenchyme/drinking/pseudobulk_",i,".Rdata"))
  save(DEG_choline_vs_noncholine, file = paste0("4.Mesenchyme/drinking/",i,"_DESeq2-DEGs.Rdata"))
  export::table2excel(DEG_choline_vs_noncholine, add.rownames = TRUE, sheetName = i,
                      paste0("4.Mesenchyme/drinking/",i,"_DESeq2-DEGs.xlsx"))
  
  print(i)
}
# Hallmark: to show inflammation pheno of NF ===================================
require(GSVA)

data <- AggregateExpression(fig, group.by = "cell_type")
dim(data$RNA)
# [1] 25414    10

count <- data$RNA %>% as.matrix

hallmarks <- qusage::read.gmt("4.Mesenchyme/HP_GASTROINTESTINAL_INFLAMMATION.v2025.1.Hs.gmt") #返回的是表格

gsvaPar <- ssgseaParam(exprData = count, 
                       geneSets = hallmarks,
                       normalize = TRUE)
gsva_data <- gsva(gsvaPar, verbose = FALSE)

# annotation_cols <- data.frame(
#   Group = c(rep("MRP", 4), rep("non-MRP", 6)) %>% 
#     factor(levels = c("non-MRP","MRP")),
#   row.names = colnames(gsva_data)
# )
pheatmap::pheatmap(gsva_data, 
                   clustering_method = "ward.D",
                   annotation_col = annotation_cols,
                   # annotation_colors = list(
                   #   Group = c(MRP = "#C00000", `non-MRP` = "black")),
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   scale = "row")

# HP_GASTROINTESTINAL_INFLAMMATION
# NF-CFD                                3.937135
# CAF-COL1A1                            4.006968
# CAF-HLA-DRA                           4.873358
# CAF-HTRA3                             3.940639
# CAF-KRT8                              4.297854
# CAF-MMP1                              4.050796
# CAF-MYH11                             3.909218
# Fibro-MKI67                           4.036532
# Pericyte-RGS5                         4.286126
# SMC-MYH11                             3.873358
# Roe---------------------------------------------------------------------------
data <- fig@meta.data
data$majorCluster = data$cell_type
data$patient = data$sample
data$loc = data$orig.ident
Roe <- calTissueDist(data,
                     byPatient = F,
                     colname.cluster = "majorCluster", # 不同细胞亚群
                     colname.patient = "patient", # 不同样本
                     colname.tissue = "loc", # 不同组织
                     method = "chisq", # "chisq", "fisher", and "freq" 
                     min.rowSum = 0) 
# sham   choline
# FB_CD74        0.6869714 1.3540566  *
# FB_Egfr        0.9420601 1.0655340
# FB_Igf1        0.9869201 1.0147943
# FB_Mfap5       1.1088832 0.8768457
# FB_Mmp2        0.9817855 1.0206018
# FB_Ogn         0.9908453 1.0103545
# FB_Postn       1.1213185 0.8627806
# Pericyte_Abcc9 1.1118635 0.8734748
# SMC_Acta2      1.1639877 0.8145188
# SMC_Cdh5       1.1569159 0.8225175
# SMC_Rgs5       1.1580446 0.8212408

# MiloR-------------------------------------------------------------------------
seurat.obj_sce <- as.SingleCellExperiment(seurat.obj)
seurat.obj_milo <- Milo(seurat.obj_sce)

# 1. Build a graph and neighbourhoods.------------------------------------------
# We need to add the KNN graph to the Milo object. This is stored in the graph slot, 
# in igraph format. The miloR package includes functionality to build and store
# the graph from the PCA dimensions stored in the reducedDim slot.

# The number of dimensions to use if the input is a matrix of cells X reduced dimensions. If this is provided, transposed should also be set=TRUE.
milo.obj <- buildGraph(seurat.obj_milo, k=20, d=20)
# 2. Defining representative neighbourhoods-------------------------------------
milo.obj <- makeNhoods(milo.obj, k=20, d=20, refined=TRUE, prop=0.2)
plotNhoodSizeHist(milo.obj)
# 3. Counting cells in neighbourhoods-------------------------------------------
require(SingleCellExperiment)
milo.obj <- countCells(milo.obj, 
                       meta.data = data.frame(colData(milo.obj)), 
                       samples="sample")
head(nhoodCounts(milo.obj))
# 4. Differential abundance testing---------------------------------------------
design <- data.frame(colData(milo.obj))[,c("sample", "condition")]
design <- distinct(design)
rownames(design) <- design$sample
## Reorder rownames to match columns of nhoodCounts(milo)
design <- design[colnames(nhoodCounts(milo.obj)), , drop=FALSE]

design

milo.obj <- calcNhoodDistance(milo.obj, d=20)
## 'as(<dgTMatrix>, "dgCMatrix")' is deprecated.
## Use 'as(., "CsparseMatrix")' instead.
## See help("Deprecated") and help("Matrix-deprecated").

# rownames(design) <- design$sample
da_results <- testNhoods(milo.obj, design = ~ condition, design.df = design)

# This calculates a Fold-change and corrected P-value for each neighbourhood, 
# which indicates whether there is significant differential abundance between 
# conditions.
da_results %>%
  arrange(- SpatialFDR) %>%
  head() 
# 5. Visualize neighbourhoods displaying DA-------------------------------------
milo.obj <- buildNhoodGraph(milo.obj)

scater::plotUMAP(milo.obj) + plotNhoodGraphDA(milo.obj, da_results, alpha=0.05) +
  patchwork::plot_layout(guides="collect")
# Score-------------------------------------------------------------------------
# The average expression of known proliferation-related genes was defined as the proliferation score
proliferation <- c('AURKA', 'BUB1', 'CCNB1', 'CCND1', 'CCNE1', 'DEK', 'E2F1', 'FEN1', 'FOXM1',
                   'H2AFZ', 'HMGB2', 'MCM2', 'MCM3', 'MCM4', 'MCM5', 'MCM6', 'MKI67', 'MYBL2',
                   'PCNA', 'PLK1', 'TOP2A', 'TYMS', 'ZWINT')
proliferation <- proliferation[proliferation %in% rownames(fig)]
proliferation_score <- fig@assays$RNA$data[proliferation,] %>% as.matrix %>% colMeans()
all(rownames(fig@meta.data) == names(proliferation_score))

fig$proliferation_score <- proliferation_score
VlnPlot(fig %>% subset(subset = cell_type != "Fibro_MKI67"), 
        features = "proliferation_score", group.by = "cell_type")
#===============================================================================
# SECTION V - Endothelial ======================================================
#===============================================================================
load("Result/scRNAseq/method4/scobj.harmony(singlet-withoutreg.cc-2000-25) anno.Rdata")

scobj.harmony@active.assay <- "RNA"

scobj.Endo <- scobj.harmony %>% subset(subset = cell_type %in% c("Endothelial"))

dim(scobj.Endo)
# [1] 21320 12520

scobj.Endo <- CreateSeuratObject(
  counts = scobj.Endo@assays$RNA$counts, 
  assay = "RNA", 
  meta.data = scobj.Endo@meta.data, 
  project = "Heart-Endothelial")
# QC----------------------------------------------------------------------------
table(scobj.Endo@meta.data$orig.ident)
# choline    sham 
# 4202    8318
# Nomalization & harmony & cluster----------------------------------------------

# 500 without regress cell cycle
if(TRUE){
  scobj.Endo.500 <- scobj.Endo %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 500) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony")
  
  scobj.Endo.500.10 <- scobj.Endo.500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:10) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Endo.500.15 <- scobj.Endo.500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:15) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Endo.500.20 <- scobj.Endo.500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Endo.500.25 <- scobj.Endo.500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Endo.500.30 <- scobj.Endo.500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Endo vst.500.h_clustree 10-30(singlet) without regression.pdf")
  ElbowPlot(scobj.Endo.500, reduction = "pca", ndims = 50)
  clustree(scobj.Endo.500.10, prefix = "RNA_snn_res.")
  clustree(scobj.Endo.500.15, prefix = "RNA_snn_res.")
  clustree(scobj.Endo.500.20, prefix = "RNA_snn_res.")
  clustree(scobj.Endo.500.25, prefix = "RNA_snn_res.")
  clustree(scobj.Endo.500.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.Endo.500, scobj.Endo.500.10, scobj.Endo.500.15,
       scobj.Endo.500.20, scobj.Endo.500.25, scobj.Endo.500.30,
       file = "scobj.Endo.500.Rdata")
}
# 1000 without regress cell cycle
if(TRUE){
  scobj.Endo.1000 <- scobj.Endo %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony")
  
  scobj.Endo.1000.10 <- scobj.Endo.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:10) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Endo.1000.15 <- scobj.Endo.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:15) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Endo.1000.20 <- scobj.Endo.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Endo.1000.25 <- scobj.Endo.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Endo.1000.30 <- scobj.Endo.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Endo vst.1000.h_clustree 10-30(singlet) without regression.pdf")
  ElbowPlot(scobj.Endo.1000, reduction = "pca", ndims = 50)
  clustree(scobj.Endo.1000.10, prefix = "RNA_snn_res.")
  clustree(scobj.Endo.1000.15, prefix = "RNA_snn_res.")
  clustree(scobj.Endo.1000.20, prefix = "RNA_snn_res.")
  clustree(scobj.Endo.1000.25, prefix = "RNA_snn_res.")
  clustree(scobj.Endo.1000.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.Endo.1000, scobj.Endo.1000.10, scobj.Endo.1000.15,
       scobj.Endo.1000.20, scobj.Endo.1000.25, scobj.Endo.1000.30,
       file = "scobj.Endo.1000.Rdata")
}
# 1500 without regress cell cycle
if(TRUE){
  scobj.Endo.1500 <- scobj.Endo %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony")
  
  scobj.Endo.1500.10 <- scobj.Endo.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:10) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Endo.1500.15 <- scobj.Endo.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:15) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Endo.1500.20 <- scobj.Endo.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Endo.1500.25 <- scobj.Endo.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Endo.1500.30 <- scobj.Endo.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Endo vst.1500.h_clustree 10-30(singlet) without regression.pdf")
  ElbowPlot(scobj.Endo.1500, reduction = "pca", ndims = 50)
  clustree(scobj.Endo.1500.10, prefix = "RNA_snn_res.")
  clustree(scobj.Endo.1500.15, prefix = "RNA_snn_res.")
  clustree(scobj.Endo.1500.20, prefix = "RNA_snn_res.")
  clustree(scobj.Endo.1500.25, prefix = "RNA_snn_res.")
  clustree(scobj.Endo.1500.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.Endo.1500, scobj.Endo.1500.10, scobj.Endo.1500.15,
       scobj.Endo.1500.20, scobj.Endo.1500.25, scobj.Endo.1500.30,
       file = "scobj.Endo.1500.Rdata")
}
# 2000 without regress cell cycle
if(TRUE){
  scobj.Endo.2000 <- scobj.Endo %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony")
  
  scobj.Endo.2000.10 <- scobj.Endo.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:10) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Endo.2000.15 <- scobj.Endo.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:15) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Endo.2000.20 <- scobj.Endo.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Endo.2000.25 <- scobj.Endo.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.Endo.2000.30 <- scobj.Endo.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Endo vst.2000.h_clustree 10-30(singlet) without regression.pdf")
  ElbowPlot(scobj.Endo.2000, reduction = "pca", ndims = 50)
  clustree(scobj.Endo.2000.10, prefix = "RNA_snn_res.")
  clustree(scobj.Endo.2000.15, prefix = "RNA_snn_res.")
  clustree(scobj.Endo.2000.20, prefix = "RNA_snn_res.")
  clustree(scobj.Endo.2000.25, prefix = "RNA_snn_res.")
  clustree(scobj.Endo.2000.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.Endo.2000, scobj.Endo.2000.10, scobj.Endo.2000.15,
       scobj.Endo.2000.20, scobj.Endo.2000.25, scobj.Endo.2000.30,
       file = "scobj.Endo.2000.Rdata")
}

scobj.Endo.h <- scobj.Endo.2000.10 %>% 
  RunUMAP(reduction = "harmony", dims = 1:10)
# Visualization-----------------------------------------------------------------
DimPlot(scobj.Endo.h, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.Endo.h, group.by = "RNA_snn_res.0.2", reduction = "umap", label = T)
# Annotation====================================================================
## Findmarkers------------------------------------------------------------------
Clist <- list()
n <- 1
for (i in 0:4) {
  Clist[[n]] <- FindMarkers(scobj.Endo.h, group.by = "RNA_snn_res.0.2", ident.1 = as.character(i))
  write.csv(Clist[[n]], paste0("D:/C",n-1,".csv"))
  
  n <- n + 1
}
## Annotation-------------------------------------------------------------------
scobj.Endo.h$cell_type <- "Unknown"
scobj.Endo.h$cell_type[scobj.Endo.h$RNA_snn_res.0.2 %in% c(0)] <- "EC_Cxcl9"
scobj.Endo.h$cell_type[scobj.Endo.h$RNA_snn_res.0.2 %in% c(1)] <- "EC_Peak1"
scobj.Endo.h$cell_type[scobj.Endo.h$RNA_snn_res.0.2 %in% c(2)] <- "EC_Hey1"
scobj.Endo.h$cell_type[scobj.Endo.h$RNA_snn_res.0.2 %in% c(3)] <- "EC_Dcn"
scobj.Endo.h$cell_type[scobj.Endo.h$RNA_snn_res.0.2 %in% c(4)] <- "EC_Myl9"

save(scobj.Endo.h, file = "scobj.Endo anno.Rdata")

table(scobj.Endo.h@meta.data$cell_type)
# EC_Cxcl9   EC_Dcn  EC_Hey1  EC_Myl9 EC_Peak1 
# 5850       94     2983       50     3543 
# Figure2=======================================================================
fig <- subset(scobj.Endo.h, subset = cell_type != "doublets")
fig$sample <- fig$orig.ident
fig$RNA_snn_res.0.2 <- factor(fig$RNA_snn_res.0.2, levels = c(0:4))
fig$cell_type <- factor(fig$cell_type,
                        levels = c("EC_Hey1","EC_Cxcl9","EC_Peak1",
                                   "EC_Dcn", "EC_Myl9"))
fig$orig.ident <- factor(fig$orig.ident, levels = c("sham", "choline"))

patient_col <- c("#c79494","#3d7d35")

p1 <- DimPlot(fig, reduction = "umap", group.by = "cell_type", cols = colors_list,
              label = TRUE, pt.size = 0.1, raster = FALSE)
p2 <- DimPlot(fig, reduction = "umap",
              group.by = "orig.ident", cols = patient_col, 
              label = TRUE, pt.size = 0.1, raster = FALSE)
p <- p1 + p2 + plot_layout(ncol = 2, nrow = 1)
ggsave(p, filename = "Endothlial UMAP.pdf", path = "Figure/", width=12, height=5, units="in")

p <- DotPlot(fig, features = c(
  # 'Thbd', 'Tek', 'Npr3',
  'Pecam1',"Eng",'Cdh5','Flt1','Rgcc',
  'Hey1','Fbln5','Fn1','Vegfc',
  'Cxcl9','Hdac9',  'Fabp5',
  'Peak1',
  'Dcn','Igfbp5','Vim', 'Vwf',
  'Myl6','Myl9',"Clu"
) %>% unique,
group.by = "cell_type") +
  scale_color_viridis() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave(p, filename = "Endothlial MARKER.pdf",
       path = "Figure/", width=6, height=6, units="in")

dat <- fig@meta.data %>% 
  dplyr::select(orig.ident, cell_type)
p <- dat %>% 
  mutate(value = 1) %>% 
  group_by(orig.ident, cell_type) %>%
  summarise(n = n()) %>% 
  group_by(orig.ident) %>%
  mutate(sum = sum(n)) %>% 
  mutate(n = n / sum) %>% 
  ggplot(aes(x = orig.ident, y = n, fill = cell_type, stratum = cell_type, alluvium = cell_type))+
  geom_stratum(width = 0.5, color='white') +
  geom_alluvium(alpha = 0.5,
                width = 0.5,
                curve_type = "linear") +
  scale_color_manual(values = colors_list) +
  scale_fill_manual(values = colors_list) +
  labs(x = "", y = "Percent", fill = "Cell type") +
  theme_bw()+
  theme(text = element_text(size = 20))
ggsave(p, filename = "Endothelial percent.pdf",
       path = "Figure/",
       width=5, height=6, units="in")
# Roe---------------------------------------------------------------------------
data <- fig@meta.data
data$majorCluster = data$cell_type
data$patient = data$sample
data$loc = data$orig.ident
Roe <- calTissueDist(data,
                     byPatient = F,
                     colname.cluster = "majorCluster", # 不同细胞亚群
                     colname.patient = "patient", # 不同样本
                     colname.tissue = "loc", # 不同组织
                     method = "chisq", # "chisq", "fisher", and "freq" 
                     min.rowSum = 0) 
# sham   choline
# EC_Hey1  1.0015962 0.9968403
# EC_Cxcl9 1.0325206 0.9356245
# EC_Peak1 0.9435454 1.1117537
# EC_Dcn   0.9447341 1.1094008
# EC_Myl9  1.2041356 0.5959067
#===============================================================================
# SECTION VI - Myeloid =========================================================
#===============================================================================
load("Result/scRNAseq/method4/scobj.harmony(singlet-withoutreg.cc-2000-25) anno.Rdata")

scobj.harmony@active.assay <- "RNA"

scobj.M <- scobj.harmony %>% subset(subset = cell_type %in% c("Myeloid"))

dim(scobj.M)
# [1] 21320   908

scobj.M <- CreateSeuratObject(
  counts = scobj.M@assays$RNA$counts, 
  assay = "RNA", 
  meta.data = scobj.M@meta.data, 
  project = "Heart-Myeloid")
# QC----------------------------------------------------------------------------
table(scobj.M@meta.data$orig.ident)
# choline    sham 
# 404     504 
# Nomalization & harmony & cluster----------------------------------------------

# 500 without regress cell cycle
if(TRUE){
  scobj.M.500 <- scobj.M %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 500) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony")
  
  scobj.M.500.10 <- scobj.M.500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:10) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.M.500.15 <- scobj.M.500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:15) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.M.500.20 <- scobj.M.500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.M.500.25 <- scobj.M.500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.M.500.30 <- scobj.M.500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Myeloid vst.500.h_clustree 10-30(singlet) without regression.pdf")
  ElbowPlot(scobj.M.500, reduction = "pca", ndims = 50)
  clustree(scobj.M.500.10, prefix = "RNA_snn_res.")
  clustree(scobj.M.500.15, prefix = "RNA_snn_res.")
  clustree(scobj.M.500.20, prefix = "RNA_snn_res.")
  clustree(scobj.M.500.25, prefix = "RNA_snn_res.")
  clustree(scobj.M.500.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.M.500, scobj.M.500.10, scobj.M.500.15,
       scobj.M.500.20, scobj.M.500.25, scobj.M.500.30,
       file = "scobj.M.500.Rdata")
}
# 1000 without regress cell cycle
if(TRUE){
  scobj.M.1000 <- scobj.M %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony")
  
  scobj.M.1000.10 <- scobj.M.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:10) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.M.1000.15 <- scobj.M.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:15) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.M.1000.20 <- scobj.M.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.M.1000.25 <- scobj.M.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.M.1000.30 <- scobj.M.1000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Myeloid vst.1000.h_clustree 10-30(singlet) without regression.pdf")
  ElbowPlot(scobj.M.1000, reduction = "pca", ndims = 50)
  clustree(scobj.M.1000.10, prefix = "RNA_snn_res.")
  clustree(scobj.M.1000.15, prefix = "RNA_snn_res.")
  clustree(scobj.M.1000.20, prefix = "RNA_snn_res.")
  clustree(scobj.M.1000.25, prefix = "RNA_snn_res.")
  clustree(scobj.M.1000.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.M.1000, scobj.M.1000.10, scobj.M.1000.15,
       scobj.M.1000.20, scobj.M.1000.25, scobj.M.1000.30,
       file = "scobj.M.1000.Rdata")
}
# 1500 without regress cell cycle
if(TRUE){
  scobj.M.1500 <- scobj.M %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony")
  
  scobj.M.1500.10 <- scobj.M.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:10) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.M.1500.15 <- scobj.M.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:15) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.M.1500.20 <- scobj.M.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.M.1500.25 <- scobj.M.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.M.1500.30 <- scobj.M.1500 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Myeloid vst.1500.h_clustree 10-30(singlet) without regression.pdf")
  ElbowPlot(scobj.M.1500, reduction = "pca", ndims = 50)
  clustree(scobj.M.1500.10, prefix = "RNA_snn_res.")
  clustree(scobj.M.1500.15, prefix = "RNA_snn_res.")
  clustree(scobj.M.1500.20, prefix = "RNA_snn_res.")
  clustree(scobj.M.1500.25, prefix = "RNA_snn_res.")
  clustree(scobj.M.1500.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.M.1500, scobj.M.1500.10, scobj.M.1500.15,
       scobj.M.1500.20, scobj.M.1500.25, scobj.M.1500.30,
       file = "scobj.M.1500.Rdata")
}
# 2000 without regress cell cycle
if(TRUE){
  scobj.M.2000 <- scobj.M %>% 
    subset(subset = dblfinder == "singlet") %>%
    NormalizeData(normalization.method = "LogNormalize") %>% # vst.flavor = 'v2', verbose = FALSE
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData %>% 
    RunPCA(npcs = 50) %>% 
    RunHarmony(
      group.by.vars = "orig.ident",
      reduction.use = "pca",
      reduction.save = "harmony")
  
  scobj.M.2000.10 <- scobj.M.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:10) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.M.2000.15 <- scobj.M.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:15) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.M.2000.20 <- scobj.M.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.M.2000.25 <- scobj.M.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:25) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  scobj.M.2000.30 <- scobj.M.2000 %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.1, 1, 0.1))
  
  pdf("Myeloid vst.2000.h_clustree 10-30(singlet) without regression.pdf")
  ElbowPlot(scobj.M.2000, reduction = "pca", ndims = 50)
  clustree(scobj.M.2000.10, prefix = "RNA_snn_res.")
  clustree(scobj.M.2000.15, prefix = "RNA_snn_res.")
  clustree(scobj.M.2000.20, prefix = "RNA_snn_res.")
  clustree(scobj.M.2000.25, prefix = "RNA_snn_res.")
  clustree(scobj.M.2000.30, prefix = "RNA_snn_res.")
  dev.off()
  
  save(scobj.M.2000, scobj.M.2000.10, scobj.M.2000.15,
       scobj.M.2000.20, scobj.M.2000.25, scobj.M.2000.30,
       file = "scobj.M.2000.Rdata")
}

scobj.M.h <- scobj.M.500.30 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)
# Visualization-----------------------------------------------------------------
DimPlot(scobj.M.h, group.by = "orig.ident", reduction = "umap", label = T)
DimPlot(scobj.M.h, group.by = "RNA_snn_res.1", reduction = "umap", label = T)
# Annotation====================================================================
## Findmarkers------------------------------------------------------------------
Clist <- list()
n <- 1
for (i in 0:5) {
  Clist[[n]] <- FindMarkers(scobj.M.h, group.by = "RNA_snn_res.1", ident.1 = as.character(i))
  write.csv(Clist[[n]], paste0("D:/C",n-1,".csv"))
  
  n <- n + 1
}
## Annotation-------------------------------------------------------------------
scobj.M.h$cell_type <- "Unknown"
scobj.M.h$cell_type[scobj.M.h$RNA_snn_res.1 %in% c(0)] <- "Macro_Cd86"
scobj.M.h$cell_type[scobj.M.h$RNA_snn_res.1 %in% c(1)] <- "doublet"
scobj.M.h$cell_type[scobj.M.h$RNA_snn_res.1 %in% c(2)] <- "Neu_S100a8"
scobj.M.h$cell_type[scobj.M.h$RNA_snn_res.1 %in% c(3)] <- "Mono_Ly6c2"
scobj.M.h$cell_type[scobj.M.h$RNA_snn_res.1 %in% c(4)] <- "Macro_Apoe"
scobj.M.h$cell_type[scobj.M.h$RNA_snn_res.1 %in% c(5)] <- "Mono_Cd300e"

save(scobj.M.h, file = "scobj.M anno.Rdata")

table(scobj.M.h@meta.data$cell_type)
# doublet  Macro_Apoe  Macro_Cd86 Mono_Cd300e  Mono_Ly6c2  Neu_S100a8 
# 186         119         218          76         154         155 
# Figure2=======================================================================
fig <- subset(scobj.M.h, subset = cell_type != "doublet")
fig$sample <- fig$orig.ident
fig$RNA_snn_res.1 <- factor(fig$RNA_snn_res.1, levels = c(0:5))
# fig$cell_type <- factor(fig$cell_type,
#                         levels = c("EC_Hey1","EC_Cxcl9","EC_Peak1",
#                                    "EC_Dcn", "EC_Myl9"))
fig$orig.ident <- factor(fig$orig.ident, levels = c("sham", "choline"))

patient_col <- c("#c79494","#3d7d35")

p1 <- DimPlot(fig, reduction = "umap", group.by = "cell_type", cols = colors_list,
              label = TRUE, pt.size = 0.1, raster = FALSE)
p2 <- DimPlot(fig, reduction = "umap",
              group.by = "orig.ident", cols = patient_col, 
              label = TRUE, pt.size = 0.1, raster = FALSE)
p <- p1 + p2 + plot_layout(ncol = 2, nrow = 1)
ggsave(p, filename = "Myeloid UMAP.pdf", path = "Figure/", width=12, height=5, units="in")

p <- DotPlot(fig, features = c(
  'Mrc1','C1qa','C1qb','C1qc','Mertk',
  'Apoe','Cd86', # 'Cd163', 
  'Cxcr4', 'Cd300e', 
  'Ly6c2', 'Ccr2', 
  'S100a8', 'S100a9', 'S100a11', 'Csf3r', 'G0s2', 'Ifitm1'
) %>% unique,
group.by = "cell_type") +
  scale_color_viridis() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave(p, filename = "Myeloid MARKER.pdf",
       path = "Figure/", width=6, height=6, units="in")

dat <- fig@meta.data %>% 
  dplyr::select(orig.ident, cell_type)
p <- dat %>% 
  mutate(value = 1) %>% 
  group_by(orig.ident, cell_type) %>%
  summarise(n = n()) %>% 
  group_by(orig.ident) %>%
  mutate(sum = sum(n)) %>% 
  mutate(n = n / sum) %>% 
  ggplot(aes(x = orig.ident, y = n, fill = cell_type, stratum = cell_type, alluvium = cell_type))+
  geom_stratum(width = 0.5, color='white') +
  geom_alluvium(alpha = 0.5,
                width = 0.5,
                curve_type = "linear") +
  scale_color_manual(values = colors_list) +
  scale_fill_manual(values = colors_list) +
  labs(x = "", y = "Percent", fill = "Cell type") +
  theme_bw()+
  theme(text = element_text(size = 20))
ggsave(p, filename = "Myeloid percent.pdf",
       path = "Figure/",
       width=5, height=6, units="in")


