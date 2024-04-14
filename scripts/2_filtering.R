#QC Filtering

#For use after processing.R

#Removes low quality cells with manually set thresholds - qc.filter
#Scores cell cycle - score.cell.cycle
#Normalises data with SCTrasnform - sct.cluster
  #regresses percent.mt, percent.rb and cell cycle score differences
#Finds doublets - doublet.finder
  #Removes doublets and saves filtered objects

#Load libraries
library(Seurat)
library(sctransform)
library(dplyr)
library(gprofiler2)
library(DoubletFinder)

#Remove low quality cells
qc.filter <- function(seurat.object, feature.count, max.mt.percent) {
  seurat.object <- subset(seurat.object, 
                          subset = percent.mt <= max.mt.percent & nFeature_RNA > feature.count)
  return(seurat.object)
}

#Filter rfp object and save filtered object
unfiltered.rfp <- readRDS("data/unfiltered.rfp.seurat.rds")
filtered.rfp <- qc.filter(unfiltered.rfp, feature.count = 700, max.mt.percent = 20)
saveRDS(filtered.rfp, "data/filtered.rfp.seurat.rds")

#Filter gfp object, save filtered object and make plots
unfiltered.gfp <- readRDS("data/unfiltered.gfp.seurat.rds")
filtered.gfp <- qc.filter(unfiltered.gfp, feature.count = 700, max.mt.percent = 20)
saveRDS(filtered.gfp, "data/filtered.gfp.seurat.rds")

#Perform default clustering (needed to score cell cycle)
cluster <- function(seurat) {
  seurat <- NormalizeData(seurat, verbose = F)
  seurat <- FindVariableFeatures(seurat)
  seurat <- ScaleData(seurat)
  seurat <- RunPCA(seurat, verbose = F)
  seurat <- RunUMAP(seurat, dims = 1:30, verbose = F)
  seurat <- FindNeighbors(seurat, dims = 1:30, verbose = F)
  seurat <- FindClusters(seurat, verbose = T)
  return(seurat)
}

#Cluster and update objects
filtered.rfp <- cluster(filtered.rfp)
filtered.gfp <- cluster(filtered.gfp)

#Score cell cycle for cells and calculate difference in cell cycle scores
#Cell cycle score
score.cell.cycle <- function(seurat) {
  mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", 
                 target_organism = "mmusculus")$ortholog_name
  mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", 
                   target_organism = "mmusculus")$ortholog_name #Get mouse cell cyle genes
  seurat <- CellCycleScoring(seurat, s.features = mmus_s, 
                             g2m.features = mmus_g2m, set.ident = FALSE)
  seurat$CC.difference <- seurat$S.Score - seurat$G2M.Score
  return(seurat)
}

filtered.rfp <- score.cell.cycle(filtered.rfp)
filtered.gfp <- score.cell.cycle(filtered.gfp)

#Perform SCT normalisation and use for clustering 
sct.cluster <- function(seurat) {
  seurat <- SCTransform(seurat, vars.to.regress = c("percent.mt", 
                                                    "percent.rb", 
                                                    "CC.difference"), 
                        verbose = T, ncells = ncol(seurat))
  seurat <- RunPCA(seurat, verbose = T)
  seurat <- RunTSNE(seurat, dims = 1:30)
  seurat <- RunUMAP(seurat, dims = 1:30, verbose = T)
  seurat <- FindNeighbors(seurat, dims = 1:30, verbose = T)
  seurat <- FindClusters(seurat, verbose = T)
  return(seurat)
}

#Cluster and update objects
filtered.rfp <- sct.cluster(filtered.rfp)
filtered.gfp <- sct.cluster(filtered.gfp)

#Use DoubletFinder to find doublets with different stringencies
doublet.finder <- function(seurat) {
  sweep.res <- paramSweep(seurat, PCs = 1:30, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bmcvn <- find.pK(sweep.stats)
  homotypic.prop <- modelHomotypic(seurat@meta.data$seurat_clusters)
  nExp_poi <- round(0.025*nrow(seurat@meta.data))  ## Assuming 2.5% doublet formation rate 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seurat <- doubletFinder(seurat, PCs = 1:10, pN = 0.25, pK = 0.09, 
                          nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
}

#Find doublets
filtered.gfp <- doublet.finder(filtered.gfp)
filtered.rfp <- doublet.finder(filtered.rfp)

#Rename doublet classifications
filtered.gfp$doublets <- filtered.gfp$DF.classifications_0.25_0.09_26
filtered.gfp$DF.classifications_0.25_0.09_26 <- NULL

#Rename doublet classifications
filtered.rfp$doublets <- filtered.rfp$DF.classifications_0.25_0.09_193
filtered.rfp$DF.classifications_0.25_0.09_193 <- NULL

#Remove doublets and save updated objects
filtered.gfp <- subset(filtered.gfp, subset = doublets == "Singlet")
filtered.rfp <- subset(filtered.rfp, subset = doublets == "Singlet")

saveRDS(filtered.gfp, "data/filtered.gfp.seurat.rds")
saveRDS(filtered.rfp, "data/filtered.rfp.seurat.rds")



