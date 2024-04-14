library(Seurat)
library(data.table)
library(dplyr)
library(Matrix)

anderson <- fread("/Users/jackgordon/Downloads/GSE125290_AAX_ALL_CLUST_RAW_COUNTS.txt.gz")
anderson <- data.frame(anderson)
anderson <- anderson %>% filter(Genotype == "CTL") 

#Create meta variable for seurat object
meta <- anderson[,1:9]
meta$V1 <- NULL

#Add variable to meta by mapping cluster to name
meta <- meta %>% mutate(Cluster_Name = 
                          case_when(Cluster %in% c(3, 10, 12) ~ "Astrocyte",
                                    Cluster %in% c(11, 21, 43, 30, 42) ~ "Endothelial", 
                                    Cluster %in% c(18, 36, 41, 1, 6, 
                                                14, 22, 28, 17, 38, 
                                                13, 29, 27, 20, 24, 25) ~ "Neuron", 
                                    Cluster %in% c(0, 5, 7, 8, 19, 
                                               32, 33, 2, 34, 40) ~ "Progenitor", 
                                    Cluster %in% c(9, 23, 39, 35) ~ "Microglia", 
                                    Cluster %in% c(4, 15, 16, 31) ~ "Oligodendrocyte"))

#Create counts matrix for seurat object
counts <- t(anderson[,10:ncol(anderson)])
(counts) <- meta$CellBarcode

#Create seurat
anderson.seurat <- CreateSeuratObject(counts = counts, meta = meta)
saveRDS(anderson.seurat, "data/anderson.seurat.rds")
anderson <- readRDS("data/anderson.seurat.rds")

#Process anderson dataset with SCT
process.reference <- function(seurat) {
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
  seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, pattern = "^Rp[ls]")
  seurat <- NormalizeData(seurat, verbose = F)
  seurat <- FindVariableFeatures(seurat)
  seurat <- ScaleData(seurat)
  seurat <- RunPCA(seurat, verbose = F)
  seurat <- RunUMAP(seurat, dims = 1:30, verbose = F)
  seurat <- FindNeighbors(seurat, dims = 1:30, verbose = F)
  seurat <- FindClusters(seurat, verbose = T)
  mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", 
                 target_organism = "mmusculus")$ortholog_name
  mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", 
                   target_organism = "mmusculus")$ortholog_name #Get mouse cell cyle genes
  seurat <- CellCycleScoring(seurat, s.features = mmus_s, 
                             g2m.features = mmus_g2m, set.ident = FALSE)
  seurat$CC.difference <- seurat$S.Score - seurat$G2M.Score
  seurat <- SCTransform(seurat, vars.to.regress = c("percent.mt", 
                                                    "percent.rb", 
                                                    "CC.difference"), 
                        ncells = ncol(seurat))
  seurat <- RunPCA(seurat, verbose = T)
  seurat <- RunUMAP(seurat, dims = 1:30, verbose = T)
  seurat <- FindNeighbors(seurat, dims = 1:30, verbose = T)
  seurat <- FindClusters(seurat, verbose = T)
  return(seurat)
}
anderson <- process.reference(anderson)
saveRDS(anderson, "data/anderson.seurat.rds")

#Load in processed anderson data
anderson <- readRDS("data/anderson.seurat.rds")

#Load in query datasets
gfp.seurat <- readRDS("data/filtered.gfp.seurat.rds")
rfp.seurat <- readRDS("data/filtered.rfp.seurat.rds")


#--------------------------------------- GFP ----------------------------------------------
#Compute anchors for transfer labels
striatum.anchors <- FindTransferAnchors(reference = anderson, features = rownames(anderson),
                                        reference.assay = "SCT", reference.reduction = "pca",
                                        query = gfp.seurat, query.assay = "SCT")

#Use anchors for prediction
predictions <- TransferData(anchorset = striatum.anchors, refdata = anderson$Cluster_Name)

anderson.preds <- data.frame(predictions$predicted.id)
rownames(anderson.preds) <- rownames(predictions)
colnames(anderson.preds) <- "anderson.preds"

#Add predictions to query object
gfp.seurat <- AddMetaData(gfp.seurat, metadata = anderson.preds)
saveRDS(gfp.seurat, "data/filtered.gfp.seurat.rds")


#--------------------------------------- RFP ----------------------------------------------
#Compute anchors for transfer labels
striatum.anchors <- FindTransferAnchors(reference = anderson, features = rownames(anderson),
                                        reference.assay = "SCT", reference.reduction = "pca",
                                        query = rfp.seurat, query.assay = "SCT")

#Use anchors for prediction
predictions <- TransferData(anchorset = striatum.anchors, refdata = anderson$Cluster_Name)

anderson.preds <- data.frame(predictions$predicted.id)
rownames(anderson.preds) <- rownames(predictions)
colnames(anderson.preds) <- "anderson.preds"

#Add predictions to query object
rfp.seurat <- AddMetaData(rfp.seurat, metadata = anderson.preds)
saveRDS(rfp.seurat, "data/filtered.rfp.seurat.rds")
