#Annotation R 
#Jack Gordon 2023

#For use after 2_filtering.R 

#Perform cell type annotation using Saunders 2018 - http://dropviz.org/
#https://satijalab.org/seurat/articles/integration_mapping.html#cell-type-classification-using-an-integrated-reference

library(Seurat)
library(DropSeq.util)
library(dplyr)
library(tidyr)
library(stringr)
library(gprofiler2)

#Read in count matrix and meta data for  reference dataset - create seurat object and save
dge.path <- "raw.data/saunders/F_GRCm38.81.P60Striatum.raw.dge (3).txt.gz"
counts <- loadSparseDge(dge.path) 

#Add cluster names to outcomes from the annotation dataframe to create metadata for seurat object
outcomes <- readRDS("raw.data/saunders/F_GRCm38.81.P60Striatum.cell_cluster_outcomes (2).RDS")
annotations <- readRDS("raw.data/saunders/annotation.BrainCellAtlas_Saunders_version_2018.04.01 (1).RDS")

annotations <- subset(annotations, tissue == "STR")
outcomes$class <- annotations$class[match(outcomes$subcluster, annotations$subcluster)]
outcomes$type_marker <- annotations$type_marker[match(outcomes$subcluster, annotations$subcluster)]
outcomes$common_name <- annotations$common_name[match(outcomes$subcluster, annotations$subcluster)]
outcomes$reason <- ifelse(is.na(outcomes$reason), "none", outcomes$reason)


#Create seurat object with counts and outcomes as meta data
saunders <- CreateSeuratObject(counts, meta = outcomes)
saunders <- subset(saunders, subset = reason %in% c(8, "none")) #Keep small cells and good quality cells
saveRDS(saunders, "data/saunders.seurat.rds")
saunders <- readRDS("data/saunders.seurat.rds")

#Run SCTransform on Saunders
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
saunders <- process.reference(saunders)
saveRDS(saunders, "data/saunders.rds")

#Read in processed reference dataset
saunders <- readRDS("data/saunders.rds")

#Load in query datasets
gfp.seurat <- readRDS("data/filtered.gfp.seurat.rds")
rfp.seurat <- readRDS("data/filtered.rfp.seurat.rds")

#--------------------------------------- GFP ----------------------------------------------
#Compute anchors for transfer labels
striatum.anchors <- FindTransferAnchors(reference = saunders, features = rownames(saunders),
                                        reference.assay = "SCT", 
                                        reference.reduction = "pca",
                                        query = gfp.seurat, query.assay = "SCT")

#Use anchors for prediction
predictions <- TransferData(anchorset = striatum.anchors, refdata = saunders$type_marker)

predictions$class <- saunders$class[match(predictions$predicted.id, saunders$type_marker)]
predictions$common_name <- saunders$common_name[match(predictions$predicted.id, saunders$type_marker)]

metadata <- predictions %>% select(predicted.id, class, common_name)

#Add predictions to query object
gfp.seurat <- AddMetaData(gfp.seurat, metadata = metadata)
saveRDS(gfp.seurat, "data/filtered.gfp.seurat.rds")


#--------------------------------------- RFP ----------------------------------------------
#Compute anchors for transfer labels
striatum.anchors <- FindTransferAnchors(reference = saunders, features = rownames(saunders),
                                        reference.assay = "SCT", 
                                        reference.reduction = "pca",
                                        query = rfp.seurat, query.assay = "SCT")

#Use anchors for prediction
predictions <- TransferData(anchorset = striatum.anchors, refdata = saunders$type_marker)

predictions$class <- saunders$class[match(predictions$predicted.id, saunders$type_marker)]
predictions$common_name <- saunders$common_name[match(predictions$predicted.id, saunders$type_marker)]

metadata <- predictions %>% select(predicted.id, class, common_name)

#Add predictions to query object
rfp.seurat <- AddMetaData(rfp.seurat, metadata = metadata)
saveRDS(rfp.seurat, "data/filtered.rfp.seurat.rds")






