#Seurat QC Filtering R 
#Jack Gordon 2023

#Imports cellranger outputs from two FACS-10X Genomics runs
  #TEL001 - GFP FACS-sorted SPNs
  #TEL002 - tdTomato FACS-sorted SPNs
#Creates seurat objects and uses SoupX to remove contamination from ambient RNA (estimated with RFP sample)
#Adds cell cycle scores and difference between cell cycle scores for each cell in each seurat

library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)
library(DropletUtils)
library(SoupX)
library(scater)
library(sctransform)
library(ggplot2)
library(dplyr)
library(DoubletFinder)
library(gprofiler2)


#Function create.seurat - reads in the data and creates a seurat object demultiplexed with sample ID stored in MUTLI_ID column
create.seurat <- function(directory) {
  data <- Read10X(data.dir = directory) #Read in data
  seurat <- CreateSeuratObject(counts = data$`Gene Expression`) #Create seurat object with counts
  seurat[["HTO"]] <- CreateAssayObject(counts = data$`Antibody Capture`) #Add antibidoy tags as an assay
  seurat <- NormalizeData(seurat, assay = "HTO", normalization.method = "CLR") #Normalise antibody tags
  seurat <- MULTIseqDemux(object = seurat, assay = "HTO", 
                          autoThresh = TRUE, maxiter = 20) #Perform demulitplexing based on antibody tags - adds $MULTI_ID
  seurat <- SetIdent(seurat, value = "orig.ident")
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
  seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, pattern = "^Rp[ls]")
  return(seurat) #Return seurat object
}

#Create rfp seurat object
rfp.dir <- "raw.data/TEL002/deep.seq/outs/filtered_feature_bc_matrix"
rfp.seurat <- create.seurat(rfp.dir)

#Use SoupX to remove ambient RNA
rfp.sc = load10X("raw.data/TEL002/deep.seq/outs") #Load in 10X output for RFP
rfp.sc = autoEstCont(rfp.sc) #Load in 10X output for RFP
rfp.out = adjustCounts(rfp.sc) #Load in 10X output for RFP

#Update seurat object with SoupX cleaned count matrix and save
rfp.seurat[["RNA"]]$counts <- rfp.out
saveRDS(rfp.seurat, "data/unfiltered.rfp.seurat.rds")

#Create gfp seurat object
gfp.dir <- "raw.data/TEL001/deep.seq/outs/filtered_feature_bc_matrix"
gfp.seurat <- create.seurat(gfp.dir)

#Use SoupX to remove ambient RNA
gfp.sc = load10X("raw.data/TEL001/deep.seq/outs") #Load in 10X output for GFP
gfp.sc = setContaminationFraction(gfp.sc, 0.14) #Set from RFP soup channel
gfp.out = adjustCounts(gfp.sc) #Load in 10X output for RFP

#Update seurat object with SoupX cleaned count matrix and save
gfp.seurat[["RNA"]]$counts <- gfp.out
saveRDS(gfp.seurat, "data/unfiltered.gfp.seurat.rds")










