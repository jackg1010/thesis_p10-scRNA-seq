#SCENIC.R 
#Jack Gordon 2023

#Runs SCENIC (https://github.com/aertslab/SCENIC) on integrated seurat object
#Follows workflow: http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html

#Install packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::version()
# If your bioconductor version is previous to 4.0, see the section beloww

## Required
BiocManager::install(c("AUCell", "RcisTarget"))
BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"))
# To export/visualize in http://scope.aertslab.org
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
#Load libraries
library(Seurat)
library(SeuratDisk)
library(SeuratObject)

#Import integrated seurat object, set assay to "SCT" and remove other assays
seurat <- readRDS("data/integrated.seurat.rds")
DefaultAssay(seurat) <- "RNA"

#Filter out genes that are expressed in less than 10 cells
seurat[["RNA"]] <- as(object = seurat[["RNA"]], Class = "Assay") #Convert SCT assay to v3
counts <- (seurat[["RNA"]]$counts)

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 1

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 5

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
filtered_counts <- Matrix::as.matrix(log2(filtered_counts+1), sparse = TRUE)

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = seurat@meta.data)
filtered_seurat[["RNA"]] <- as(object = filtered_seurat[["RNA"]], Class = "Assay")

#Export seurat object to loom for use in pySCENIC
dir.create("pySCENIC/RNA/")
seurat.loom <- as.loom(x = filtered_seurat, filename = "pySCENIC/RNA/seurat.loom", verbose = TRUE)
