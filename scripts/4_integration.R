#Integration of GFP and RFP samples

#For use after annotation.R 

#Load libraries
library(Seurat)
library(ggplot2)

#Set seed
set.seed(13)

#Load seurat objects
gfp.seurat <- readRDS("data/filtered.gfp.seurat.rds")
rfp.seurat <- readRDS("data/filtered.rfp.seurat.rds")

#Set active idents to class
gfp.seurat@active.ident <- factor(gfp.seurat$class)
rfp.seurat@active.ident <- factor(rfp.seurat$class)

#Subset seurats for neurons, neurogenesis, astrocyte classes
gfp.seurat <- subset(gfp.seurat, idents = c("NEURON", "NEUROGENESIS"))
rfp.seurat <- subset(rfp.seurat, idents = c("NEURON", "NEUROGENESIS"))

#Function integrate.seurats - integrates two seurat objects
integrate.seurats <- function(seurat1, seurat2, seurat1.label, seurat2.label) {
  DefaultAssay(seurat1) <- "SCT"
  DefaultAssay(seurat2) <- "SCT"
  seurat1$progenitor.origin <- seurat1.label
  seurat2$progenitor.origin <- seurat2.label
  seurats <- list(seurat1, seurat2)
  features <- SelectIntegrationFeatures(object.list = seurats, nfeatures = 5000)
  seurats <- PrepSCTIntegration(seurats, anchor.features = features)
  anchors <- FindIntegrationAnchors(object.list = seurats, anchor.features = features, 
                                    normalization.method = "SCT")
  integrated.seurat <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  dir.create("integration")
  saveRDS(features, "integration/integration.features.rds")
  saveRDS(anchors, "integration/integration.anchors.rds")
  saveRDS(integrated.seurat, "data/filtered.integrated.seurat.rds")
  return(integrated.seurat)
}


#Create integrated seurat object
integrated.seurat <- integrate.seurats(rfp.seurat, gfp.seurat, 
                                       seurat1.label = "OP-derived", 
                                       seurat2.label = "aIP-derived")


#Function - seurat.dimred - performs scaling, PCA, tSNE and UMAP 
seurat.dimred <- function(seurat, group.by, cols, filename) {
  DefaultAssay(seurat) <- "integrated"
  seurat <- ScaleData(seurat, verbose = TRUE)
  seurat <- RunPCA(seurat, npcs = 30, verbose = TRUE)
  seurat <- RunTSNE(seurat, reduction = "pca", dims = 1:30)
  seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:30)
  seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:30)
  seurat <- FindClusters(seurat, group.singletons = TRUE)
  tsne <- DimPlot(seurat, reduction = "tsne", group.by = group.by, cols = cols, 
                  alpha = 0.5, pt.size = 2) + theme_classic(base_size = 40)
  umap <- DimPlot(seurat, reduction = "umap", group.by = group.by, cols = cols, 
                  alpha = 0.5, pt.size = 2) + theme_classic(base_size = 40)
  ggsave(filename = paste("plots/integration/tsne", filename, sep="."), tsne, 
         "svg", width = 10, height = 10, units = "in", dpi = "retina")
  ggsave(filename = paste("plots/integration/umap", filename, sep="."), umap, 
         "svg", width = 10, height = 10, units = "in", dpi = "retina")
  saveRDS(seurat, "data/filtered.integrated.seurat.rds")
  return(seurat)
}

#Run dimensionality reduction on integrated object and save plots
integrated.seurat <- seurat.dimred(integrated.seurat, 
                                   group.by = "progenitor.origin", 
                                   cols = c("forestgreen", "firebrick3"), 
                                   filename = "integrated.svg")

#Visualise integration results
DimPlot(integrated.seurat, group.by = "integrated_snn_res.0.8", label = TRUE)
DimPlot(integrated.seurat, group.by = "progenitor.origin", label = TRUE)
DimPlot(integrated.seurat, group.by = "class", label = TRUE)

#Remove hypoxic cluster from integrated object
integrated.seurat <- subset(integrated.seurat, 
                            subset = integrated_snn_res.0.8 != 15 &
                              integrated_snn_res.0.8 != 10)
#Set default assay
DefaultAssay(integrated.seurat) <- "SCT"

#Split object back into progenitor origin datasets
gfp.seurat <- subset(integrated.seurat, 
                     subset = progenitor.origin == "aIP-derived")
rfp.seurat <- subset(integrated.seurat, 
                     subset = progenitor.origin == "OP-derived")

#Re-run sct and integration
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

gfp.seurat <- sct.cluster(gfp.seurat)
rfp.seurat <- sct.cluster(rfp.seurat)

#Re-run integration
integrated.seurat <- integrate.seurats(rfp.seurat, gfp.seurat, 
                                       seurat1.label = "OP-derived", 
                                       seurat2.label = "aIP-derived")

#Re-run dimensionality reduction
integrated.seurat <- seurat.dimred(integrated.seurat, 
                                   group.by = "progenitor.origin", 
                                   cols = c("forestgreen", "firebrick3"), 
                                   filename = "integrated.svg")


#Save and load updated object
saveRDS(integrated.seurat, "data/integrated.seurat.rds")
integrated.seurat <- readRDS("data/integrated.seurat.rds")


#Plot join UMAP on integrated object
p <- DimPlot(integrated.seurat, reduction = "umap", group.by = "progenitor.origin", 
        cols = c("forestgreen", "firebrick3"), pt.size = 1) + ggtitle(NULL) +
  theme(legend.position = "none")
  ggsave(filename = "plots/integration/umap.integrated.svg", p, 
       "svg", width = 10, height = 10, units = "in", 
       dpi = "retina")


