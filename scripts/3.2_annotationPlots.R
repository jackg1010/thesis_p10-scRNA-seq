library(Seurat)
library(dplyr)
library(ggpubr)
library(paletteer)
library(ggplot2)

set.seed(13)

#Plot tsne and umaps for cell type annotation
#Plot cell type composition

#-------------------------------- PLOT CLASS tSNEs + MARKER GENES ------------------------
gfp.seurat <- readRDS("data/filtered.gfp.seurat.rds")
rfp.seurat <- readRDS("data/filtered.rfp.seurat.rds")

#Subset for shared cell types
rfp.seurat <- subset(rfp.seurat, subset = class %in% (c("ASTROCYTE", 
                                                        "MACROPHAGE", 
                                                        "MICROGLIA",
                                                        "NEUROGENESIS",
                                                        "NEURON",
                                                        "OLIGODENDROCYTE", 
                                                        "POLYDENDROCYTE")))

gfp.seurat <- subset(gfp.seurat, subset = class %in% (c("ASTROCYTE", 
                                                        "MACROPHAGE", 
                                                        "MICROGLIA",
                                                        "NEUROGENESIS",
                                                        "NEURON",
                                                        "OLIGODENDROCYTE", 
                                                        "POLYDENDROCYTE")))

#Plot by class
p <- DimPlot(gfp.seurat, group.by = "class", reduction = "tsne") + 
  theme(legend.position = "none")
ggsave("plots/annotation/gfp.class.svg", p, width = 10, height = 10, units = "in", dpi = "retina")
p1 <- DimPlot(rfp.seurat, group.by = "class", reduction = "tsne") + 
  theme(legend.position = "none")
ggsave("plots/annotation/rfp.class.svg", p1, width = 10, height = 10, units = "in", dpi = "retina")


#Plot marker genes from anderson 2020
p <- FeaturePlot(gfp.seurat, reduction = "tsne", features = c("Ppp1r1b", "Sox4", "Gfap", "Olig1"))
ggsave("plots/annotation/gfp.markers.png", p, width = 10, height = 10, units = "in", dpi = "retina")
p <- FeaturePlot(rfp.seurat, reduction = "tsne", features = c("Ppp1r1b", "Sox4", "Gfap", "Olig1"))
ggsave("plots/annotation/rfp.markers.png", p, width = 10, height = 10, units = "in", dpi = "retina")


#--------------------------------------- Plot Cell Type Composition ----------------------------------------------
#set colour palettes
gfp.cols <- paletteer::paletteer_dynamic("cartography::green.pal", 7)
rfp.cols <- paletteer::paletteer_dynamic("cartography::red.pal", 7)

#Access meta data and calculate % cell type composition
gfp.meta <- data.frame(gfp.seurat@meta.data)
gfp.meta <- gfp.meta %>% select(predicted.id, class, common_name) %>% group_by(class) %>%
  dplyr::summarise(n = n()) %>% mutate(percent = n/sum(n)*100)

#Plot piechart of cell type composition (%)
p <- ggpie(gfp.meta, x = "percent", fill  = "class", palette = gfp.cols, label = "class", 
           label.pos = "out", lab.font = c(0, "bold", "white")) +
  theme(legend.position = "none") + xlab("") + ylab("") 
ggsave("plots/annotation/gfp.pie.svg", p, width = 10, height = 10, units = "in", dpi = "retina")


#Plot piechart of cell type composition (%)
rfp.meta <- data.frame(rfp.seurat@meta.data)
rfp.meta <- rfp.meta %>% select(predicted.id, class, common_name) %>% group_by(class) %>%
  dplyr::summarise(n = n()) %>% mutate(percent = n/sum(n)*100)

p <- ggpie(rfp.meta, x = "percent", fill  = "class", palette = rfp.cols, label = "class", 
           label.pos = "out", lab.font = c(0, "bold", "white")) + 
  theme(legend.position = "none") + xlab("") + ylab("")
ggsave("plots/annotation/rfp.pie.svg", p, width = 10, height = 10, units = "in", dpi = "retina")


#Access anderson 2020 meta data and calculate % cell type composition
gfp.meta <- data.frame(gfp.seurat@meta.data)
gfp.meta <- gfp.meta %>% select(anderson.preds) %>% group_by(anderson.preds) %>%
  dplyr::summarise(n = n()) %>% mutate(percent = n/sum(n)*100)

#Access anderson2020 meta data and calculate % cell type composition
rfp.meta <- data.frame(rfp.seurat@meta.data)
rfp.meta <- rfp.meta %>% select(anderson.preds) %>% group_by(anderson.preds) %>%
  dplyr::summarise(n = n()) %>% mutate(percent = n/sum(n)*100)
