library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(ggplot2)
library(ggpubr)

set.seed(13)

#https://rpubs.com/mahima_bose/Seurat_and_Monocle3_p

#Load in integrated seurat
integrated.seurat <- readRDS("data/integrated.seurat.rds")

#Visualise clustering and expression of markers
#Visualise UMAP of clusters - find Mki67+ root clusters
#Cluster no4 is root in SCT_snn_res.0.8
DefaultAssay(integrated.seurat) <- "SCT"
FeaturePlot(integrated.seurat, features = c("Sox4", "Mki67", "Snap25", "Foxp1"))
p <- DimPlot(integrated.seurat, group.by = c("class")) + theme(legend.position = "none")
p1 <- DimPlot(integrated.seurat, group.by = c("Phase")) + theme(legend.position = "none")
both <- p + p1


#Set idents to SCT clusters 
Idents(integrated.seurat) <- "SCT_snn_res.0.8"

#Convert integrated seurat object into monocle object
cds <- SeuratWrappers::as.cell_data_set(integrated.seurat)

#Update rowdata in integrated object
fData(cds)$gene_short_name <- rownames(fData(cds))
head(fData(cds))

#Assign partitions to cds object set all cells as 1st partition
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

#Manually add clustering and UMAP info to monocle object - clusters, co-ordinates
list.cluster <- integrated.seurat@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- 
integrated.seurat@reductions$umap@cell.embeddings

#Plot pre.trajectory
cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", 
                                 label_groups_by_cluster = F, 
                                 group_label_size = 5) + 
  theme(legend.position = "right")

#Learn trajectory and plot
cds <- learn_graph(cds, use_partition = F)

plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)

#Order cells in pseudo-time
cds <- order_cells(cds, reduction_method = "UMAP", 
                  root_cells = colnames(cds[, clusters(cds) == 4]))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F)

#Plot pseudo-time for neurogenesis vs. neuronal cells
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

p <- ggbarplot(data = data.pseudo, 
                x = "class", y = "monocle3_pseudotime", 
                color = "black", fill = "class",
                order = c("NEUROGENESIS", "NEURON"),
                add = "mean_se",
                error.plot = "upper_errorbar", 
                add.params = list(color = "black", alpha = 0.5)) +
  theme_classic(base_size = 40) + 
  theme(legend.position = "none") + xlab(NULL) + ylab("pseudotime") +
  stat_compare_means(method = "wilcox", label = "p.signif")
ggsave("plots/pseudo-time/class.barplot.svg", p, "svg", 
       width = 10, height = 10, units = "in", dpi = "retina")


#Find genes that vary as a function of pseudotime
deg <- graph_test(cds, neighbor_graph = "principal_graph")
deg %>% arrange(q_value) %>% filter(status == "OK") %>% 
  filter(q_value < 0.05) %>% head(10)

#Subset for selected genes and plot change over pseudotime
cds <- estimate_size_factors(cds)
my_genes <- row.names(subset(fData(cds), 
                             gene_short_name %in% c("Foxp1")))
cds_subset <- cds[my_genes,]

p <- plot_genes_in_pseudotime(cds_subset, cell_size = 0.1,
                              vertical_jitter = 0.5) + 
  theme_classic(base_size = 40) + 
  theme(legend.position = "none") + ylab("Expression") 
ggsave("plots/pseudo-time/foxp1.new.svg", p, "svg", 
       width = 10, height = 10, units = "in", dpi = "retina")


#Add variable combining class and progenitor.origin
integrated.seurat$celltype.origin <- paste(integrated.seurat$class, 
                                           integrated.seurat$progenitor.origin, 
                                           sep = "_")

integrated.seurat@active.ident <- factor(integrated.seurat$celltype.origin)

#Add pseudotime as metadata column to seurat object
integrated.seurat$pseudotime <- pseudotime(cds)
saveRDS(integrated.seurat, "data/integrated.seurat")

#Comparisons of aIP and OP-derived pseudo-time
my.comparisons <- list(c("aIP-derived", "OP-derived"))
meta <- integrated.seurat@meta.data
meta$pseudotime <- as.numeric(meta$pseudotime)
meta$progenitor.origin <- as.factor(meta$progenitor.origin)

neurogenesis.meta <- meta %>% filter(celltype.origin %in% 
                                       c("NEUROGENESIS_OP-derived", 
                                         "NEUROGENESIS_aIP-derived"))

neuron.meta <- meta %>% filter(celltype.origin %in% 
                                 c("NEURON_OP-derived", 
                                   "NEURON_aIP-derived"))


#Calculate results
str(meta$pseudotime)
pseudo.res <- meta %>% select(progenitor.origin, pseudotime) %>%
group_by(progenitor.origin) %>% 
  summarise(mean = mean(pseudotime), 
  sem = plotrix::std.error(pseudotime), 
  median = median(pseudotime), mad = mad(pseudotime))


p <- ggdensity(meta, x = "pseudotime",
          color = "progenitor.origin", fill = "progenitor.origin", 
          palette = c("forestgreen", "firebrick3"), 
          add = "median", add.params = list(size = 1.5, linetype = "dashed"), rug = TRUE) + 
  theme_classic(base_size = 50) + theme(legend.position = "none") +
  ylab("Density") 
ggsave("plots/pseudo-time/density.plot.svg", p, "svg", 
                           width = 25, height = 15, units = "in", dpi = "retina")

result <- meta %>% dplyr::select(progenitor.origin, pseudotime) %>%
  summarise(across(!progenitor.origin, ~wilcox.test(.x ~ progenitor.origin)$p.value)) 

