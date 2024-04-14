#SCENIC.R 
#Jack Gordon 2023

#https://rawcdn.githack.com/aertslab/SCENIC/0a4c96ed8d930edd8868f07428090f9dae264705/inst/doc/SCENIC_Running.html#exploringinterpreting_the_results
#https://r.igraph.org/
#https://www.sc-best-practices.org/mechanisms/gene_regulatory_networks.html
#https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/scenic-differential-regulons.html#load-auc-socre-and-binary-mat

library(Seurat)
library(SeuratDisk)
library(SCENIC)
library(AUCell)
library(SCopeLoomR)
library(loomR)
library(ggplot2)
library(viridis)
library(dplyr)
library(ggpubr)
library(pheatmap)

#Read information from RNA loom file:
loom <- open_loom('pySCENIC/RNA/pyscenic_output.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name = 'Regulons')
RNA_regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
AUCmat <- AUCell::getAUC(regulonsAUC)
rownames(AUCmat) <- gsub("[(+)]", "", rownames(AUCmat))

#Load seurat object, add AUCmat as a new assay and save updated object
seurat <- readRDS("data/integrated.seurat.rds")
seurat[["RNA_AUC"]] <- CreateAssayObject(data = AUCmat)
saveRDS(seurat, "data/integrated.seurat.rds")

#Repeat for SCT loom file
loom <- open_loom('pySCENIC/SCT/pyscenic_output.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name = 'Regulons')
SCT_regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
AUCmat <- AUCell::getAUC(regulonsAUC)
rownames(AUCmat) <- gsub("[(+)]", "", rownames(AUCmat))

#Load seurat object, add AUCmat as a new assay and save updated object
seurat <- readRDS("data/integrated.seurat.rds")
seurat[["SCT_AUC"]] <- CreateAssayObject(data = AUCmat)
saveRDS(seurat, "data/integrated.seurat.rds")

#Read in seurat object with AUC matrix
seurat <- readRDS("data/integrated.seurat.rds")
seurat@active.ident <- factor(seurat$progenitor.origin)

#USE SCENIC TO PLOT HEATMAP
cellClusters <- seurat@meta.data
# Split the cells by integrated cluster to plot heatmap
cellsPerCluster <- split(rownames(cellClusters), cellClusters$integrated_snn_res.0.8)
regulonAUC <- regulonsAUC[onlyNonDuplicatedExtended(rownames(regulonsAUC)),]

# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

# Clean row names by removing '(+)'
cleaned_row_names <- gsub("\\(\\+\\)", "", rownames(regulonActivity_byCellType_Scaled))
rownames(regulonActivity_byCellType_Scaled) <- cleaned_row_names

# Plot heatmap
heatmap <- pheatmap(regulonActivity_byCellType_Scaled)

# Get the order of rows after clustering
row_order <- heatmap$tree_row$order

# Select labels for certain rows
selected_labels <- c("Dlx5", "Atf4", "Lhx8", "Lhx9", "Foxf2", "Neurod2", "Foxp2", 
                     "Otx2", "Pax6", "Setdb1", "E2f1", "Stat6", "Zfhx3", "Hdac2")

# Create a modified selected_row_names vector with empty space for non-displayed row names
selected_row_names <- ifelse(cleaned_row_names %in% selected_labels, cleaned_row_names, "")

# Draw the updated heatmap with modified row labels
heatmap <- pheatmap::pheatmap(
  regulonActivity_byCellType_Scaled,
  labels_row = selected_row_names,
  fontsize_row = 9  # Add legend title
)


########## Plot barplot of percent.change for different GRNs
#Calculate average AUC for each grn by progenitor origin
cellClusters <- seurat@meta.data
# Split the cells by integrated cluster to plot heatmap
cellsPerCluster <- split(rownames(cellClusters), cellClusters$progenitor.origin)
regulonAUC <- regulonsAUC[onlyNonDuplicatedExtended(rownames(regulonsAUC)),]

# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType <- data.frame(regulonActivity_byCellType)
rownames(regulonActivity_byCellType) <- gsub("[(+)]", "", rownames(regulonActivity_byCellType))
regulonActivity_byCellType$Difference <- regulonActivity_byCellType$aIP.derived - 
                                          regulonActivity_byCellType$OP.derived
regulonActivity_byCellType$percent.change <- ((regulonActivity_byCellType$aIP.derived - 
                                                 regulonActivity_byCellType$OP.derived) / 
                                                regulonActivity_byCellType$OP.derived) * 100

regulonActivity_byCellType$TF <- rownames(regulonActivity_byCellType)

#Perform DE testing - to filter for signifcant GRNs based on adjusted p_values
seurat@active.ident <- as.factor(seurat$progenitor.origin)
grns <- FindMarkers(seurat, ident.1 = "aIP-derived", 
                    ident.2 = "OP-derived", assay = "RNA_AUC", 
                    logfc.threshold = 0)

grns <- grns %>% filter(p_val_adj < 0.05)
grns.to.plot <- regulonActivity_byCellType %>% filter(TF %in% selected_TFs) %>%
  filter(TF %in% rownames(grns)) %>%
  mutate(progenitor.bias = ifelse(percent.change > 0, 
                                  "aIP-derived", "OP-derived")) %>%
  arrange(desc(TF))

#Plot percent.change for the DE GRNs
ggbarplot(grns.to.plot, x = "TF", y = "percent.change", 
          fill = "progenitor.bias", color = "black", 
          palette = c("forestgreen", "firebrick3"),
          xlab = "Gene Regulatory Network (GRN)", ylab = "Difference in GRN Activity (%)") + 
  theme_classic(base_size = 15) + 
  theme(legend.position = "none") +
  coord_flip()

###### TO DO #########
#Find top 50 GRNs in aIP and OP-derived SPNs
aIP.grns <- regulonActivity_byCellType %>% arrange(desc(aIP.derived)) %>% slice_head(n = 100)
OP.grns <- regulonActivity_byCellType %>% arrange(desc(OP.derived)) %>% slice_head(n = 100)

shared_row_names <- intersect(rownames(aIP.grns), rownames(OP.grns))
unique.aIP <- setdiff(rownames(aIP.grns), rownames(OP.grns))
unique.OP <- setdiff(rownames(OP.grns), rownames(aIP.grns))

selected_TFs <- c(shared_row_names, unique.aIP, unique.OP)
plot_grn_names <- intersect(rownames(grns.to.plot), selected_TFs)



###### TO DO #########
#IGRAPH w/ regulons
#Regulons stores each TF as a column and regulated genes in the column
#Creates network graph with igraph, with TF as nodes and genes as edges
library(igraph)

regulons <- RNA_regulons
selected_TFs <- c(shared_row_names, unique.aIP, unique.OP)

# Create an empty edgelist
edgelist <- matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c("TF", "Target")))

# Populate the edgelist with TF-gene connections
for (tf in names(regulons)) {
  tf_name <- gsub("\\(\\+\\)", "", tf)  # Remove '(+)' from TF names
  if (tf_name %in% selected_TFs) {
    for (gene in regulons[[tf]]) {
      edgelist <- rbind(edgelist, c(tf_name, gene))
    }
  }
}

# Create a directed graph from the edgelist
g <- graph_from_edgelist(as.matrix(edgelist), directed = TRUE)

# Get TFs that are regulated by other TFs
regulated_TFs <- V(g)$name[degree(g, mode = "in") > 0 & degree(g, mode = "out") > 0]

# Filter the edgelist to include only regulated TFs and their connections
filtered_edgelist <- edgelist[edgelist[,1] %in% regulated_TFs & edgelist[,2] %in% regulated_TFs, ]

# Create a new graph from the filtered edgelist
subgraph <- graph_from_edgelist(as.matrix(filtered_edgelist), directed = TRUE, )

# Get the induced subgraph of the selected TFs and their connections
selected_TF_names <- names(V(subgraph))  # Get names of vertices in the subgraph
selected_TF_indices <- which(selected_TF_names %in% selected_TFs)  # Find indices of selected TFs
subgraph <- induced_subgraph(subgraph, vids = selected_TF_indices)

# Plot using tkplot
plot <- plot(subgraph, vertex.label = V(subgraph)$name, vertex.size = 20, vertex.color = "lightblue", edge.arrow.size = 0.5)

tkplot(subgraph, vertex.label.cex = 2, label.font = "Arial", label.color = "black", edge.width = 0.5, 
       arrow.size = 0.5, vertex.size = 20)
tk_off()


#### Regulon specificity scores
AUC <- getAUC(regulonsAUC)

#Clean AUC regulon_names
current_regulon_names <- dimnames(AUC)[[1]]
cleaned_regulon_names <- sub("\\(\\+\\)$", "", current_regulon_names)
dimnames(AUC)[[1]] <- cleaned_regulon_names

#Filter for regulons of interest
regulon_names <- rownames(AUC)

# Identify indices of regulons that belong to selected transcription factors
selected_indices <- which(grepl(paste(selected_TFs, collapse = "|"), regulon_names))

# Filter the AUC matrix to include only regulons from selected transcription factors
filtered_AUC <- AUC[selected_indices, , drop = FALSE]


rss <- calcRSS(filtered_AUC, cellAnnotation=cellClusters[colnames(regulonsAUC), 
                                                                    "progenitor.origin"])

plotRSS_oneSet(rss, setName = "aIP-derived", n = 0) +
  xlab("GRN Rank") + ylab("GRN Specificity (AU)") + 
  scale_y_continuous(breaks = c(0.08, 0.09, 0.1, 0.11, 0.12),
                     labels = c("0.3", "0.4", "0.5", "0.6", "0.7"))
  theme_classic(base_size = 20)
sorted_indices <- order(rss[, "OP-derived"], decreasing = TRUE)
sorted_rss <- rss[sorted_indices, ]
