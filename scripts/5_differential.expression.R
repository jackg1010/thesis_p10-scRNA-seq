#Differential expression
library(Seurat)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)
library(VennDiagram)

#Read in Seurat object
integrated.seurat <- readRDS("data/integrated.seurat.rds")

#Assign to pseudo-bulk - update multi_ID to reflect the combinations of LMOs
htos <- t(Matrix::as.matrix(integrated.seurat@assays$HTO$data))
htos <- as.data.frame(htos)
htos$sample1 <- htos$`LMO-B1` + htos$`LMO-B2`
htos$sample2 <- htos$`LMO-B3`
htos$sample3 <- htos$`LMO-B4` + htos$`LMO-B5`
htos$sample4 <- htos$`LMO-B6` + htos$`LMO-B7`
htos <- htos %>% mutate(MULTI_ID = case_when(
  sample1 > sample2 & sample1 > sample3 & sample1 > sample4 ~ "Sample1", 
  sample2 > sample1 & sample2 > sample3 & sample2 > sample4 ~ "Sample2",
  sample3 > sample1 & sample3 > sample2 & sample3 > sample4 ~ "Sample3", 
  sample4 > sample1 & sample4 > sample2 & sample4 > sample3 ~ "Sample4"))

integrated.seurat$MULTI_ID <- htos$MULTI_ID
saveRDS(integrated.seurat, "data/integrated.seurat.rds")

#Read in seurat object
integrated.seurat <- readRDS("data/integrated.seurat.rds")

#Set progenitor origin as active identity
integrated.seurat@active.ident <- as.factor(integrated.seurat$progenitor.origin)

#Set RNA as default assay for DESeq2
DefaultAssay(integrated.seurat) <- "RNA"

#Create pseudo-bulk based on MULTI_ID, class and progenitor origin
pseudo <- AggregateExpression(integrated.seurat, return.seurat = T, 
                              group.by = c("MULTI_ID", "class", "progenitor.origin"))
pseudo@active.ident <- factor(pseudo$progenitor.origin)

#Create variable combining celltype (class) and progenitor origin for DE
pseudo$celltype.origin <- paste(pseudo$class, 
                              pseudo$progenitor.origin, sep = "_")

#Convert RNA counts to integers for DESeq2
DefaultAssay(pseudo) <- "RNA"
pseudo[["RNA"]]$counts <- round(pseudo[["RNA"]]$counts)

#Set variable as active ident
pseudo@active.ident <- factor(pseudo$celltype.origin)


#Pseudo-bulk differential expression testing with deseq2
bulk.de <- FindMarkers(pseudo, ident.1 = "NEURON_aIP-derived", 
                       ident.2 = "NEURON_OP-derived",
                       test.use = "DESeq2")

#Volcano plot
plot <-  EnhancedVolcano(bulk.de, lab = rownames(bulk.de),
                         x = 'avg_log2FC', y = 'p_val_adj', max.overlaps = Inf,
                         FCcutoff = 0.5) +
  theme(legend.position = "none")

#Arange as descending log fold change and add gene name as variable
bulk.de <- bulk.de %>% arrange(desc(avg_log2FC))
bulk.de$gene <- rownames(bulk.de)
bulk.de <- bulk.de %>% filter(p_val_adj < 0.05)

#Pseudo-bulk differential expression testing with deseq2
bulk.de.neurogenesis <- FindMarkers(pseudo, ident.1 = "NEUROGENESIS_aIP-derived", 
                       ident.2 = "NEUROGENESIS_OP-derived",
                       test.use = "DESeq2")

#Volcano plot
to.lab <- c("Il1rapl2", "Bhlhe22", "Eif1b", "Cfap54", "Snca", "Syne1", 
            "Pgm2l1", "Jph4", "Tubb2a", "Eef1a2", "Spef2", "Rprml", 
            "Crmp1", "Caly", "Aplp1", "Epha5", "L1cam", "Map1b", "Gabra1", 
            "Nrg3", "Mapt", "Grin2b", "Dner", "Syt4", "Scn2a", "Tbr1")

plot <-  EnhancedVolcano(bulk.de.neurogenesis, lab = rownames(bulk.de.neurogenesis),
                         x = 'avg_log2FC', y = 'p_val_adj', FCcutoff = 1, pCutoff = 0.05, 
                         selectLab = to.lab) +
  theme(legend.position = "none")

#Arange as descending log fold change and add gene name as variable
bulk.de.neurogenesis <- bulk.de.neurogenesis %>% arrange(desc(avg_log2FC))
bulk.de.neurogenesis$gene <- rownames(bulk.de.neurogenesis)
bulk.de.neurogenesis <- bulk.de.neurogenesis %>% filter(p_val_adj < 0.05, 
                                                        avg_log2FC >  1 | avg_log2FC < 1)


#ClusterProfiler ------------------------------------------------------------------------------
#Combine bulk.des 
bulk.de.merged <- merge(bulk.de, bulk.de.neurogenesis, by = "gene", all = TRUE)
bulk.de.merged <- bulk.de.merged %>% arrange(desc(bulk.de.merged))

# Calculate the average of avg_log2FC columns
bulk.de.merged$avg_log2FC <- rowMeans(bulk.de.merged[, c("avg_log2FC.x", "avg_log2FC.y")], na.rm = TRUE)

# Remove the redundant columns
bulk.de.merged <- bulk.de.merged[, !(names(bulk.de.merged) %in% c("avg_log2FC.x", "avg_log2FC.y"))]
write.csv(bulk.de.merged, "results/differential.expression.csv")

#Arange bulk.de by descending log2FC
bulk.de.merged <- bulk.de.merged %>% arrange(desc(bulk.de.merged$avg_log2FC))

#Convert gene symbols to entrez IDs for use in clusterprofiler
gene.ids <- bitr(bulk.de.merged$gene, fromType="SYMBOL", 
          toType= c("ENTREZID", "ENSEMBL"), OrgDb="org.Mm.eg.db")

#Create geneList - logFC change, name and sorted
geneList <- bulk.de.merged$avg_log2FC #numeric vector of logFC
names(geneList) <- gene.ids$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

#Select genes
gene <- names(geneList)

#Run gene ontology representation analysis
representation <- enrichGO(gene = gene,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "bonferroni",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

#Select significant terms to plot
rep.results <- representation@result
signif.terms <- rep.results$Description[rep.results$p.adjust < 0.05]
write.csv(signif.terms, "results/signif.GO.terms.csv")

#Plot neurogenesis results
p <- dotplot(representation)

ggsave("plots/GO/neurogenesis.dotplot.svg", p, "svg", 
       width = 10, height = 10, units = "in", dpi = "retina")

#Run gene set enrichment representation analysis
enrichment <- gseGO(geneList = geneList,
                        OrgDb        = org.Mm.eg.db,
                        ont          = "BP",
                        minGSSize    = 10,
                        maxGSSize    = 500,
                        pvalueCutoff = 0.05,
                        verbose      = FALSE)

#Plot gene ontology enrichment results
#Select significant terms to plot
enrich.results <- enrichment@result
signif.sets <- enrich.results$Description[enrich.results$p.adjust < 0.05]
write.csv(signif.sets, "results/signif.GSEA.sets.csv")


p <- dotplot(enrichment, showCategory = 
               c("neuron differentiation", 
                 "neuron projection morphogenesis", 
                 "synaptic signaling", 
                 "chemical synaptic transmission", 
                 "axon guidance", 
                 "regulation of synaptic plasticity", 
                 "synapse organization", 
                 "regulation of membrane potential", 
                 "vesicle-mediated transport in synapse", 
                 "dendrite development"))

ggsave("plots/GO/dotplot.enrichment.svg", p, "svg", 
       width = 10, height = 10, units = "in", dpi = "retina")


#Differentially expressed transcription factors
TFs <- read.table("pySCENIC/SCT/1.1_inputTFs.txt")
TFs <- TFs$V1
bulk_genes <- bulk.de.merged$gene

# Identify common genes - DEGs + TFS
common_genes <- intersect(TFs, bulk_genes)

#Venn diagram
library(RColorBrewer)
myCol <- brewer.pal(1, "Pastel2")

venn.plot <- venn.diagram(
  x = list(bulk_genes, TFs),
  category.names = c("" , 
                     ""),
  filename = "plots/venn.png", 
  lwd = 4,
  lty = 2,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)), 
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff'))
