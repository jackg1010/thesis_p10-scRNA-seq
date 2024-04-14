#qc plots

library(Seurat)
library(ggplot2)
library(ggpubr)
library(plotrix)
library(dplyr)

#Load seurat object
integrated.seurat <- readRDS("data/integrated.seurat.rds")

#Plot and save QC metrics
p <- VlnPlot(integrated.seurat, group.by = "progenitor.origin", features = "nCount_RNA", 
             pt.size = 0.1, log = TRUE, cols = c("forestgreen", "firebrick3")) +
  ylab("Mapped Reads (log)") + ggtitle(NULL) + xlab("") +
  theme_classic(base_size = 40) + theme(legend.position = "none")
ggsave("plots/QC/mapped.reads.svg", p, "svg", width = 10, height = 10, units = "in", dpi = "retina")

p <- VlnPlot(integrated.seurat, group.by = "progenitor.origin", features = "nFeature_RNA", 
             pt.size = 0.1, log = TRUE, cols = c("forestgreen", "firebrick3")) +
  ylab("Detected Genes (log)") + ggtitle(NULL) + xlab("") +
  theme_classic(base_size = 40) + theme(legend.position = "none")
ggsave("plots/QC/detected.genes.svg", p, "svg", width = 10, height = 10, units = "in", dpi = "retina")

p <- VlnPlot(integrated.seurat, group.by = "progenitor.origin", features = "percent.mt", 
             pt.size = 0.1, log = FALSE, cols = c("forestgreen", "firebrick3")) +
  ylab("Mitochondrial Reads (%)") + ggtitle(NULL) + xlab("") +
  theme_classic(base_size = 40) + theme(legend.position = "none")
ggsave("plots/QC/percent.mt.svg", p, "svg", width = 10, height = 10, units = "in", dpi = "retina")


#Calculate QC metrics and save as .csv
metadata <- as.data.frame(integrated.seurat@meta.data)
write.csv(metadata, "results/qc/metadata.csv")
metadata <- as.data.frame(read.csv("results/qc/metadata.csv"))

qcmetrics <- metadata %>% dplyr::select(progenitor.origin, nCount_RNA, nFeature_RNA, percent.mt, percent.rb) %>% 
  group_by(progenitor.origin) %>%
  summarize(across(everything(), list(mean = mean, sem = std.error)))
write.csv(qcmetrics, "results/qc/qcmetrics.csv")

#Stats testing
result <- metadata %>% dplyr::select(progenitor.origin, nCount_RNA, nFeature_RNA, percent.mt, percent.rb) %>%
  summarise(across(!progenitor.origin, ~wilcox.test(.x ~ progenitor.origin)$p.value)) 


#Cell cycle
cell.cycle <- metadata %>% group_by(progenitor.origin, Phase) %>% 
  dplyr::summarise(n = n()) %>% dplyr::mutate(percent = n/sum(n)*100)
  
plot <- ggbarplot(cell.cycle, x = "progenitor.origin", 
                  order = c("aIP-derived", "OP-derived"), y = "percent", fill = "Phase", merge = TRUE) +
  ylab("% of Cells") + xlab("") + theme_classic(base_size = 40) + theme(legend.position = "none")
ggsave("plots/QC/cell.cycle.svg", plot, "svg", width = 10, height = 10, units = "in", dpi = "retina")

