#### Pipeline for Analysis of GFP+ and RFP+ P10 sequencing data ###

#1_processing.R
# - Creates seurat objects from Cell Ranger outputs
# - Uses SoupX to remove RNA contam. and updates count matrices
# - Saves seurat objects 

#2_filtering.R
# - Removes low quality cells
# - Scores cell cycle
# - Runs SCTransform 
# - Removes doublets 
# - Saves filtered seurat objects

#2.1_qc.plots.R
# - Plots qc metrics

#3_annotation.R
# - Annotates cell types w/ Saunders 2018 data
# - Adds cell type predictions to metadata
# - Saves updated seurat objects
# - Plots w/ class labels and cell type composition

#3.1_annotation2.R
# - Annotates cell types w/ Anderson 2020 data
# - Adds cell type predictions to metadata
# - Saves updated seurat objects
# - Plots w/ class labels and cell type composition

#4_integration.R
# - Removes non-neuronal cells from both objects
# - Integrates GFP and RFP seurat objects
# - Runs integrated dimensionality reduction and clustering

#5.1_DEandGO.R
# - Finds markers for GFP and RFP neurons
# - Runs gene ontology and gene set enrichment (cluster profiler)

#5.2_pseudobulk.R
# - Runs pseudobulk differential expression find markers
# - Runs gene ontology and gene set enrichment (cluster profiler)

#6.1_seurat to loom.R 
# - Converts seurat to loom object for use with pyscenic

#6.2_psyscenic.py
# - Commands for command line pySCENIC run (done on cluster)

#6.3_gene.reg.networks.R
# - Analysis of AUC matrix from pySCENIC output





