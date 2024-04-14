# thesis_p10-scRNA-seq
Repository hosts the analysis code for P10 scRNA-seq data of aIP and OP-derived SPNs.

Pipeline:

  1. Processing
  - Creates seurat objects from Cell Ranger outputs
  - Uses SoupX to remove RNA contam. and updates count matrices

  2. Filtering
   - Removes low quality cells
   - Scores cell cycle
   - Runs SCTransform 
   - Removes doublets 
   - Plots qc metrics

  3. Cell Type Annotation
   - Annotates cell types w/ Saunders 2018 data
   - Adds cell type predictions to metadata
   - Saves updated seurat objects
   - Plots w/ class labels and cell type composition

  4. Dataset Integration
   - Removes non-neuronal cells from both objects
   - Integrates GFP and RFP seurat objects
   - Runs integrated dimensionality reduction and clustering

  5. Differential Expression
   - Pseudo-bulk DESEq2 for GFP and RFP neurons
   - Runs gene ontology and gene set enrichment

  6. SCENIC
   - Converts seurat to loom object for use with pyscenic
   - Runs pySCENIC with CLI
   - Analysis of AUC matrix from pySCENIC output
 
 7. Pseudo-time
   - Trajectory with Monocle3
   - Scoring

