#pySCENIC.py 

import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import seaborn as sns
import matplotlib.pyplot as plt
import pyscenic

# set a working directory and filepaths
wdir = "/Users/jackgordon/Documents/DPhil/RNA-Sequencing/P10-seq/pySCENIC"
os.chdir( wdir )

MOTIFS_MGI_FNAME = 'motifs/motifs-v9-nr.mgi-m0.001-o0.0.tbl'
OUT_TFS_MGI_FNAME = 'SCT/inputTFs.txt'

#GRN inference with GRNBoost2
!arboreto_with_multiprocessing.py  '/RNA/seurat.loom' 'SCT/inputTFs.txt' -o adj.csv --num_workers 8 -m grnboost2

#Regulon prediction with cisTarget
!pyscenic ctx 'adj.csv' 'mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather' 'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather' --annotations_fname "motifs-v9-nr.mgi-m0.001-o0.0.tbl" --expression_mtx_fname 'seurat.loom' --output 'reg.csv' --mask_dropouts --num_workers 8

#Cellular enrichment with AUCell
!pyscenic aucell 'RNA/seurat.loom' 'RNA/reg.csv' --output 'RNA/pyscenic_output.loom' --num_workers 8

