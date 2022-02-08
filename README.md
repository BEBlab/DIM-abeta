# Deep Indel Mutagenesis of amyloid beta
Pipeline to replicate the analysis and figures in: 
<a href="https://www.biorxiv.org/content/10.1101/2022.01.18.476804v1">An atlas of amyloid aggregation: the impact of substitutions, insertions, deletions and truncations on amyloid beta fibril nucleation</a>
(Mireia Seuma, Ben Lehner, Benedetta Bolognesi, 2022)

# System requirements
Software and packages to run the DIM-abeta pipeline :

- R v3.6.3 (ggplot2, dplyr, reshape2, weigths, stringr, readxl, ggpubr, pROC, ggrepel, data.table, RColorBrewer, grid, ggdendro, ggrepel, DescTools)
- DiMSum (https://github.com/lehner-lab/DiMSum; Faure, A.J., Schmiedel, J.M., Baeza-Centurion, P., Lehner B. DiMSum: an error model and pipeline for analyzing deep mutational scanning data and diagnosing common experimental pathologies. Genome Biol 21, 207 (2020))


# Installation
No installation required

# Demo and Instructions for use
- Raw sequencing data and the processed data table required for running the pipeline are deposited in NCBI's Gene Expression Omnibus (GEO) as GSE193837
- Script 01_processed data file.R outputs the processed data table from DiMSum output
- Script 02 outputs the complete data table for each mutation type subset for downstream analysis. All other scripts (03-16) can be run independently 

The run-time is <1h in total on a normal desktop computer.
