# Spinal cord paper

This repository contains the code to recreate the results and figures of our current publication " " xxx.

We calculate the tSNE embeddings with the FFT (Fast Fourier Transform) accerelated Interpolation-based t-SNE [FIt-SNE](https://github.com/KlugerLab/FIt-SNE) developed by the Kluger lab. Install it to a directory to your machine following their instruction.

The Initial QC filtering and seurat object creation markdown files are rendered using an array job in the scripts/ directory. Use the following order:

-   ***QCfilter_array.sh***: Initial removal of ambient RNA, quality control and filtering. Follows the methods described in [Feregrino et al. 2019](https://doi.org/10.1186/s12864-019-5802-2).
-   ***Seurat_Gg_NT_array.sh***: create and analyse [Seurat](https://github.com/satijalab/seurat) objects for all the individual samples, run dimensionality reduction, marker detection and plot markers for cluster annotation.
The chickent and mouse ortholog table (~/spinal_cord_paper/data/ortho_gg_mm_v102.rds) is retrieved from [BioMart](http://nov2020.archive.ensembl.org/biomart/martview/) (release 102).
-   ***Seurat_integration_array.sh***: Integration step of the 4 integrated data sets
    -   Devel: D05_ctrl, D07_ctrl, and ctrl_1 (day 10) brachial samples
    -   brachial: ctrl_1 and ctrl_2
    -   lumbar: lumb_1 and lumb_2
    -   polydactyl: poly_1 and poly_2
-   ***Seurat_Gg_NT_int_array.sh***: Create and analyse the integrated seurat objects, run dimensionality reduction, marker detection and plot markers for cluster annotation.
-   ***Seurat_Gg_ctrl_lumb_ctrl_poly_integration_array.sh***: Integration and analysis of the brachial & lumbar (ctrl_1, ctrl_2, lumb_1, & lumb_2), and brachial & poly (ctrl_1, ctrl_2, poly_1, & poly_2). They get an individual, combined pipeline since the time and memory allocation are higher.

Cluster annotation is guided by marker gene expression:

-   ***Broad_clusters_and_DV_domain_plots_render.sh***: Plots marker gene expression of all data sets and exports them to annotations/figures. Slurmed by *annotations/Broad_clusters_and_DV_domain_plots_render.Rmd*.

## scWGNCA

Module construction and co-expression analysis is done with [scWGCNA](https://github.com/CFeregrino/scWGCNA). The markdown scripts are again called with and array job:

-   ***scWGCNA_array.sh***: Run scWGCNA on the integrated ctrl, lumb, poly, and devel data sets. Produces the Gg_*_int_scWGCNA_modules and _modules_expression.pdf plots.

Next we run the comparative scWGCNA pipeline to calculate module conservation between the data sets:

-   ***comparative_scWGCNA_array.sh***: compares each data set's modules with the other two data sets. Produces the "Gg_*_int_comp_scWGCNA.pdf" plots.

## DA (differential abundance)

Differential abundance is calculated with [DAseq](https://github.com/KlugerLab/DAseq) and [miloR](https://github.com/MarioniLab/miloR). Again the scripts are called with separate array jobs:

-   ***DAseq_Gg_ctrl_lumb_cltr_poly_int_array.sh***: Runs DAseq to calculate differential abundance. Requires a python3 verison available (install with anaconda).
-   ***miloR_Gg_ctrl_lumb_cltr_poly_int_array.sh***: Runs miloR to calculate differential abundance.

## augur Cell Type Prioritisation

[augur](https://github.com/neurorestore/Augur) is used to prioritise cell types based on ther response to an experimental condition:

-   ***augur_Gg_ctrl_lumb_ctrl_poly_int_array.sh***:  Runs augur between ctrl and lumb, as well as ctrl and poly. 

## Sister Pair DE analysis

This script runs DE analysis of closest neuronal clusters between ctrl & lumb, as well as ctrl & poly. Not only strict sister pairs, but also the next closest cluster from the other sample was considered.

-   ***Sister_pair_DE_analysis_render.sh***: Loads the ctrl, lumbar and poly int data sets, transfers the cluster labels to the ctrl_lumb_int and ctrl_poly_int respectively. Then on the integrated data, DE between the sister pairs (or adjacent clusters) as by the spearman cor heatmaps ordering.

## Plotting

The following scripts are used to create the plots for the figures, if they are not already produced by the pipelines above:

-   ***dotplots_broad_render.sh***: Plots the dotplots that show average marker gene expression by cluster.

-   ***tsne_and_bar_plots_render.sh***: Plots the tsne embeddings colored by broad clustering as well as stacked bar plot by data set to show the relative cell type contributions.

-   ***heatmap_spearman_ctrl_lumb_poly_int.sh*** and ***heatmap_spearman_devel.sh***: Produces the heatmaps of cluster and module correlation.

-   ***Vln_plots_render.sh***: Plots diagnostic violin/beeswarm plots to quickly reference mitochondrial and ribosomal gene expression by cluster.

-   ***Gg_devel_scWGCNA_module_analysis_render.sh***: All the analysis of the scWGCNA modules on the development (D5, D7, and D10/ctrl1) data set.
