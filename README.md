# Spinal cord paper

This repository contains the code to recreate the results and figures of our current publication " " xxx.

We calculate the tSNE embeddings with the FFT (Fast Fourier Transform) accerelated Interpolation-based t-SNE [FIt-SNE](https://github.com/KlugerLab/FIt-SNE) developed by the Kluger lab. Install it to a directory to your machine following their instruction.

The Initial QC filtering and seurat object creation markdown files are rendered using an array job in the scripts/ directory. Use the following order:

-   ***QCfilter_array.sh***: Initial removal of ambient RNA, quality control and filtering. Follows the methods described in [Feregrino et al. 2019](https://doi.org/10.1186/s12864-019-5802-2).
-   ***Seurat_Gg_NT_array.sh***: create and analyse [Seurat](https://github.com/satijalab/seurat) objects for all the individual samples, run dimensionality reduction, marker detection and plot markers for cluster annotation.
The chickent and mouse ortholog table (~/spinal_cord_paper/data/ortho_gg_mm_v102.rds) is retrieved from [BioMart](http://nov2020.archive.ensembl.org/biomart/martview/) (release 102).
-   ***Seurat_integration_array.sh***: Integration step of the 4 integrated data sets
    -   Devel: B05, B07, and B10<sub>1</sub> (day 5, 7, and 10) brachial samples
    -   brachial: B10<sub>1</sub> and B10<sub>2</sub>
    -   lumbar: L10<sub>1</sub> and L10<sub>2</sub>
    -   polydactyl: P10<sub>1</sub> and P10<sub>2</sub>
-   ***Seurat_Gg_NT_int_array.sh***: Create and analyse the integrated seurat objects, run dimensionality reduction, marker detection and plot markers for cluster annotation.
-   ***Seurat_Gg_ctrl_lumb_ctrl_poly_integration_array.sh***: Integration and analysis of the brachial & lumbar (B/L10<sub>int</sub>: B10<sub>1</sub>, B10<sub>2</sub>, L10<sub>1</sub>, & L10<sub>2</sub>), and brachial & poly (B/P10<sub>int</sub>:  B10<sub>1</sub>, B10<sub>2</sub>, P10<sub>1</sub>, & P10<sub>2</sub>). They get an individual, combined pipeline since the time and memory allocation are higher.
-   ***Seurat_Gg_all_int_render.sh***: Integration and analysis of all 8 samples for Supplementary figure 1 gets it own job due to time and memory demands.

Cluster annotation is guided by marker gene expression:

-   ***Broad_clusters_and_DV_domain_plots_render.sh***: Plots marker gene expression of all data sets and exports them to annotations/figures. Slurmed by *annotations/Broad_clusters_and_DV_domain_plots_render.Rmd*.

## scWGNCA

Module construction and co-expression analysis is done with [scWGCNA](https://github.com/CFeregrino/scWGCNA). The markdown scripts are again called with and array job:

-   ***scWGCNA_array.sh***: Run scWGCNA on the integrated ctrl, lumb, poly, and devel data sets. Produces the Gg_*_int_scWGCNA_modules and _modules_expression.pdf plots.

Next we run the comparative scWGCNA pipeline to calculate module conservation between the data sets:

-   ***comparative_scWGCNA_array.sh***: compares each data set's modules with the other two data sets. Produces the "Gg_*_int_comp_scWGCNA.pdf" plots.

## miloR (differential abundance)

Differential abundance is calculated with [miloR](https://github.com/MarioniLab/miloR). Again the script is called with a  separate array job:

-   ***miloR_Gg_ctrl_lumb_cltr_poly_int_array.sh***: Runs miloR to calculate differential abundance.

## motor neuron counting

Preparing the manual count data for the motor neuron counting, line plots and bar plots (Figure 5).

- ***MN_counting.Rmd***

## Plotting

The following scripts are used to create the plots for the figures, if they are not already produced by the pipelines above:

-   ***dotplots_broad_render.sh***: Plots the dotplots that show average marker gene expression by cluster.

-   ***tsne_and_bar_plots_render.sh***: Plots the tsne embeddings colored by broad clustering as well as stacked bar plot by data set to show the relative cell type contributions.

-   ***heatmap_spearman_ctrl_lumb_poly_int.sh*** and ***heatmap_spearman_devel.sh***: Produces the heatmaps of cluster and module correlation.

-   ***Vln_plots_render.sh***: Plots diagnostic violin/beeswarm plots to quickly reference mitochondrial and ribosomal gene expression by cluster.

-   ***Gg_devel_scWGCNA_module_analysis_render.sh***: All the analysis of the scWGCNA modules on the development (D5, D7, and D10/ctrl1) data set.

-   ***Fig_1_plots.R***: Produces the PCA of the Gg_all_int pseudobulk data shown in Supp. Fig. 1D.
-   ***Fig_2_plots.R***: Produces the module network plots of Fig. 2A-C.
-   ***Fig_3_plots.R***: Produces the additional plots for Figure 3. 
-   ***Fig_4_plots.R***: Produces the additional plots for Figure 4 and Supp. Fig 4.
-   ***Fig_5_plots.R***: Produces the additional plots for Figure 5 and Supp. Fig 5.
