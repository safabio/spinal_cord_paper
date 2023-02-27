# Spinal cord paper

This repository contains the code to recreate the findings in our current publication " " xxx.

We calculate the tSNE embeddings with the FTF (Fast Fourier Transform) accerelated Interpolation-based t-SNE [FIt-SNE](https://github.com/KlugerLab/FIt-SNE) developed by the Kluger lab. Install it to a directory to your machine following their instruction.

The Initial QC filtering and seurat object creation markdown files are rendered using an array job in the scripts/ directory. Use the following order:

-   ***QCfilter_array.sh***: initial removal of ambient RNA, quality control and filtering with SoupX and ddqcR.
-   ***Seurat_Gg_NT_array.sh***: create and analyse seurat objects for all the individual samples, run dimensionality reduction, marker detection and plot markers for cluster annotation.
-   ***Seurat_integration_array.sh***: Integration step of the 4 integrated data sets
    -   Devel: D05_ctrl, D07_ctrl, and ctrl_1 (day 10) brachial samples
    -   brachial: ctrl_1 and ctrl_2
    -   lumbar: lumb_1 and lumb_2
    -   polydactyl: poly_1 and poly_2
-   ***Seurat_Gg_NT_int_array.sh***: Create and analyse the integrated seurat objects, run dimensionality reduction, markter detection and plot markers for cluster annotation.
-   ***Seurat_ctrl_lumb_ctrl_poly_integration_array.sh***: Integration and analysis of the brachial & lumbar (ctrl_1, ctrl_2, lumb_1, & lumb_2), and brachial & poly (ctrl_1, ctrl_2, poly_1, & poly_2). They get an individual, combined pipeline since the time and memory allocation are higher.

Module construction and co-expression analysis is done with [scWGCNA](https://github.com/CFeregrino/scWGCNA). The markdown scripts are again called with and array job:

-   ***scWGCNA_array.sh***: Run scWGCNA on the integrated ctrl, lumb, poly, and devel data sets.

Differential abundance is calculated with [DAseq](https://github.com/KlugerLab/DAseq) and [miloR](https://github.com/MarioniLab/miloR). Again the scripts are called with separate array jobs:

-   ***DAseq_Gg_ctrl_lumb_cltr_poly_int_array.sh***: Runs DAseq to calculate differential abundance. Requires a python3 verison available (install with anaconda).
-   ***miloR_Gg_ctrl_lumb_cltr_poly_int_array.sh***: Runs miloR to calculate differential abundance.
