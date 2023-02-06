# Spinal cord paper

This repository contains the code to recreate the findings in our current publication " " xxx.

The Initial QC filtering and seurat object creation markdown files are rendered using an array job in the scripts/ directory. Use the following order:

-   ***QCfilter_array.sh***: initial removal of ambient RNA, quality control and filtering with SoupX and ddqcR.
-   ***Seurat_Gg_NT_array.sh***: create and analyse seurat objects for all the individual samples, run dimensionality reduction, marker detection and plot markers for cluster annotation.
-   ***Seurat_integration_array.sh***: Integration step of the 4 integrated data sets
    -   Devel: D05_ctrl, D07_ctrl, and ctrl_1 (day 10) brachial samples
    -   brachial: ctrl_1 and ctrl_2
    -   lumbar: lumb_1 and lumb_2
    -   polydactyl: poly_1 and poly_2
-   ***Seurat_Gg_NT_int_array.sh***: create and analyse the integrated seurat objects, run dimensionality reduction, markter detection and plot markers for cluster annotation.
