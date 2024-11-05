# ./scripts

## *_.cmd files

Contains the commands for the array jobs to render the markdown files and are slurmed with the corresponding *_array.sh bash files.

## *_render.sh

Submits the corresponding .Rmd (Rmarkdown file) to slurm to render it.

## Fig_*_plots.R

Additional plots and code to inspect data for the individual figures.Those plots are not incorporated into the markdown to prevent rerunning big parts of the pipelines.

## marker_tsnes.R

Feature tSNE plots of different marker genes on all data sets for any questions and ideas during the manuskript writing.

## REVIGO_Gg_devel_modules.R

Script for the semantic GO term analysis. Loads and prepares the data for the upload to [REVIGO](http://revigo.irb.hr/), contains the results of the analysis, and plots them.

## source code

  ***scWGCNA_modplots.R***: Contains modified code from scWGCNA to plot the comparative scWGCNA outputs.
 
## Supplementary_tables.R

Script to prepare and export the supplementary tables for the manuscript.


