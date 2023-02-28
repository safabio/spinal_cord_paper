#!/bin/bash


#SBATCH --job-name=Broad_clusters_and_DV_domain_plots
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --mem=35G
#SBATCH --mail-user=fabio.sacher@unibas.ch
#SBATCH --mail-type=ALL
#SBATCH --output ./slurm_out/%x.out


module purge
ml load R/4.1.0-foss-2018b
ml Pandoc/2.7.3

Rscript -e "rmarkdown::render('~/spinal_cord_paper/annotations/Broad_clusters_and_DV_domain_plots.Rmd')"