#!/bin/bash


#SBATCH --job-name=heatmap_spearman_ctrl_lumb_poly_int_render
#SBATCH --time=00:02:00
#SBATCH --qos=30min
#SBATCH --mem=10G
#SBATCH --mail-user=fabio.sacher@unibas.ch
#SBATCH --mail-type=ALL
#SBATCH --output ./slurm_out/%x-%j.out


module purge
ml load R/4.1.0-foss-2018b
ml Pandoc/2.7.3

Rscript -e "rmarkdown::render('~/spinal_cord_paper/markdown/heatmap_spearman_ctrl_lumb_poly_int.Rmd')"