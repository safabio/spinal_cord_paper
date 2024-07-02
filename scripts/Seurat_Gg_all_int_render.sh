#!/bin/bash


#SBATCH --job-name=Seurat_all_integration
#SBATCH --time=02:00:00
#SBATCH --qos=6hours
#SBATCH --mem=80G
#SBATCH --mail-user=fabio.sacher@unibas.ch
#SBATCH --mail-type=ALL
#SBATCH --output ./slurm_out/%x-%j.out


module purge
ml load R/4.1.0-foss-2018b
ml Pandoc/2.7.3

Rscript -e "rmarkdown::render('~/spinal_cord_paper/markdown/Gg_all_int.Rmd')"
