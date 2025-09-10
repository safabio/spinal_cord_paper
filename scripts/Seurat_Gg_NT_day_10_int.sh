#!/bin/bash


#SBATCH --job-name=Seurat_Gg_NT_day10
#SBATCH --time=03:00:00
#SBATCH --qos=6hours
#SBATCH --mem=40G
#SBATCH --mail-user=fabio.sacher@unibas.ch
#SBATCH --mail-type=ALL
#SBATCH --output ./slurm_out/%x-%j.out


module purge
ml load R/4.4.1-foss-2023b
ml Pandoc/2.13

Rscript -e "rmarkdown::render('~/spinal_cord_paper/markdown/Gg_day_10_int.Rmd')"