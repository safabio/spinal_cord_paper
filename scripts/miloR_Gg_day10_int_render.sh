#!/bin/bash


#SBATCH --job-name=Milo_Gg_day10_int
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --mem=10G
#SBATCH --mail-user=fabio.sacher@unibas.ch
#SBATCH --mail-type=ALL
#SBATCH --output ./slurm_out/%x-%j.out


module purge
ml load R/4.4.1-foss-2023b
ml Pandoc/2.13

Rscript -e "rmarkdown::render('~/spinal_cord_paper/markdown/miloR_Gg_day10_int.Rmd')"