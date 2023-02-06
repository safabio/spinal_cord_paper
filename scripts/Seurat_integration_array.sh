#!/bin/bash


#SBATCH --job-name=Seurat_integration
#SBATCH --time=00:10:00
#SBATCH --qos=30min
#SBATCH --mem=10G
#SBATCH --mail-user=fabio.sacher@unibas.ch
#SBATCH --mail-type=ALL
#SBATCH --output ./slurm_out/%x-%A_%a.out
#SBATCH --array=1-4


module purge
ml load R/4.1.0-foss-2018b
ml Pandoc/2.7.3

$(head -$SLURM_ARRAY_TASK_ID Seurat_integration.cmd | tail -1) 
