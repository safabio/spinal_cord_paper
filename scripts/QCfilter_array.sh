#!/bin/bash


#SBATCH --job-name=QC_filter_array
#SBATCH --time=00:10:00
#SBATCH --qos=30min
#SBATCH --mem=10G
#SBATCH --mail-user=fabio.sacher@unibas.ch
#SBATCH --mail-type=ALL
#SBATCH --output ./slurm_out/%x-%A_%a.out
#SBATCH --array=1-8


module purge
ml load R/4.1.0-foss-2018b
ml Pandoc/2.7.3

$(head -$SLURM_ARRAY_TASK_ID QCfilter.cmd | tail -1) 
