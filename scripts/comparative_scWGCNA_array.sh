#!/bin/bash


#SBATCH --job-name=comparative_scWGCNA_ctrl_lumb_poly_int
#SBATCH --time=00:20:00
#SBATCH --qos=30min
#SBATCH --mem=20G
#SBATCH --mail-user=fabio.sacher@unibas.ch
#SBATCH --mail-type=ALL
#SBATCH --output ./slurm_out/%x-%A_%a.out
#SBATCH --array=1-3


module purge
ml load R/4.1.0-foss-2018b
ml Pandoc/2.7.3

$(head -$SLURM_ARRAY_TASK_ID comparative_scWGCNA.cmd | tail -1) 
