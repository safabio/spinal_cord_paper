#!/bin/bash


#SBATCH --job-name=DAseq_Gg_ctrl_lumb_ctrl_poly_int_array
#SBATCH --time=01:00:00
#SBATCH --qos=6hours
#SBATCH --mem=15G
#SBATCH --mail-user=fabio.sacher@unibas.ch
#SBATCH --mail-type=ALL
#SBATCH --output ./slurm_out/%x-%A_%a.out
#SBATCH --array=1-2


module purge
ml load R/4.1.0-foss-2018b
ml Pandoc/2.7.3

$(head -$SLURM_ARRAY_TASK_ID DAseq_Gg_ctrl_lumb_ctrl_poly_int.cmd | tail -1) 