#!/bin/bash
#SBATCH --partition=medium
#SBATCH --job-name=sim
#SBATCH --output=sim_%A.out
#SBATCH --error=sim_%A.err
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=15G
#SBATCH --array=1-13

module add apps/R/4.2.1

TAG=$1  # passed in as first argument
START=$(( ($SLURM_ARRAY_TASK_ID - 1) * 8 + 1 ))
END=$(( $SLURM_ARRAY_TASK_ID * 8 ))

Rscript Step1_estTransitivity.R $START $END $TAG
