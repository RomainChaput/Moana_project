#!/bin/bash -e
#SBATCH --job-name=ArrayJob_Nine_Beach_forward_1_run_June    #run for June 2017 Ninety Miles Beach
#SBATCH --account=vuw03295                          #Project number
#SBATCH --time=16:00:00                             #10 hours was too short
#SBATCH --mem=64000                                 #40Gb of memory gave me a out of memory error (49 sites x 10 pts x 250 particles)
#SBATCH --array=1-4

module purge
module load Miniconda3
source activate opendrift

python ./inputs_June/mussel_90milebeach_foreward_1_june_${SLURM_ARRAY_TASK_ID}.py
