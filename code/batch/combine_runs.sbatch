#!/bin/bash
#SBATCH --job-name=disc_combine_runs
#SBATCH --output=/home/livingstonb/GitHub/Discrete_HA/output/combine_runs.out
#SBATCH --error=/home/livingstonb/GitHub/Discrete_HA/output/combine_runs.err
#SBATCH --partition=broadwl
#SBATCH --time=00:20:00
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=8000

module load matlab/2019b

matlab -nodisplay < "/home/livingstonb/GitHub/Discrete_HA/code/+statistics/combine_runs.m"
