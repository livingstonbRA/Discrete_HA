#!/bin/bash
#SBATCH --job-name=disc
#SBATCH --output=/home/livingstonb/GitHub/Discrete_HA/output/run%a.out
#SBATCH --error=/home/livingstonb/GitHub/Discrete_HA/output/run%a.err
#SBATCH --partition=broadwl
#SBATCH --time=12:00:00
#SBATCH --array=1-49
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10000

module load matlab/2019b

matlab -nodisplay < /home/livingstonb/GitHub/Discrete_HA/master.m
