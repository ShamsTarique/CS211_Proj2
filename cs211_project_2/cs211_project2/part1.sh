#!/bin/bash
#BATCH --job-name=part1
#SBATCH --output=part1.txt
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 01:00:00

./part1
