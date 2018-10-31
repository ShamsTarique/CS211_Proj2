#!/bin/bash
#SBATCH --job-name=part2
   #SBATCH --output=part2.txt
   #SBATCH -N 1
   #SBATCH -n 32
   #SBATCH -t 01:00:00
   
   ./part2

