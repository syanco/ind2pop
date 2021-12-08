#!/bin/bash

#SBATCH -t 12:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH -c 1
#SBATHC --mem-per-cpu 50G
#SBATCH -J annotate_storks_and_elephants

cd ~/project/ind2pop

module load miniconda
conda activate niche_mix

Rscript ./src/annotate_data.r