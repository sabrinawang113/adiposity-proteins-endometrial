#!/bin/bash

#SBATCH --job-name=001_MR_proteins-cancer
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M
#SBATCH --partition=low_p

cd /home/wangs/_PROJECTS/adiposity-proteins-endometrial
Rscript code/007_MR_proteins-cancer/001_MR_proteins-cancer.R