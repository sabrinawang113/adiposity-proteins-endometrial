#!/bin/bash

#SBATCH --job-name=instruments-adiposity
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M
#SBATCH --partition=low_p

cd /home/wangs/_PROJECTS/adiposity-proteins-endometrial
Rscript code/009_colocalization/003_coloc_proteins-ukb_ECAC-noUKBB.R