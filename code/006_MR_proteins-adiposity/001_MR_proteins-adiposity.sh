#!/bin/bash

#SBATCH --job-name=001_MR_proteins-adiposity
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M
#SBATCH --partition=low_p

cd /home/wangs/_PROJECTS/adiposity-proteins-endometrial
Rscript code/006_MR_proteins-adiposity/001_MR_proteins-adiposity.R