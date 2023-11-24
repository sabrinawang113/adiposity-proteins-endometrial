#!/bin/bash

#SBATCH --job-name=instruments-adiposity
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10-1:0:00
#SBATCH --mem=100000M
#SBATCH --partition=low_p

cd /home/wangs/_PROJECTS/adiposity-proteins-endometrial
Rscript code/001_instruments/001_instruments-adiposity.R