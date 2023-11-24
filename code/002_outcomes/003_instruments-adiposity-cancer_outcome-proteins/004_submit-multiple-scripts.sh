#!/bin/bash

#SBATCH --job-name=make-multiple-submissionscripts
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M
#SBATCH --partition=low_p

# set directory ====
cd ~/_PROJECTS/adiposity-proteins-endometrial/code/002_outcomes/003_instruments-adiposity-cancer_outcome-proteins

# submit the filelist-*.sh scripts  ====
for file in filelist-*.sh; do
sbatch $file
sleep 2
done
