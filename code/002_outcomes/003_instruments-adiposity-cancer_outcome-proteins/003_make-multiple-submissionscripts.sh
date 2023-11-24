#!/bin/bash

#SBATCH --job-name=make-multiple-submissionscripts
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M
#SBATCH --partition=low_p

# set directory ====
cd ~/_PROJECTS/adiposity-proteins-endometrial/code/002_outcomes/003_instruments-adiposity-cancer_outcome-proteins
rm filelist*
rm filenames
  
# create a file with the names of all files within a directory  ====
ls ~/_PROJECTS/adiposity-proteins-endometrial/data/filelist/outcome/filelist_proteins > filenames.txt 

# create multiple .sh scripts from a single file with names in based on a master script  ====
# line (NR == 17) holds the variable name
# line (NR == 3) holds the job name
cat filenames.txt | while read i; do echo ${i}; 
awk '{ 
  if (NR == 3) 
    print "#SBATCH --job-name=extract-outcome-proteins_'${i}'";
  else if (NR == 17)
    print "VAR1='${i}'";
  else
    print $0
}' 002_master-filelist.sh > ${i}.sh; done 
rm filenames.txt