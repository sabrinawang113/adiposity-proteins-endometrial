#!/bin/bash

#SBATCH --job-name=extract-outcome-proteins_filelist-47
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M
#SBATCH --partition=low_p

# set directory ====
export SNP_LIST=/home/wangs/_PROJECTS/adiposity-proteins-endometrial/data/filelist/instruments/instruments_adiposity-cancer-combined.txt
export OUTCOMES=/home/wangs/_PROJECTS/adiposity-proteins-endometrial/data/filelist/outcome/
  
cd /home/wangs/_PROJECTS/adiposity-proteins-endometrial/data/filelist/outcome/filelist_proteins
tmp=$(mktemp) || { ret="$?"; printf 'Failed to create temp file\n'; exit "$ret"; }

VAR1=filelist-47
  
for file in $(cat "$VAR1"); do
  echo "running $(basename "$file")"
  zgrep -w -F -f "$SNP_LIST" "$file" > "$tmp" &&
  mv -- "$tmp" "${OUTCOMES}$(basename "$file")"
  echo "finished $(basename "$file")"
done
