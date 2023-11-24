#!/bin/bash

#SBATCH --job-name=extract-outcome-proteins_master
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M
#SBATCH --partition=low_p

cd /home/wangs/_PROJECTS/adiposity-proteins-endometrial/data/filelist/outcome/filelist_proteins
ls /data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS/8942_2_MRPL21_RM21.txt.gz.annotated.gz.exclusions.gz.alleles.gz /data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS/9003_99_MARCO_MARCO.txt.gz.annotated.gz.exclusions.gz.alleles.gz /data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS/9380_2_PLA2G12B_sPLA_2__XIII.txt.gz.annotated.gz.exclusions.gz.alleles.gz > filelist-discrepencies
wc -l filelist-discrepencies

# set directory ====
export SNP_LIST=/home/wangs/_PROJECTS/adiposity-proteins-endometrial/data/filelist/instruments/instruments_adiposity-cancer-combined.txt
export OUTCOMES=/home/wangs/_PROJECTS/adiposity-proteins-endometrial/data/filelist/outcome/
  
  cd /home/wangs/_PROJECTS/adiposity-proteins-endometrial/data/filelist/outcome/filelist_proteins
tmp=$(mktemp) || { ret="$?"; printf 'Failed to create temp file\n'; exit "$ret"; }

VAR1=filelist-discrepencies
  
  for file in $(cat "$VAR1"); do
echo "running $(basename "$file")"
zgrep -w -F -f "$SNP_LIST" "$file" > "$tmp" &&
  mv -- "$tmp" "${OUTCOMES}$(basename "$file")"
echo "finished $(basename "$file")"
done
