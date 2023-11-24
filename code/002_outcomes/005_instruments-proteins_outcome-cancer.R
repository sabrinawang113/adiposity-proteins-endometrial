rm(list=ls())
set.seed(821)

# environment ====
library(TwoSampleMR)
library(ieugwasr)
library(genetics.binaRies)
library(data.table)

# data ====
data <- fread("analysis/001_instruments/proteins/instruments_proteins.txt")
list_files <- list.files("/data/GWAS_data/work/omara_2018_PMID30093612/", pattern = "ECAC", full.names = T)
list_data <- list()

# extract adiposity SNPs from cancer GWAS
for (FILE in list_files){
  data_outcome <- read_outcome_data(filename = FILE, 
                                    sep = "\t",
                                    snps = data$SNP,
                                    phenotype_col = "phenotype", 
                                    snp_col = "SNP", 
                                    beta_col = "BETA", 
                                    se_col = "SE", 
                                    eaf_col = "EAF", 
                                    effect_allele_col = "EA", 
                                    other_allele_col = "OA", 
                                    pval_col = "PVALUE", 
                                    chr_col = "CHR", 
                                    pos_col = "POS")
  data_outcome$id.outcome <- data_outcome$outcome
  
  list_data <- append(list_data, list(data_outcome))
}
data_all <- bind_rows(list_data)

# save ====
write.table(data_all, "analysis/002_outcomes/cancer/instruments-proteins_outcome-cancer.txt",
            row.names = FALSE, col.names = T, quote = FALSE, sep = "\t")

