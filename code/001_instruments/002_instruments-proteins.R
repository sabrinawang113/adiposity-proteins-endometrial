rm(list=ls())
set.seed(821)

# environment ====
library(TwoSampleMR)
library(data.table)
library(dplyr)

# data : OLINK ====
list_files <- list.files("/data/GWAS_data/work/UKB_PPP/cis/", ".txt", all.files = T, full.names = T)
for (FILE in list_files){
  # make label
  label <- gsub("/data/GWAS_data/work/UKB_PPP/cis//", "", FILE)
  label <- gsub(".txt", "", label)
  # data
  data <- read_exposure_data(filename = FILE, sep = "\t", 
                             clump = F, 
                             min_pval = 5e-08,
                             phenotype_col = "phenotype", 
                             id_col = "phenotype", 
                             snp_col = "SNP", 
                             beta_col = "BETA", 
                             se_col = "SE", 
                             pval_col = "P", log_pval = TRUE,
                             eaf_col = "EAF", 
                             effect_allele_col = "EA", 
                             other_allele_col = "OA", 
                             chr_col = "CHR", 
                             pos_col = "POS",
                             samplesize_col = "N")
  data$id.exposure <- paste0(data$exposure, "_", label)
  ## save
  write.table(data, paste0("analysis/001_instruments/proteins/instruments_", label, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

# data : SOMALOGIC ====
data <- fread("/data/GWAS_data/work/ferkingstad_2021_PMID34857953/cis_snps.txt")
data <- format_data(data, type="exposure", snps = NULL,  header = TRUE,  
                    phenotype_col = "phenotype",  
                    snp_col = "rsID",  
                    beta_col = "BETA",  
                    se_col = "SE",  
                    eaf_col = "EAF",  
                    effect_allele_col = "EffectAllele",  
                    other_allele_col = "OtherAllele",  
                    pval_col = "Pval", 
                    id_col = "ID", 
                    chr_col = "CHR",  
                    pos_col = "BP",  
                    log_pval = FALSE, min_pval = 1.8e-9)
data$id.exposure <- paste0(data$exposure, "_ferkingstad")
## save
write.table(data, "analysis/001_instruments/proteins/instruments_ferkingstad.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# combine ====
list_files <- list.files("analysis/001_instruments/proteins/", pattern = "txt", full.names = T)
list_data <- lapply(list_files, fread)
data <- bind_rows(list_data, .id = "id")
write.table(data, "analysis/001_instruments/proteins/instruments_proteins_all.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# combine UKB_discovery-EC & ferkingstad only ====
list_files <- list.files("analysis/001_instruments/proteins/", pattern = "txt", full.names = T)
list_files <- list_files[!grepl("combined", list_files)]
list_files <- list_files[!grepl("nonEU", list_files)]
list_files <- list_files[!grepl("proteins_all", list_files)]

list_data <- lapply(list_files, fread)
data <- bind_rows(list_data, .id = "id")
write.table(data, "analysis/001_instruments/proteins/instruments_proteins.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
