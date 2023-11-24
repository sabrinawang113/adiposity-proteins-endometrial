rm(list=ls())
set.seed(821)

# environment ====
library(TwoSampleMR)
library(ieugwasr)
library(genetics.binaRies)
library(data.table)

# data ====
list_files <- list.files("/data/GWAS_data/work/omara_2018_PMID30093612/", pattern = "gz", full.names = T)
list_files <- list_files[grepl("ECAC2018", list_files)]
list_files <- list_files[!grepl("nonendometrioid", list_files)] ## no instrument with p < 5e-8, run separately with p < 5e-7

# select instruments ====
for (FILE in list_files){
  # make label
  label <- gsub("/data/GWAS_data/work/omara_2018_PMID30093612//", "", FILE)
  label <- gsub(".txt.gz", "", label)
  # data
  data <- read_exposure_data(filename = FILE, sep = "\t", 
                             clump = F, 
                             phenotype_col = "phenotype", 
                             snp_col = "SNP", 
                             beta_col = "BETA", 
                             se_col = "SE", 
                             pval_col = "P",
                             eaf_col = "EAF", 
                             effect_allele_col = "EA", 
                             other_allele_col = "OA", 
                             chr_col = "CHR", 
                             pos_col = "POS",
                             samplesize_col = "N")
  data$id.exposure <- label
  colnames(data)[colnames(data) == "SNP"] <- "rsid"
  colnames(data)[colnames(data) == "pval.exposure"] <- "pval"
  ## clump
  data_clumped <- ld_clump(dat = data,
                           clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8,
                           pop = "EUR",
                           access_token = NULL,
                           bfile = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR",
                           plink_bin = genetics.binaRies::get_plink_binary())
  colnames(data_clumped)[colnames(data_clumped) == "rsid"] <- "SNP"
  colnames(data_clumped)[colnames(data_clumped) == "pval"] <- "pval.exposure"
  ## save
  nrow(data)
  nrow(data_clumped)
  write.table(data_clumped, paste0("analysis/001_instruments/cancer/instruments_", label, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}
  
# combine ====
list_files <- list.files("analysis/001_instruments/cancer/", pattern = "txt", full.names = T)
list_data <- lapply(list_files, fread)
data <- bind_rows(list_data, .id = "id")
write.table(data, "analysis/001_instruments/cancer/instruments_cancer.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
