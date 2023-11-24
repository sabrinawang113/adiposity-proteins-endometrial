rm(list=ls())

# environment ====
#if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
#remotes::install_github("MRCIEU/genetics.binaRies", force = F)
#remotes::install_github("explodecomputer/plinkbinr", force = F)
#remotes::install_github("chr1swallace/coloc@main", force = F)
#remotes::install_github("sjmgarnier/viridis", force = F)
library(genetics.binaRies)
library(plinkbinr)
library(coloc)
library(viridis)
library(data.table)
library(ieugwasr)
library(dplyr)
library(TwoSampleMR)
library(tidyverse)
#detach("package:functions", unload = TRUE)
source("/home/wangs/_PROJECTS/adiposity-proteins-endometrial/tools/my_coloc_chriswallace.R")

# file paths ====
DIRECTORY_OUT <- "/home/wangs/_PROJECTS/adiposity-proteins-endometrial/analysis/009_colocalization/ECAC-noUKBB/"

# data_exposure ====
list_files_exposure <- list.files("/data/GWAS_data/work/UKB_PPP/cis-snps/european/", full.names = T)
list_exposure <- lapply(list_files_exposure, fread, sep = " ", header = F)
## identify missing
length(list_exposure)
list_exposure_missing <- list_exposure[sapply(list_exposure, nrow) == 0] # identify proteins with missing - probably need to extract window manually
list_exposure <- list_exposure[sapply(list_exposure, nrow) > 0] # remove proteins with no window
length(list_exposure)
# list_exposure <- purrr::discard(list_exposure, ~any(.x$CHR == "chrX")) # X CHR not available in outcome
length(list_exposure)
## label_exposure
label_exposure <- sub("/data/GWAS_data/work/UKB_PPP/cis-snps/european//", "", list_files_exposure)
label_exposure <- sub(".gz", "", label_exposure)
## format ====
print("starting: for (i in seq_along(list_exposure))")
for (i in seq_along(list_exposure)) {
  list_exposure[[i]] <- separate(list_exposure[[i]], V16, into = c("ID", "REF", "ALT", "rsid", "POS19POS38"), sep = "\t")
  colnames(list_exposure[[i]]) <- c("CHR", "POS", "SNPID", "other_allele.exposure", "effect_allele.exposure",
                                    "eaf.exposure", "INFO", "samplesize.exposure", "TEST", 
                                    "beta.exposure", "se.exposure", "CHISQ", "pval.exposure",
                                    "EXTRA", "exposure", "ID", "REF", "ALT", "SNP", "POS19POS38")
  list_exposure[[i]]$id.exposure <- label_exposure[i]}
## save
exposure <- bind_rows(list_exposure)
write.table(exposure, paste0(DIRECTORY_OUT, "data_exposure.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# data_outcome ====
list_files_outcome <- list.files("/data/GWAS_data/work/omara_2018_PMID30093612/", pattern = "ECAC", full.names = T)
list_files_outcome <- list_files_outcome[3]
## label_outcome
label_outcome <- sub("/data/GWAS_data/work/omara_2018_PMID30093612//", "", list_files_outcome)
label_outcome <- sub(".txt.gz", "", label_outcome)
## extract exposure SNPs from outcome
list_outcome <- list()
for (i in 1:length(list_exposure)){
  list_outcome[i] <- lapply(list_files_outcome, 
                            read_outcome_data,
                            snps = list_exposure[[i]]$SNP,
                            sep = "\t",
                            snp_col = "SNP",
                            beta_col = "BETA",
                            se_col = "SE",
                            eaf_col = "EAF",
                            effect_allele_col = "EA",
                            other_allele_col = "OA",
                            pval_col = "P",
                            min_pval = 1e-200,
                            log_pval = FALSE,
                            chr_col = "CHR",
                            pos_col = "POS",
                            phenotype_col = "phenotype")
  
  list_outcome[[i]]$outcome <- label_outcome
  list_outcome[[i]]$id.outcome <- paste0(label_exposure[[i]], "_", list_outcome[[i]]$outcome)}
## save
outcome <- bind_rows(list_outcome)
write.table(outcome, paste0(DIRECTORY_OUT, "data_outcome.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# harmonise ====
harmonise_data <- harmonise_data(exposure, outcome, action = 2)
## remove duplicates
harmonise_data$remove_duplicates <- paste0(harmonise_data$SNP, "_", harmonise_data$id.exposure)
harmonise_data <- harmonise_data[!duplicated(harmonise_data$remove_duplicates),]
## make a list of each harmonised data frame
list_harmonise <- split(harmonise_data, harmonise_data$id.exposure)
## save
write.table(harmonise_data, paste0(DIRECTORY_OUT, "data_harmonised.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# coloc ====
table_master <- data.frame() # make empty dataframe for final results
## LD matrix ====
i = 1
for (i in 1:length(list_harmonise)){
  # make label  
  label <- paste0(label_exposure[i], ";", label_outcome)
  
  # make ld matrix ====
  ld <- ld_matrix_local(
    list_harmonise[[i]]$SNP,
    with_alleles = FALSE, 
    bfile = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR",
    plink_bin = get_plink_exe())
  
  # format LD matrix and harmonised list ====
  ld <- ld[which(rownames(ld) %in% list_harmonise[[i]]$SNP), which(colnames(ld) %in% list_harmonise[[i]]$SNP)]
  list_harmonise[[i]] <- list_harmonise[[i]][which(list_harmonise[[i]]$SNP %in% rownames(ld)),]
  ld <- ld[match(list_harmonise[[i]]$SNP,rownames(ld)),]
  ld <- ld[,match(list_harmonise[[i]]$SNP, colnames(ld))]
  list_harmonise[[i]] <- list_harmonise[[i]][match(rownames(ld), list_harmonise[[i]]$SNP),]
  
  # make lists for coloc ====
  coloc_data_exposure <- list(beta = list_harmonise[[i]]$beta.exposure, varbeta = list_harmonise[[i]]$se.exposure^2, MAF = list_harmonise[[i]]$eaf.exposure, type = "quant", N = 35559, snp = rownames(ld), LD = ld, position = list_harmonise[[i]]$POS)
  coloc_data_outcome <- list(beta = list_harmonise[[i]]$beta.outcome, varbeta = list_harmonise[[i]]$se.outcome^2, MAF = list_harmonise[[i]]$eaf.outcome, type = "cc", N = 120328, snp = rownames(ld), LD = ld, position = list_harmonise[[i]]$POS)
  
  # coloc ====  
  coloc_results <- coloc.abf(dataset1 = coloc_data_exposure, dataset2 = coloc_data_outcome, p1 = 1E-6, p2 = 1E-6, p12 = 1E-7)
  
  pdf(paste0(DIRECTORY_OUT, "figures/", label, ".pdf"), 
      height = 10, width = 10)
  coloc_sensitivity <- my_sensitivity(coloc_results, "H4 > 0.8", 
                                      trait1_title = label_exposure[i], trait2_title = label_outcome)
  dev.off()
  
  # save ====
  saveRDS(coloc_results, paste0(DIRECTORY_OUT, "RData/", label, ".RData"))
  
  # make table ====
  table <- data.frame(
    exposure = label_exposure[i],
    outcome = label_outcome,
    id = label,
    nsnps = coloc_results["summary"][[1]][1],
    h0 = coloc_results["summary"][[1]][2],
    h1 = coloc_results["summary"][[1]][3],
    h2 = coloc_results["summary"][[1]][4],
    h3 = coloc_results["summary"][[1]][5],
    h4 = coloc_results["summary"][[1]][6])
  
  table_master <- rbind(table_master, table)
  
}

write.table(table_master, paste0(DIRECTORY_OUT, "table_master.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
