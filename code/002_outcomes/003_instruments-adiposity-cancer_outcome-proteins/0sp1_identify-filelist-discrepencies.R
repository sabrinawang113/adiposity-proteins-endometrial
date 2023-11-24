## identify discrepancies between list_original-proteins and list_extracted-proteins ======

# cr list_extracted-proteins ======
rm(list=ls())
list_proteins <- data.frame()
list_files <- list.files("~/_PROJECTS/adiposity-proteins-endometrial/data/filelist/outcome", pattern = "gz", full.names = T)
for (FILE in list_files){
label <- gsub("/home/wangs/_PROJECTS/adiposity-proteins-endometrial/data/filelist/outcome/", "", FILE)
label <- gsub(".txt.gz.annotated.gz.exclusions.gz.alleles.gz", "", label)
list_proteins <- rbind(list_proteins, label)
}
colnames(list_proteins)[colnames(list_proteins) == "X.10000_28_CRYBB2_CRBB2."] <- "V1"


# cr list_original-proteins ======
list_ukb <- data.frame()
list_files_ukb <- list.files("/data/GWAS_data/work/UKB_PPP/european", pattern = "gz", full.names = T)
#FILE <- list_files_ukb[[1]]
for (FILE in list_files_ukb){
  label <- gsub("/data/GWAS_data/work/UKB_PPP/european/", "", FILE)
  #label <- gsub(".gz", "", label)
  list_ukb <- rbind(list_ukb, label)
}
colnames(list_ukb)[colnames(list_ukb) == "X.A1BG_P04217_OID30771.gz."] <- "V1"

list_ferkingstad <- data.frame()
list_files_ferkingstad <- list.files("/data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS/", pattern = "gz", full.names = T)
FILE <- list_files_ferkingstad[[1]]
for (FILE in list_files_ferkingstad){
  label <- gsub("/data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS//", "", FILE)
  label <- gsub(".txt.gz.annotated.gz.exclusions.gz.alleles.gz", "", label)
  list_ferkingstad <- rbind(list_ferkingstad, label)
}
colnames(list_ferkingstad)[colnames(list_ferkingstad) == "X.10000_28_CRYBB2_CRBB2."] <- "V1"
list_proteins_all <- rbind(list_ukb, list_ferkingstad)


# find discrepencies ====== 
list_diff <- setdiff(list_proteins_all, list_proteins)
# save ====
write.table(list_diff, "data/filelist/discrepencies.txt",
            row.names = FALSE, col.names = T, quote = FALSE, sep = "\t")

# check proteins not extracted =====
#data <- fread("/data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS/5008_51_SOD2_Mn_SOD.txt.gz.annotated.gz.exclusions.gz.alleles.gz")
#data <- fread("/data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS/5715_4_ABHD14A_ABHEA.txt.gz.annotated.gz.exclusions.gz.alleles.gz")
#data <- fread("/data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS/8942_2_MRPL21_RM21.txt.gz.annotated.gz.exclusions.gz.alleles.gz")

#data <- fread("/data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS/10000_28_CRYBB2_CRBB2.txt.gz.annotated.gz.exclusions.gz.alleles.gz")
