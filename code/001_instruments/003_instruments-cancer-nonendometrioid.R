rm(list=ls())
set.seed(821)

# environment ====
library(TwoSampleMR)
library(ieugwasr)
library(genetics.binaRies)
library(data.table)

## manually prepare exposure data for non-endometrioid GWAS ======
## no instrument with p < 5e-8, run separately with p < 5e-7

# data ====
list_files <- list.files("/data/GWAS_data/work/omara_2018_PMID30093612/", pattern = "gz", full.names = T)
list_files <- list_files[grepl("ECAC2018", list_files)]
FILE <- list_files[[2]]
data1 <- fread(FILE)

# set min pval 5e-7 ====
min_value <- min(data1$P, na.rm = TRUE)
max_value <- max(data1$P, na.rm = TRUE)
data <- subset(data1, P < 5e-7)

# format exposure data ====
# make label
label <- gsub("/data/GWAS_data/work/omara_2018_PMID30093612//", "", FILE)
label <- gsub(".txt.gz", "", label)
# format
data <- format_data(data, type="exposure", snps = NULL,  header = TRUE,  
                    phenotype_col = "phenotype", 
                    id_col = "phenotype",
                    snp_col = "SNP", 
                    beta_col = "BETA", 
                    se_col = "SE", 
                    pval_col = "P",
                    eaf_col = "EAF", 
                    effect_allele_col = "EA", 
                    other_allele_col = "OA", 
                    chr_col = "CHR", 
                    pos_col = "POS",  
                    log_pval = FALSE)
data$id.exposure <- label
# prep for clump
colnames(data)[colnames(data) == "SNP"] <- "rsid"
colnames(data)[colnames(data) == "pval.exposure"] <- "pval"

# clump ====
data_clumped <- ld_clump(dat = data,
                         clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-7,
                         pop = "EUR",
                         access_token = NULL,
                         bfile = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR",
                         plink_bin = genetics.binaRies::get_plink_binary())
colnames(data_clumped)[colnames(data_clumped) == "rsid"] <- "SNP"
colnames(data_clumped)[colnames(data_clumped) == "pval"] <- "pval.exposure"

## save ====
nrow(data)
nrow(data_clumped)
write.table(data_clumped, paste0("analysis/001_instruments/cancer/instruments_", label, ".txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
