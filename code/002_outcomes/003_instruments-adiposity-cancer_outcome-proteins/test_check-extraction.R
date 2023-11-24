
rm(list=ls())

data_outcome <- fread("data/filelist/outcome/A1BG_P04217_OID30771.gz")
data_outcome <- fread("data/filelist/outcome/10000_28_CRYBB2_CRBB2.txt.gz.annotated.gz.exclusions.gz.alleles.gz")
data_instruments <- fread("data/filelist/instruments/instruments_adiposity-cancer-combined.txt", header = F)

outcome <- data_outcome$V4
instrument <- data_instruments$V1

snp_common <- intersect(instrument, outcome)
snp_diff <- setdiff(instrument, outcome)

