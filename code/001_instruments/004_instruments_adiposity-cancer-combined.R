rm(list=ls())
set.seed(821)

data1 <- fread("analysis/001_instruments/adiposity/instruments_adiposity.txt")
data2 <- fread("analysis/001_instruments/cancer/instruments_cancer.txt")

data <- data_frame(c(data1$SNP, data2$SNP))
data <- data[!duplicated(data), ]


write.table(data, paste0("data/filelist/instruments/instruments_adiposity-cancer-combined.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
