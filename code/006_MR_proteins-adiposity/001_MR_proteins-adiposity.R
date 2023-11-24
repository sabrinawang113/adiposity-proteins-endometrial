rm(list=ls())
set.seed(821)

# environment ====
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(RadialMR)
#remotes::install_github("mattlee821/functions")
#remotes::install_github("NightingaleHealth/ggforestplot")
library(functions)

## methods ====
methods <- mr_method_list()
methods_heterogeneity <- subset(methods, heterogeneity_test == TRUE)$obj
methods <- c(
  "mr_wald_ratio")

## exposure ====
data_exposure <- fread("analysis/001_instruments/proteins/instruments_proteins.txt")
data_exposure$f_stats <- (data_exposure$b / data_exposure$se)^2 
data_exposure %>%
  group_by(id.exposure) %>%
  summarise(mean = mean(f_stats))
data_exposure_na <- subset(data_exposure, is.na(data_exposure$other_allele.exposure))
length(unique(data_exposure_na$exposure))
data_exposure <- subset(data_exposure, !is.na(data_exposure$other_allele.exposure))
length(unique(data_exposure$exposure))
write.table(data_exposure, "analysis/006_MR_proteins-adiposity/data_exposure.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## outcome ====
data_outcome <- fread("analysis/002_outcomes/adiposity/instruments-proteins_outcome-adiposity.txt")
write.table(data_outcome, "analysis/006_MR_proteins-adiposity/data_outcome.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## harmonise ====
data_harmonise <- harmonise_data(data_exposure, data_outcome, action=2)
write.table(data_harmonise, "analysis/006_MR_proteins-adiposity/data_harmonise.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## MR ====
data_mr <- mr(data_harmonise, method_list = methods)
write.table(data_mr, "analysis/006_MR_proteins-adiposity/data_mr.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## sensitivity analyses ====
