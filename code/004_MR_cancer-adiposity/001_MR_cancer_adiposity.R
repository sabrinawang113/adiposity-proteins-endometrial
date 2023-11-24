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
  "mr_egger_regression",
  "mr_weighted_median",
  "mr_ivw_mre",
  "mr_weighted_mode")

## exposure ====
data_exposure <- fread("analysis/001_instruments/cancer/instruments_cancer.txt")
data_exposure$f_stats <- (data_exposure$b / data_exposure$se)^2 
data_exposure %>%
  group_by(id.exposure) %>%
  summarise(mean = mean(f_stats))
write.table(data_exposure, "analysis/004_MR_cancer-adiposity/data_exposure.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## outcome ====
data_outcome <- fread("analysis/002_outcomes/adiposity/instruments-cancer_outcome-adiposity.txt")
write.table(data_outcome, "analysis/004_MR_cancer-adiposity/data_outcome.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## harmonise ====
data_harmonise <- harmonise_data(data_exposure, data_outcome, action=2)
write.table(data_harmonise, "analysis/004_MR_cancer-adiposity/data_harmonise.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## MR ====
data_mr <- mr(data_harmonise, method_list = methods)
write.table(data_mr, "analysis/004_MR_cancer-adiposity/data_mr.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## sensitivity analyses ====
data_sensitivity <- mr_singlesnp(data_harmonise)
write.table(data_sensitivity, "analysis/004_MR_cancer-adiposity/data_mr-singlesnp.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
plot_sensitivity <- mr_forest_plot(data_sensitivity)
            pdf("analysis/004_MR_cancer-adiposity/figures/plot_mr-singlesnp.pdf")
            for (i in 1:length(plot_sensitivity)) {
              print(plot_sensitivity[[i]])
            }
            dev.off()
plot_sensitivity <- mr_funnel_plot(data_sensitivity)
            pdf("analysis/004_MR_cancer-adiposity/figures/plot_mr-funnelplot.pdf")
            for (i in 1:length(plot_sensitivity)) {
              print(plot_sensitivity[[i]])
            }
            dev.off()

data_sensitivity <- mr_heterogeneity(data_harmonise, method_list = methods_heterogeneity)
write.table(data_sensitivity, "analysis/004_MR_cancer-adiposity/data_mr-heterogeneity.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data_sensitivity <- mr_pleiotropy_test(data_harmonise)
write.table(data_sensitivity, "analysis/004_MR_cancer-adiposity/data_mr-pleiotropy.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data_sensitivity <- mr_leaveoneout(data_harmonise)
write.table(data_sensitivity, "analysis/004_MR_cancer-adiposity/data_mr-leaveoneout.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
plot_sensitivity <- mr_leaveoneout_plot(data_sensitivity)
            pdf("analysis/004_MR_cancer-adiposity/figures/plot_mr-leaveoneout.pdf")
            for (i in 1:length(plot_sensitivity)) {
              print(plot_sensitivity[[i]])
            }
            dev.off()

plot_sensitivity <- functions::mr_scatter_plot(data_mr, data_harmonise)
            pdf("analysis/004_MR_cancer-adiposity/figures/plot_mr-scatter.pdf")
            for (i in 1:length(plot_sensitivity)) {
              print(plot_sensitivity[[i]])
            }
            dev.off()

