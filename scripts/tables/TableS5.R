# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for Table S5 of manuscript

# load libraries
library(sigProcessor)
library(tidyverse)

# set working directory
setwd("~/doxorubicin_mutSig/data/ID30/")

# load in matrices
data <- process_mutation_matrix("./doxorubicin.ID30.all")

# change order of samples to match table S3 and S4
data <- data %>% select(BT16_doxorubicin, PEDS0005_T1_doxorubicin, `CHLA-06-ATRT_doxorubicin`, PEDS9001_T1_doxorubicin, HepG2_doxorubicin_clone1,
                        HepG2_doxorubicin_clone2, MCF10A_doxorubicin_clone1, MCF10A_doxorubicin_clone2)


# count number of mutations for each sample
total <- colSums(data) |> as.data.frame()

# format the table 
total$Sample <- rownames(total)
colnames(total)[1] <- "Total number of single-base substitutions"
total_final <- total[,c(2,1)]
rownames(total_final) <- NULL

# save table
write_tsv(x = total_final, "../tables/TableS5.tsv")

rm(data, total, total_final); gc() # clean
