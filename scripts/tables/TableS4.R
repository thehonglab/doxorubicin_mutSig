# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for Table S4 of manuscript

# load libraries
library(sigProcessor)
library(tidyverse)

# set working directory
setwd("~/doxorubicin_mutSig/data/DBS78/")

# load in matrices
data <- process_mutation_matrix(c("./cisplatin/BT16_cisplatin.DBS78.all",
                                  "./doxorubicin/BT16_doxorubicin.DBS78.all",
                                  "./cisplatin/PEDS0005_T1_cisplatin.DBS78.all",
                                  "./doxorubicin/PEDS0005_T1_doxorubicin.DBS78.all",
                                  "./doxorubicin/CHLA-06-ATRT_doxorubicin.DBS78.all",
                                  "./doxorubicin/PEDS9001_T1_doxorubicin.DBS78.all",
                                  "./cisplatin/HepG2_clone1.DBS78.all",
                                  "./cisplatin/HepG2_clone2.DBS78.all",
                                  "./cisplatin/HepG2_clone3.DBS78.all",
                                  "./cisplatin/HepG2_clone4.DBS78.all",
                                  "./doxorubicin/HepG2_doxorubicin_clone1.DBS78.all",
                                  "./doxorubicin/HepG2_doxorubicin_clone2.DBS78.all",
                                  "./cisplatin/MCF10A_05uM_cis_4wks.DBS78.all",
                                  "./cisplatin/MCF10A_1uM_cis_4wks.DBS78.all",
                                  "./cisplatin/MCF10A_05uM_cis_8wks_clone1.DBS78.all",
                                  "./cisplatin/MCF10A_05uM_cis_8wks_clone2.DBS78.all",
                                  "./cisplatin/MCF10A_1uM_cis_8wks_clone1.DBS78.all",
                                  "./cisplatin/MCF10A_1uM_cis_8wks_clone2.DBS78.all",
                                  "./doxorubicin/MCF10A_doxorubicin_clone1.DBS78.all",
                                  "./doxorubicin/MCF10A_doxorubicin_clone2.DBS78.all"))

# count number of mutations for each sample
total <- colSums(data) |> as.data.frame()

# format the table 
total$Sample <- rownames(total)
colnames(total)[1] <- "Total number of single-base substitutions"
total_final <- total[,c(2,1)]
rownames(total_final) <- NULL

# save table
write_tsv(x = total_final, "../tables/TableS4.tsv")

rm(data, total, total_final); gc() # clean
