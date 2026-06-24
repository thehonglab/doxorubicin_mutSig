# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for figures presented in Supplemental Figure 6 of manuscript

# load libraries
library(sigProcessor)
library(tidyverse)
library(indelsig.tools.lib)

# load in additional functions
source("~/doxorubicin_mutSig/scripts/additional_functions/ID30_functions.R")

# set working directory
setwd("~/doxorubicin_mutSig/data/")

## Figure S6A-H
# read in INDEL VCFs
indels <- rbind(prepare_indel_vcf("/path/to/BT16_doxorubicin_INDEL_VCF"),
                prepare_indel_vcf("/path/to/CHLA-06-ATRT_doxorubicin_INDEL_VCF"),
                prepare_indel_vcf("/path/to/PEDS0005_T1_doxorubicin_INDEL_VCF"),
                prepare_indel_vcf("/path/to/PEDS9001_T1_doxorubicin_INDEL_VCF"),
                prepare_indel_vcf("/path/to/HepG2_doxorubicin_clone1_INDEL_VCF"),
                prepare_indel_vcf("/path/to/HepG2_doxorubicin_clone2_INDEL_VCF"),
                prepare_indel_vcf("/path/to/MCF10A_doxorubicin_clone1_INDEL_VCF"),
                prepare_indel_vcf("/path/to/MCF10A_doxorubicin_clone2_INDEL_VCF")) # data not included due to controlled access

# prepare indel sequence context information for classification
indels_prepared <- indelsig.tools.lib::prepare_indels(indels, "doxorubicin", "hg38")

# from ID30_functions.R script (loaded above with packages), derived from functions shared with me by Gene Ching Chiek Koh (Nik Zainal lab)
# call legacy indels (ID30 classification as found in the Nature Genetics 2023 manuscript)
indels_called <- call_legacy_indels(indels_prepared)

# now generate the full matrix
indels_cat <- get_matrix_ID30(indels_called)

# save indels matrix
indels_final <- indels_cat
indels_final$MutationType <- rownames(indels_final)
rownames(indels_final) <- NULL
indels_final <- indels_final[,c(9, 1:8)]

# save data 
write_tsv(indels_final, "./ID30/doxorubicin.ID30.all")

# plotting
# define subpanel letters 
panel <- c("A", "B", "C", "D", "E", "F", "G", "H")

# re-order samples based on subpanel order
indels_ordered <- indels_cat %>% dplyr::select(BT16_doxorubicin, PEDS0005_T1_doxorubicin, HepG2_doxorubicin_clone2, MCF10A_doxorubicin_clone1,
                                               MCF10A_doxorubicin_clone2, `CHLA-06-ATRT_doxorubicin`, PEDS9001_T1_doxorubicin, 
                                               HepG2_doxorubicin_clone1)

# lastly, set up a counter
counter = 0

# in a loop, plot each sample and save
for (s in colnames(indels_ordered)) {

  # increase a counter
  counter <- counter + 1
  
  # subset data frame for the sample
  data <- indels_ordered %>% dplyr::select(all_of(s))

  # plot
  p <- gen_plot_catalouge_legacy_single(data, 6, s)

  # save
  filename <- paste0("../figures/FigS5/Figure_S5", panel[counter], ".png")
  ggsave(plot = p, filename = filename, width = 12, height = 6)

}

# clean environment
rm(list = ls()); gc()
