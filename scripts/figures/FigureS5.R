# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for figures presented in Supplemental Figure 5 of manuscript

# load libraries
library(sigProcessor)
library(tidyverse)

# set working directory
setwd("~/doxorubicin_mutSig/data/")

## Figure S5A
  # created in python

## Figure S5B
  # created in python

## Figure S5C
# load doxorubicin and SMARCB1-cisplatin SBS mutational matrices
data <- process_mutation_matrix(c("./SBS96/doxorubicin/BT16_doxorubicin.SBS96.all",
                                  "./SBS96/doxorubicin/CHLA-06-ATRT_doxorubicin.SBS96.all",
                                  "./SBS96/doxorubicin/PEDS0005_T1_doxorubicin.SBS96.all",
                                  "./SBS96/doxorubicin/PEDS9001_T1_doxorubicin.SBS96.all",
                                  "./SBS96/doxorubicin/HepG2_doxorubicin_clone1.SBS96.all",
                                  "./SBS96/doxorubicin/HepG2_doxorubicin_clone2.SBS96.all",
                                  "./SBS96/doxorubicin/MCF10A_doxorubicin_clone1.SBS96.all",
                                  "./SBS96/doxorubicin/MCF10A_doxorubicin_clone2.SBS96.all",
                                  "./SBS96/cisplatin/BT16_cisplatin.SBS96.all",
                                  "./SBS96/cisplatin/PEDS0005_T1_cisplatin.SBS96.all"))

# convert to percentages
data_percent <- freq2percent(data)

# generate cosine similarity matrix
data_cosim <- get_cosine_sim_matrix(data_percent)

# plot matrix 
p <- plot_cosine_sim(data_cosim, "SBS - Doxorubicin compared to cisplatin")
print(p)

# save plot
ggsave(plot = p, filename = "../figures/FigS4/Figure_S4C.png", height = 7, width = 9)

# clean environment
rm(data, data_cosim, data_percent, p); gc()

## Figure S5D
  # created in python
