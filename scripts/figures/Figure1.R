# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for figures presented in Figure 1 of manuscript

# load libraries
library(sigProcessor)
library(tidyverse)

# set working directory
setwd("~/doxorubicin_mutSig/data/")

## Figure 1A
  # created in biorender

## Figure 1B
  # created in python

## Figure 1C
# load cisplatin SBS mutational matrices
cisplatin <- process_mutation_matrix(c("./SBS96/cisplatin/BT16_cisplatin.SBS96.all",
                                       "./SBS96/cisplatin/PEDS0005_T1_cisplatin.SBS96.all",
                                       "./SBS96/cisplatin/HepG2_clone1.SBS96.all",
                                       "./SBS96/cisplatin/MCF10A_05uM_cis_4wks.SBS96.all"))

# convert to percentages
cisplatin_percent <- freq2percent(cisplatin)

# generate cosine similarity matrix
cisplatin_cosim <- get_cosine_sim_matrix(cisplatin_percent)

# plot matrix 
p <- plot_cosine_sim(cisplatin_cosim, "SBS - Cisplatin")
print(p)

# save plot
ggsave(plot = p, filename = "../figures/Fig1/Figure_1C.png", height = 4, width = 6)

# clean environment
rm(cisplatin, cisplatin_cosim, cisplatin_percent, p); gc()

## Figure 1D
  # created in python

## Figure 1E
# load cisplatin DBS mutational matrices
cisplatin <- process_mutation_matrix(c("./DBS78/cisplatin/BT16_cisplatin.DBS78.all",
                                       "./DBS78/cisplatin/PEDS0005_T1_cisplatin.DBS78.all",
                                       "./DBS78/cisplatin/HepG2_clone1.DBS78.all",
                                       "./DBS78/cisplatin/MCF10A_05uM_cis_4wks.DBS78.all"))

# convert to percentages
cisplatin_percent <- freq2percent(cisplatin)

# generate cosine similarity matrix
cisplatin_cosim <- get_cosine_sim_matrix(cisplatin_percent)

# plot matrix 
p <- plot_cosine_sim(cisplatin_cosim, "DBS - Cisplatin")
print(p)

# save plot
ggsave(plot = p, filename = "../figures/Fig1/Figure_1E.png", height = 4, width = 6)

# clean environment
rm(cisplatin, cisplatin_cosim, cisplatin_percent, p); gc()
