# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for figures presented in Supplemental Figure 7 of manuscript

# load libraries
library(sigProcessor)
library(tidyverse)
library(data.table)

# set working directory
setwd("~/doxorubicin_mutSig/data/")

## Figure S7A-E
# read in SBS VCFs
sig1 <- prepare_sbs96_vcf(c("/path/to/BT16_doxorubicin_SBS_VCF", 
                            "/path/to/PEDS0005_T1_doxorubicin_SBS_VCF",
                            "/path/to/HepG2_doxorubicin_clone2_SBS_VCF",
                            "/path/to/MCF10A_doxorubicin_clone1_SBS_VCF",
                            "/path/to/MCF10A_doxorubicin_clone2_SBS_VCF")) # data not included due to controlled access

# prepare data for plotting
sig1_plot_df <- prepare_rainfall_plot_df(sig1)

# check whether any of the log10(intermutational distances) are -Inf (multiple mutations at the same genomic loci)
sig1_plot_df %>% filter(log10_dist == "-Inf") # none

# acquire a list of samples
sig1_samples <- sig1_plot_df$sample |> unique()

# rearrange samples
sig1_samples <- sig1_samples[c(1,5,2,3,4)]

# define subpanel letters
panel <- c("A", "B", "C", "D", "E")

# iterate through samples and plot
for (i in 1:length(sig1_samples)) {
  
  # subset df
  df <- sig1_plot_df %>% filter(sample == sig1_samples[i])
  
  # plot
  p <- plot_rainfall_plot(df, plot.title = paste((df$sample |> unique()),"(SBS-Doxorubicin-1)"))
  
  # save
  filename <- paste0("../figures/FigS6/Figure_S6", panel[i], ".png")
  ggsave(filename = filename, plot = p, height = 3, width = 12)
  
}

# clean environment
rm(list = ls()); gc()