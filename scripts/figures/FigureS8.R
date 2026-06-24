# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for figures presented in Supplemental Figure 8 of manuscript

# load libraries
library(sigProcessor)
library(tidyverse)
library(data.table)

# set working directory
setwd("~/doxorubicin_mutSig/data/")

## Figure S8A-C
# read in SBS VCFs
sig2 <- prepare_sbs96_vcf(c("/path/to/PEDS9001_T1_doxorubicin_SBS_VCF",
                            "/path/to/CHLA-06-ATRT_doxorubicin_SBS_VCF",
                            "/path/to/HepG2_doxorubicin_clone1_SBS_VCF")) # data not included due to controlled access

# prepare data for plotting
sig2_plot_df <- prepare_rainfall_plot_df(sig2)

# check whether any of the log10(intermutational distances) are -Inf (multiple mutations at the same genomic loci)
sig2_plot_df %>% filter(log10_dist == "-Inf") # 12

# inspect the mutations that have log10(intermutational distance) == -Inf
positions <- sig2_plot_df %>% filter(Position %in% (sig2_plot_df %>% filter(log10_dist == "-Inf") %>% select(Position))[,1])
table(positions$sample) # dinucleotide substitutions that have been split into two single nucleotide substitutions
  
# filter out mutations that have log10(intermutational distance == -Inf)
# they would just be on top of one another on the plot anyway
sig2_plot_df <- sig2_plot_df %>% filter(!log10_dist == "-Inf")
sig2_plot_df %>% filter(log10_dist == "-Inf") # sanity check, now none

# acquire a list of samples
sig2_samples <- sig2_plot_df$sample |> unique()

# rearrange samples
sig2_samples <- sig2_samples[c(3,1,2)]

# define subpanel letters
panel <- c("A", "B", "C")

# iterate through samples and plot
for (i in 1:length(sig2_samples)) {
  
  # subset df
  df <- sig2_plot_df %>% filter(sample == sig2_samples[i])
  
  # plot
  p <- plot_rainfall_plot(df, plot.title = paste((df$sample |> unique()),"(SBS-Doxorubicin-2)"))
  
  # save
  filename <- paste0("../figures/FigS7/Figure_S7", panel[i], ".png")
  ggsave(filename = filename, plot = p, height = 3, width = 12)
  
}

# clean environment
rm(list = ls()); gc()