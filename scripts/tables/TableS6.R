# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for Table S6 of manuscript

# load libraries
library(sigProcessor)
library(tidyverse)

# set working directory
setwd("~/doxorubicin_mutSig/data/")

# import the activities table
data <- read.table("./SigProfilerAssignment/output/average_doxorubicin_activities.txt", header = T)

# remove columns that are all zero
data <- data %>% select(where(~ any(.x != 0)))

# change samples to signature
colnames(data)[1] <- "Doxorubicin signature"

# save
write_tsv(x = data, file = "./tables/TableS6.tsv")

rm(data); gc()
