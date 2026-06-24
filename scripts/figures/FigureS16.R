# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for figures presented in Figure S16 of manuscript

# load libraries
library(sigProcessor)
library(tidyverse)

# set working directory
setwd("~/doxorubicin_mutSig/data/")

# load in doxorubicin signatures
doxorubicin <- cbind(process_mut_types_sbs("./SBS96/doxorubicin/SBS-Doxorubicin-1.SBS96.all", "SBS-Doxorubicin-1"),
                     process_mut_types_sbs("./SBS96/doxorubicin/SBS-Doxorubicin-2.SBS96.all", "SBS-Doxorubicin-2")) |> as.data.frame()

# load in Mitchell et al patient data
# acquired from: https://github.com/emily-mitchell/chemotherapy/blob/858f01318acca0c97c9e106296577613f071f648/5_Mutational_signature_analysis/mutational_signatures_analysis/trinuc_mut_mat.txt
# pre-subset for the immune populations of the following 3 patients: PD47540, PD47695, PD47696
# for this figure, subset for the first patient
patient <- read.table("./SBS96/patient_data/Mitchell_2025.SBS96.all", header = T)
patient <- cbind(patient %>% dplyr::select(MutationType, all_of(colnames(patient)[grepl("PD47695", colnames(patient))])),
                 patient %>% dplyr::select(all_of(colnames(patient)[grepl("PD47696", colnames(patient))])))

# double check rows are in the same order
(rownames(doxorubicin) == convert_mut_types(patient$MutationType)) |> table()

# combine dfs
data <- cbind(patient[,2:7], doxorubicin)

# convert to %
data_percent <- freq2percent(data)

# create cosine similarities matrix
data_cosim <- get_cosine_sim_matrix(data_percent)

# plot
p <- plot_cosine_sim(data_cosim, "SBS doxorubicin in PD47695 & PD47696 immune populations\n                                 (Mitchell, et al 2025)")
print(p)

# save plot
ggsave(plot = p, filename = "../figures/FigS15/Figure_S15.png", height = 6, width = 8)

# clean environment
rm(data, data_cosim, data_percent, doxorubicin, patient, p); gc()
