# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for figures presented in Figure 3 of manuscript

# load libraries
library(sigProcessor)
library(tidyverse)

# set working directory
setwd("~/doxorubicin_mutSig/data/")

## Figure 3A
  # created in biorender

## Figure 3B
  # created in python

## Figure 3C
# load in PEDS0005_T2A/B patient profiles
patient <- process_mutation_matrix(c("./SBS96/patient_data/PEDS0005_T2A.SBS96.all", 
                                     "./SBS96/patient_data/PEDS0005_T2B.SBS96.all"))

# load in samples comprising SBS-Doxorubicin-1, minus PEDS0005_T1_doxorubicin
doxo1 <- process_mutation_matrix(c("./SBS96/doxorubicin/BT16_doxorubicin.SBS96.all",
                                   "./SBS96/doxorubicin/HepG2_doxorubicin_clone2.SBS96.all",
                                   "./SBS96/doxorubicin/MCF10A_doxorubicin_clone1.SBS96.all",
                                   "./SBS96/doxorubicin/MCF10A_doxorubicin_clone2.SBS96.all"))

# calculate the average of these 4 samples
doxo1_av <- rowMeans(doxo1) |> as.data.frame()

# rename column
colnames(doxo1_av) <- "SBS-Doxorubicin-1\n(-PEDS0005_T1_doxorubicin)"

# load in second doxorubicin signature
doxo2 <- process_mutation_matrix("./SBS96/doxorubicin/SBS-Doxorubicin-2.SBS96.all")

# convert to %
doxo_percent <- freq2percent(cbind(doxo1_av, doxo2))

# convert patient data to %
patient_percent <- freq2percent(patient)

# combine dfs
data <- cbind(doxo_percent, patient_percent)

# create cosine similarity matrix
data_cosim <- get_cosine_sim_matrix(data)

# plot
p <- plot_cosine_sim(data_cosim, "Doxorubicin SBS in PEDS0005 metastatic samples")
print(p)

# save plot
ggsave(plot = p, filename = "../figures/Fig3/Figure_3C.png", height = 5, width = 7.5)

# clean environment
rm(data, data_cosim, doxo_percent, doxo1, doxo1_av, doxo2, patient, patient_percent, p); gc()

## Figure 3D
# load in PEDS0005_T2A/B patient profiles
patient <- process_mutation_matrix(c("./SBS96/patient_data/PEDS0005_T2A.SBS96.all", 
                                     "./SBS96/patient_data/PEDS0005_T2B.SBS96.all"))

# load in cisplatin spectra: BT16, MCF10A_c1, and HepG2_c1
cisplatin <- process_mutation_matrix(c("./SBS96/cisplatin/BT16_cisplatin.SBS96.all",
                                       "./SBS96/cisplatin/HepG2_clone1.SBS96.all",
                                       "./SBS96/cisplatin/MCF10A_05uM_cis_4wks.SBS96.all"))

# combine dfs
data <- cbind(cisplatin, patient)

# convert to %
data_percent <- freq2percent(data)

# create cosine similarity matrix
data_cosim <- get_cosine_sim_matrix(data_percent)

# plot
p <- plot_cosine_sim(data_cosim, "Cisplatin SBS in PEDS0005 metastatic samples")
print(p)

# save plot
ggsave(plot = p, filename = "../figures/Fig3/Figure_3D.png", height = 5, width = 7.5)

# clean environment
rm(cisplatin, data, data_cosim, data_percent, patient, p); gc()

## Figure 3E
  # created in biorender

## Figure 3F
  # created in python

## Figure 3G
# load in doxorubicin signatures
doxorubicin <- cbind(process_mut_types_sbs("./SBS96/doxorubicin/SBS-Doxorubicin-1.SBS96.all", "SBS-Doxorubicin-1"),
                     process_mut_types_sbs("./SBS96/doxorubicin/SBS-Doxorubicin-2.SBS96.all", "SBS-Doxorubicin-2")) |> as.data.frame()

# load in Mitchell et al patient data
# acquired from: https://github.com/emily-mitchell/chemotherapy/blob/858f01318acca0c97c9e106296577613f071f648/5_Mutational_signature_analysis/mutational_signatures_analysis/trinuc_mut_mat.txt
# pre-subset for the immune populations of the following 3 patients: PD47540, PD47695, PD47696
# for this figure, subset for the first patient
patient <- read.table("./SBS96/patient_data/Mitchell_2025.SBS96.all", header = T)
patient <- patient %>% dplyr::select(MutationType, all_of(colnames(patient)[grepl("PD47540", colnames(patient))]))

# double check rows are in the same order
(rownames(doxorubicin) == convert_mut_types(patient$MutationType)) |> table()

# combine dfs
data <- cbind(patient[,2:5], doxorubicin)

# convert to %
data_percent <- freq2percent(data)

# create cosine similarities matrix
data_cosim <- get_cosine_sim_matrix(data_percent)

# plot
p <- plot_cosine_sim(data_cosim, "         SBS doxorubicin in PD47540\nimmune populations (Mitchell, et al 2025)")
print(p)

# save plot
ggsave(plot = p, filename = "../figures/Fig3/Figure_3G.png", height = 5, width = 7.5)

# clean environment
rm(data, data_cosim, data_percent, doxorubicin, patient, p); gc()



