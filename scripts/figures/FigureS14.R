# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for figures presented in Supplemental Figure S14 of manuscript

# load libraries
library(sigProcessor)

# set working directory
setwd("~/doxorubicin_mutSig/data/")

# load in SBS96 matrices
data <- process_mutation_matrix(c("./SBS96/cisplatin/BT16_cisplatin.SBS96.all",
                                  "./SBS96/cisplatin/PEDS0005_T1_cisplatin.SBS96.all",
                                  "./SBS96/cisplatin/HepG2_clone1.SBS96.all",
                                  "./SBS96/cisplatin/HepG2_clone2.SBS96.all",
                                  "./SBS96/cisplatin/HepG2_clone3.SBS96.all",
                                  "./SBS96/cisplatin/HepG2_clone4.SBS96.all",
                                  "./SBS96/cisplatin/MCF10A_05uM_cis_4wks.SBS96.all",
                                  "./SBS96/cisplatin/MCF10A_1uM_cis_4wks.SBS96.all",
                                  "./SBS96/cisplatin/MCF10A_05uM_cis_8wks_clone1.SBS96.all",
                                  "./SBS96/cisplatin/MCF10A_05uM_cis_8wks_clone2.SBS96.all",
                                  "./SBS96/cisplatin/MCF10A_1uM_cis_8wks_clone1.SBS96.all",
                                  "./SBS96/cisplatin/MCF10A_1uM_cis_8wks_clone2.SBS96.all",
                                  "./SBS96/doxorubicin/BT16_doxorubicin.SBS96.all",
                                  "./SBS96/doxorubicin/CHLA-06-ATRT_doxorubicin.SBS96.all",
                                  "./SBS96/doxorubicin/PEDS0005_T1_doxorubicin.SBS96.all",
                                  "./SBS96/doxorubicin/PEDS9001_T1_doxorubicin.SBS96.all",
                                  "./SBS96/doxorubicin/HepG2_doxorubicin_clone1.SBS96.all",
                                  "./SBS96/doxorubicin/HepG2_doxorubicin_clone2.SBS96.all",
                                  "./SBS96/doxorubicin/MCF10A_doxorubicin_clone1.SBS96.all",
                                  "./SBS96/doxorubicin/MCF10A_doxorubicin_clone2.SBS96.all",
                                  "./SBS96/patient_data/PEDS0005_T2A.SBS96.all",
                                  "./SBS96/patient_data/PEDS0005_T2B.SBS96.all"))

# convert to %
data_percent <- freq2percent(data)

# calculate cosine similarities
data_cosim <- get_cosine_sim_matrix(data_percent)

# plot
p <- plot_cosine_sim(data_cosim, "SBS doxorubicin and cisplatin as compared to PEDS0005 metastatic samples")
print(p)

# save
ggsave(plot = p, filename = "../figures/FigS13/Figure_S13.png", width = 12, height = 10)

# clean environment
rm(data, data_cosim, data_percent, p); gc()
 