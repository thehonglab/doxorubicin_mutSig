# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for figures presented in Supplemental Figure 3 of manuscript

# load libraries
library(sigProcessor)

# set working directory
setwd("~/doxorubicin_mutSig/data/")

## Figure S3A
  # created in python

## Figure S3B
# load cisplatin SBS mutational matrices
cisplatin <- cbind(process_mut_types_sbs("./SBS96/cisplatin/BT16_cisplatin.SBS96.all", "BT16_cisplatin"),
                   process_mut_types_sbs("./SBS96/cisplatin/PEDS0005_T1_cisplatin.SBS96.all", "PEDS0005_T1_cisplatin"),
                   process_mut_types_sbs("./SBS96/cisplatin/HepG2_clone1.SBS96.all", "HepG2_cisplatin_1"),
                   process_mut_types_sbs("./SBS96/cisplatin/HepG2_clone2.SBS96.all", "HepG2_cisplatin_2"),
                   process_mut_types_sbs("./SBS96/cisplatin/HepG2_clone3.SBS96.all", "HepG2_cisplatin_3"),
                   process_mut_types_sbs("./SBS96/cisplatin/HepG2_clone4.SBS96.all", "HepG2_cisplatin_4"),
                   process_mut_types_sbs("./SBS96/cisplatin/MCF10A_05uM_cis_4wks.SBS96.all", "MCF10A_cisplatin_1"),
                   process_mut_types_sbs("./SBS96/cisplatin/MCF10A_1uM_cis_4wks.SBS96.all", "MCF10A_cisplatin_2"),
                   process_mut_types_sbs("./SBS96/cisplatin/MCF10A_05uM_cis_8wks_clone1.SBS96.all", "MCF10A_cisplatin_3"),
                   process_mut_types_sbs("./SBS96/cisplatin/MCF10A_05uM_cis_8wks_clone2.SBS96.all", "MCF10A_cisplatin_4"),
                   process_mut_types_sbs("./SBS96/cisplatin/MCF10A_1uM_cis_8wks_clone1.SBS96.all", "MCF10A_cisplatin_5"),
                   process_mut_types_sbs("./SBS96/cisplatin/MCF10A_1uM_cis_8wks_clone2.SBS96.all", "MCF10A_cisplatin_6"))

# convert to percentages
cisplatin_percent <- freq2percent(cisplatin)

# generate cosine similarity matrix
cisplatin_cosim <- get_cosine_sim_matrix(cisplatin_percent)

# plot matrix 
p <- plot_cosine_sim(cisplatin_cosim, "SBS - Cisplatin")
print(p)

# save plot
ggsave(plot = p, filename = "../figures/FigS3/Figure_S3B.png", height = 7, width = 9)

# clean environment
rm(cisplatin, cisplatin_cosim, cisplatin_percent, p); gc()

## Figure S3C
  # created in python

## Figure S3D
# load cisplatin DBS mutational matrices
cisplatin <- cbind(process_mut_types_dbs("./DBS78/cisplatin/BT16_cisplatin.DBS78.all", "BT16_cisplatin"),
                   process_mut_types_dbs("./DBS78/cisplatin/PEDS0005_T1_cisplatin.DBS78.all", "PEDS0005_T1_cisplatin"),
                   process_mut_types_dbs("./DBS78/cisplatin/HepG2_clone1.DBS78.all", "HepG2_cisplatin_1"),
                   process_mut_types_dbs("./DBS78/cisplatin/HepG2_clone2.DBS78.all", "HepG2_cisplatin_2"),
                   process_mut_types_dbs("./DBS78/cisplatin/HepG2_clone3.DBS78.all", "HepG2_cisplatin_3"),
                   process_mut_types_dbs("./DBS78/cisplatin/HepG2_clone4.DBS78.all", "HepG2_cisplatin_4"),
                   process_mut_types_dbs("./DBS78/cisplatin/MCF10A_05uM_cis_4wks.DBS78.all", "MCF10A_cisplatin_1"),
                   process_mut_types_dbs("./DBS78/cisplatin/MCF10A_1uM_cis_4wks.DBS78.all", "MCF10A_cisplatin_2"),
                   process_mut_types_dbs("./DBS78/cisplatin/MCF10A_05uM_cis_8wks_clone1.DBS78.all", "MCF10A_cisplatin_3"),
                   process_mut_types_dbs("./DBS78/cisplatin/MCF10A_05uM_cis_8wks_clone2.DBS78.all", "MCF10A_cisplatin_4"),
                   process_mut_types_dbs("./DBS78/cisplatin/MCF10A_1uM_cis_8wks_clone1.DBS78.all", "MCF10A_cisplatin_5"),
                   process_mut_types_dbs("./DBS78/cisplatin/MCF10A_1uM_cis_8wks_clone2.DBS78.all", "MCF10A_cisplatin_6"))

# convert to percentages
cisplatin_percent <- freq2percent(cisplatin)

# generate cosine similarity matrix
cisplatin_cosim <- get_cosine_sim_matrix(cisplatin_percent)

# plot matrix 
p <- plot_cosine_sim(cisplatin_cosim, "DBS - Cisplatin")
print(p)

# save plot
ggsave(plot = p, filename = "../figures/FigS3/Figure_S3D.png", height = 7, width = 9)

# clean environment
rm(cisplatin, cisplatin_cosim, cisplatin_percent, p); gc()
