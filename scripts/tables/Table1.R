# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for Table 1 of manuscript

# load libraries
library(sigProcessor)
library(tidyverse)
library(mSigAct)
library(cosmicsig)

# set working directory
setwd("~/doxorubicin_mutSig/data/")

# load in PEDS0005_T2A/B patient profiles
patient <- cbind(process_mut_types_sbs("./SBS96/patient_data/PEDS0005_T2A.SBS96.all", "PEDS0005_T2A"),
                 process_mut_types_sbs("./SBS96/patient_data/PEDS0005_T2B.SBS96.all", "PEDS0005_T2B")) |> as.data.frame()

# load in samples comprising SBS-Doxorubicin-1, minus PEDS0005_T1_doxorubicin
doxo1 <- cbind(process_mut_types_sbs("./SBS96/doxorubicin/BT16_doxorubicin.SBS96.all", "BT16_doxorubicin"),
               process_mut_types_sbs("./SBS96/doxorubicin/HepG2_doxorubicin_clone2.SBS96.all", "HepG2_doxorubicin_clone2"),
               process_mut_types_sbs("./SBS96/doxorubicin/MCF10A_doxorubicin_clone1.SBS96.all", "MCF10A_doxorubicin_clone1"),
               process_mut_types_sbs("./SBS96/doxorubicin/MCF10A_doxorubicin_clone2.SBS96.all", "MCF10A_doxorubicin_clone2"))

# calculate the average of these 4 samples
doxo1_av <- rowMeans(doxo1) |> as.data.frame()

# rename column
colnames(doxo1_av) <- "SBS-Doxorubicin-1"

# load in second doxorubicin signature
doxo2 <- process_mut_types_sbs("./SBS96/doxorubicin/SBS-Doxorubicin-2.SBS96.all", "SBS-Doxorubicin-2")

# convert to %
doxo_percent <- freq2percent(cbind(doxo1_av, doxo2))

# import cisplain spectra
cis <- cbind(process_mut_types_sbs("./SBS96/cisplatin/MCF10A_05uM_cis_4wks.SBS96.all", "MCF10A_cisplatin_1"),
             process_mut_types_sbs("./SBS96/cisplatin/MCF10A_1uM_cis_4wks.SBS96.all", "MCF10A_cisplatin_2"),
             process_mut_types_sbs("./SBS96/cisplatin/MCF10A_05uM_cis_8wks_clone1.SBS96.all", "MCF10A_cisplatin_3"),
             process_mut_types_sbs("./SBS96/cisplatin/MCF10A_05uM_cis_8wks_clone2.SBS96.all", "MCF10A_cisplatin_4"),
             process_mut_types_sbs("./SBS96/cisplatin/MCF10A_1uM_cis_8wks_clone1.SBS96.all", "MCF10A_cisplatin_5"),
             process_mut_types_sbs("./SBS96/cisplatin/MCF10A_1uM_cis_8wks_clone2.SBS96.all", "MCF10A_cisplatin_6"),
             process_mut_types_sbs("./SBS96/cisplatin/HepG2_clone1.SBS96.all", "HepG2_cisplatin_1"),
             process_mut_types_sbs("./SBS96/cisplatin/HepG2_clone2.SBS96.all", "HepG2_cisplatin_2"),
             process_mut_types_sbs("./SBS96/cisplatin/HepG2_clone3.SBS96.all", "HepG2_cisplatin_3"),
             process_mut_types_sbs("./SBS96/cisplatin/HepG2_clone4.SBS96.all", "HepG2_cisplatin_4"))

# take average of all spectra
cis_av <- rowMeans(cis) |> as.data.frame() 

# normalise to percentages
cis_percent <- freq2percent(cis_av) 

# rename column
colnames(cis_percent) <- "SBS-Cisplatin" 

# get COSMIC signatures
sigs_all <- COSMIC_v3.3$signature$GRCh38$SBS96

# double check column names are correct
rownames(sigs_all) == rownames(doxo_percent)
rownames(sigs_all) == rownames(cis_percent)
rownames(sigs_all) == rownames(patient)

# combine into one df
sigs <- cbind(sigs_all[,c("SBS1", "SBS5", "SBS18")], doxo_percent, cis_percent)

# clean environment
rm(cis, cis_av, cis_percent, doxo_percent, doxo1, doxo1_av, doxo2, sigs_all); gc()

# test presence of cisplatin
res <- SignaturePresenceTest(patient, 
                             sigs = as.matrix(sigs), 
                             target.sig.index = "SBS-Cisplatin", 
                             seed = 1234,
                             mc.cores = 2)

# save in a table
res.table <- data.frame(Sample = c("PEDS0005_T2A", "PEDS0005_T2B"),
                        Signature_tested = rep("SBS-Cisplatin", 2),
                        chisq = c(signif(res$PEDS0005_T2A$chisq.p, digits = 3),
                                  signif(res$PEDS0005_T2B$chisq.p, digits = 3)))

# test presence of first doxorubicin signature
res <- SignaturePresenceTest(patient, 
                             sigs = as.matrix(sigs[,c(1:4,6)]), 
                             target.sig.index = "SBS-Doxorubicin-1", 
                             seed = 1234,
                             mc.cores = 2)

# save in a table
res.table <- rbind(res.table, 
                   c("PEDS0005_T2A", "SBS-Doxorubicin-1 (-PEDS0005_T1)", signif(res$PEDS0005_T2A$chisq.p, digits = 3)),
                   c("PEDS0005_T2B", "SBS-Doxorubicin-1 (-PEDS0005_T1)", signif(res$PEDS0005_T2B$chisq.p, digits = 3)))

# test presence of second doxorubicin signature
res <- SignaturePresenceTest(patient, 
                             sigs = as.matrix(sigs[,c(1:3,5,6)]), 
                             target.sig.index = "SBS-Doxorubicin-2", 
                             seed = 1234,
                             mc.cores = 2)

# save in a table
res.table <- rbind(res.table, 
                   c("PEDS0005_T2A", "SBS-Doxorubicin-2", signif(res$PEDS0005_T2A$chisq.p, digits = 3)),
                   c("PEDS0005_T2B", "SBS-Doxorubicin-2", signif(res$PEDS0005_T2B$chisq.p, digits = 3)))

# rework table for manuscript
final_table <- data.frame(`Sample ID` = c("PEDS0005_T2A", "PEDS0005_T2B"),
                          `Cancer type` = rep("RMC", 2),
                          `Cisplatin mSigAct p-value` = (res.table %>% filter(Signature_tested == "SBS-Cisplatin") %>% select(chisq))$chisq,
                          `SBS-Doxorubicin-1 (-PEDS0005_T1) mSigAct p-value` = (res.table %>% filter(Signature_tested == "SBS-Doxorubicin-1 (-PEDS0005_T1)") %>% select(chisq))$chisq,
                          `SBS-Doxorubicin-2 mSigAct p-value` = (res.table %>% filter(Signature_tested == "SBS-Doxorubicin-2") %>% select(chisq))$chisq,
                          check.names = F)

# save table
write_tsv(final_table, "./tables/Table1.tsv")


# clean environment
rm(final_table, patient, res, res.table, sigs); gc()
