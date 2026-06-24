# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for Table 2 of manuscript

# load libraries
library(sigProcessor)
library(tidyverse)
library(mSigAct)
library(cosmicsig)

# set working directory
setwd("~/doxorubicin_mutSig/data/")

# load in Mitchell et al patient data
# acquired from: https://github.com/emily-mitchell/chemotherapy/blob/858f01318acca0c97c9e106296577613f071f648/5_Mutational_signature_analysis/mutational_signatures_analysis/trinuc_mut_mat.txt
# pre-subset for the immune populations of the following 3 patients: PD47540, PD47695, PD47696
# for this figure, subset for the first patient
patient <- read.table("./SBS96/patient_data/Mitchell_2025.SBS96.all", header = T)
rownames(patient) <- convert_mut_types(patient$MutationType)
patient$MutationType <- NULL

# load in doxorubicin signatures
doxorubicin <- cbind(process_mut_types_sbs("./SBS96/doxorubicin/SBS-Doxorubicin-1.SBS96.all", "SBS-Doxorubicin-1"),
                     process_mut_types_sbs("./SBS96/doxorubicin/SBS-Doxorubicin-2.SBS96.all", "SBS-Doxorubicin-2")) |> as.data.frame()

# load in cyclophosphamide from COSMIC experimental signatures database
cyclo <- process_mut_types_sbs("./SBS96/exposures/experimental_cyclophosphamide.SBS96.all", "SBS-Cyclophosphamide") |> as.data.frame()

# load in canonical signatures found to be active in patients 
sigs_all <- cbind(process_mut_types_sbs("./SBS96/exposures/Mitchell_2025_SBS1_SBS5.SBS96.all", "SBS1+SBS5"),
                  process_mut_types_sbs("./SBS96/exposures/Mitchell_2025_SBSBlood.SBS96.all", "SBSBlood")) |> as.data.frame()

# get COSMIC signatures
cosmic <- COSMIC_v3.3$signature$GRCh38$SBS96 |> as.data.frame()

# ensure all rows are in the same order
rownames(doxorubicin) == rownames(patient)
rownames(doxorubicin) == rownames(cyclo)
rownames(doxorubicin) == rownames(sigs_all)
rownames(doxorubicin) == rownames(cosmic)

# load in other exposure for these samples from manuscript
exposures <- read.table("./SBS96/exposures/Mitchell_2025_patient_exposures.txt", header = T, sep = "\t", check.names = F)

# check which signatures need to be converted to %
colSums(doxorubicin) 
colSums(cyclo)
colSums(sigs_all)
colSums(cosmic)
  # only doxorubicin

doxorubicin <- freq2percent(doxorubicin)
colSums(doxorubicin)

# subset the COSMIC list for COSMIC signatures found to be active in the patient samples
# not including de novo extracted signatures from the manuscript as we want to test the presence of our own signatures
cosmic <- cosmic %>% select(SBS7a, SBS9)

# combine all exposure signatures together into one df
sigs <- cbind(sigs_all, cosmic, doxorubicin, cyclo)

# one last check to ensure they are all %
colSums(sigs) # all good

# clean environment
rm(cosmic, cyclo, doxorubicin, sigs_all); gc()

# determine which samples have which exposures (that are not doxorubicin or cyclophosphamide)
# subset exposures table for columns of interest
exposures_subset <- exposures %>% select(Sample, `SBS1+SBS5`, SBSBlood, SBS7a, SBS9)
sample_exposures <- list() # open an empty list to append exposures lists to

for (i in 1:nrow(exposures_subset)) {
  
  # find which signatures are active in each sample
  exp <- exposures_subset[i,2:5] %>% select(where(~ any(.x != 0)))
  
  # add to list
  sample_exposures[[i]] <- colnames(exp)
  
  # name list with sample ID
  names(sample_exposures)[i] <- exposures_subset[i,1]
  
}

# save additional patient info
patient_info <- exposures %>% select(PDid, Sample, cell_type)

rm(exp, exposures_subset, exposures, i);gc() # clean

# now we have all the information we need, we can perform the signature presence tests
# first, make sure patient table and exposures table are in the same order
# same with patient info
colnames(patient) == names(sample_exposures) 
colnames(patient) == patient_info$Sample # yes to both

res_table <- c()
for (i in 1:ncol(patient)) {
  
  # subset patient
  data <- patient[,i] |> as.data.frame()
  colnames(data) <- colnames(patient)[i]
  
  # acquire exposures
  exp1 <- c(sample_exposures[[i]], "SBS-Doxorubicin-1", "SBS-Cyclophosphamide")
  exp2 <- c(sample_exposures[[i]], "SBS-Doxorubicin-2", "SBS-Cyclophosphamide")
  
  # test SBS-Doxor
  res1 <- SignaturePresenceTest(data,
                                as.matrix(sigs[,exp1]),
                                "SBS-Doxorubicin-1",
                                seed = 1234,
                                mc.cores = 2)
  
  res2 <- SignaturePresenceTest(data,
                                as.matrix(sigs[,exp2]),
                                "SBS-Doxorubicin-2",
                                seed = 1234,
                                mc.cores = 2)
  
  # save results to table
  r <- c(patient_info$PDid[i], colnames(patient)[i], patient_info$cell_type[i], signif(res1[[1]][[4]], 3), signif(res2[[1]][[4]], 3))
  res_table <- rbind(res_table, r)
  
}

# convert results table to df
res_table <- res_table |> as.data.frame()

# add column names
colnames(res_table) <- c("Patient ID", "Sample ID", "Cell type", "SBS-Doxorubicin-1 mSigAct p-value", "SBS-Doxorubicin-2 mSigAct p-value")

# save
write_tsv(x = res_table, file = "../data/tables/Table2.tsv")

# clean environment
rm(data, patient, patient_info, res1, res2, sample_exposures, sigs, exp1, exp2, i, r, res_table); gc()
