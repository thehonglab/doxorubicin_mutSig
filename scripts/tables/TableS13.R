# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for Table S13 of manuscript

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

rm(exp, exposures_subset, exposures, i);gc() # clean

# now we have all the information we need, we can perform the signature presence tests
# first, make sure patient table and exposures table are in the same order
colnames(patient) == names(sample_exposures) # yes

# import Table 2 results
res_table <- read.table("./tables/Table2.tsv", header = T, sep = "\t", check.names = F)

# select significant results
doxo1 <- res_table %>% filter(`SBS-Doxorubicin-1 mSigAct p-value` <= 0.05)
doxo2 <- res_table %>% filter(`SBS-Doxorubicin-2 mSigAct p-value` <= 0.05)

# select the corresponding samples from the exposures list
doxo1_exp <- sample_exposures[doxo1$`Sample ID`]
doxo2_exp <- sample_exposures[doxo2$`Sample ID`]

# double check both table and lists are in the same order
doxo1$`Sample ID` == names(doxo1_exp)
doxo2$`Sample ID` == names(doxo2_exp)

# for each sample we want to check how many of the total mutations can be attributed to each signature
# define full list of all exposures
all_exp <- colnames(sigs)

# SBS-Doxorubicin-1
burden_df <- c()

for (i in 1:nrow(doxo1)) {
  
  # subset patient
  data <- patient %>% select(all_of(doxo1$`Sample ID`[i]))
  
  # acquire exposures
  exp <- c(doxo1_exp[[i]], "SBS-Doxorubicin-1", "SBS-Cyclophosphamide")
  
  # get the burden
  burden <- as.data.frame(optimize_exposure_QP(spectrum = data, signatures = as.matrix(sigs[,exp])) |> round() |> as.data.frame() |> t())
  
  # reset rownames
  rownames(burden) <- NULL
  
  # find missing signatures (for table formatting purposes)
  missing <- setdiff(all_exp[c(1:5,7)], colnames(burden))
  
  # fill those in as zero
  burden[missing] <- 0
  
  # add sample name
  burden$sample <- colnames(data)
  
  # reformat
  burden <- burden %>% select(sample, `SBS1+SBS5`, SBSBlood, SBS7a, SBS9, `SBS-Doxorubicin-1`, `SBS-Cyclophosphamide`)
  
  # append to df
  burden_df <- rbind(burden_df, burden)
  
  rm(burden, data, exp, missing); gc() # clean
  
}

# rename object
burden_df1 <- burden_df

rm(burden_df, i); gc() # clean

# SBS-Doxorubicin-2
burden_df <- c()

for (i in 1:nrow(doxo1)) {
  
  # subset patient
  data <- patient %>% select(all_of(doxo2$`Sample ID`[i]))
  
  # acquire exposures
  exp <- c(doxo2_exp[[i]], "SBS-Doxorubicin-2", "SBS-Cyclophosphamide")
  
  # get the burden
  burden <- as.data.frame(optimize_exposure_QP(spectrum = data, signatures = as.matrix(sigs[,exp])) |> round() |> as.data.frame() |> t())
  
  # reset rownames
  rownames(burden) <- NULL
  
  # find missing signatures (for table formatting purposes)
  missing <- setdiff(all_exp[c(1:4,6,7)], colnames(burden))
  
  # fill those in as zero
  burden[missing] <- 0
  
  # add sample name
  burden$sample <- colnames(data)
  
  # reformat
  burden <- burden %>% select(sample, `SBS1+SBS5`, SBSBlood, SBS7a, SBS9, `SBS-Doxorubicin-2`, `SBS-Cyclophosphamide`)
  
  # append to df
  burden_df <- rbind(burden_df, burden)
  
  rm(burden, data, exp, missing); gc() # clean
  
}

# rename object
burden_df2<- burden_df

rm(burden_df, doxo1, doxo1_exp, doxo2, doxo2_exp, patient, sample_exposures, sigs, all_exp, i); gc() # clean

# reformat table
burden_df1$ID <- (res_table %>% filter(`Sample ID` %in% burden_df1$sample))$`Patient ID`
burden_df2$ID <- (res_table %>% filter(`Sample ID` %in% burden_df2$sample))$`Patient ID`

burden_df1_final <- burden_df1 %>% select(ID, sample, `SBS1+SBS5`, SBSBlood, SBS7a, SBS9, `SBS-Doxorubicin-1`, `SBS-Cyclophosphamide`)
burden_df2_final <- burden_df2 %>% select(ID, sample, `SBS1+SBS5`, SBSBlood, SBS7a, SBS9, `SBS-Doxorubicin-2`, `SBS-Cyclophosphamide`)

colnames(burden_df1_final)[1:2] <- c("Patient ID", "Sample ID")
colnames(burden_df2_final)[1:2] <- c("Patient ID", "Sample ID")

# save
write_tsv(x = as.data.frame(burden_df1_final), file = "./tables/TableS13_1.tsv")
write_tsv(x = burden_df2_final, file = "./tables/TableS13_2.tsv")

# clean environment
rm(burden_df1, burden_df1_final, burden_df2, burden_df2_final, res_table); gc()
