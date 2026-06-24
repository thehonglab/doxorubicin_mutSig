# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for figures presented in Figure 2 of manuscript

# load libraries
library(sigProcessor)
library(tidyverse)
library(ggplotify)

# set working directory
setwd("~/doxorubicin_mutSig/data/")

## Figure 2A
  # created in biorender

## Figure 2B
# load doxorubicin SBS mutational matrices
doxorubicin <- process_mutation_matrix(c("./SBS96/doxorubicin/BT16_doxorubicin.SBS96.all",
                                         "./SBS96/doxorubicin/CHLA-06-ATRT_doxorubicin.SBS96.all",
                                         "./SBS96/doxorubicin/PEDS0005_T1_doxorubicin.SBS96.all",
                                         "./SBS96/doxorubicin/PEDS9001_T1_doxorubicin.SBS96.all",
                                         "./SBS96/doxorubicin/HepG2_doxorubicin_clone1.SBS96.all",
                                         "./SBS96/doxorubicin/HepG2_doxorubicin_clone2.SBS96.all",
                                         "./SBS96/doxorubicin/MCF10A_doxorubicin_clone1.SBS96.all",
                                         "./SBS96/doxorubicin/MCF10A_doxorubicin_clone2.SBS96.all"))

# convert to percentages
doxorubicin_percent <- freq2percent(doxorubicin)

# set a seed
set.seed(1234)

# scale the data
doxorubicin_scaled <- scale(doxorubicin_percent)

# build distance matrix
dist_mat <- dist(t(doxorubicin_scaled), method = "euclidean") # transpose to cluster samples

# build dendrogram
hclust_dend <- hclust(dist_mat, method = 'complete')
plot(hclust_dend) 

# save plot
p <- as.ggplot(~plot(hclust_dend))
ggsave(plot = p, "../figures/Fig2/Figure_2B.png", height = 6, width = 6)

# clean environment
rm(doxorubicin, doxorubicin_scaled, doxorubicin_percent, hclust_dend, dist_mat, p); gc()

## Figure 2C
# import doxorubicin SBS96 spectra
doxorubicin <- process_mutation_matrix(c("./SBS96/doxorubicin/BT16_doxorubicin.SBS96.all",
                                         "./SBS96/doxorubicin/CHLA-06-ATRT_doxorubicin.SBS96.all",
                                         "./SBS96/doxorubicin/PEDS0005_T1_doxorubicin.SBS96.all",
                                         "./SBS96/doxorubicin/PEDS9001_T1_doxorubicin.SBS96.all",
                                         "./SBS96/doxorubicin/HepG2_doxorubicin_clone1.SBS96.all",
                                         "./SBS96/doxorubicin/HepG2_doxorubicin_clone2.SBS96.all",
                                         "./SBS96/doxorubicin/MCF10A_doxorubicin_clone1.SBS96.all",
                                         "./SBS96/doxorubicin/MCF10A_doxorubicin_clone2.SBS96.all"))

# create two data frames that represent our two groups from the cluster dendrogram
doxo1 <- doxorubicin %>% dplyr::select(BT16_doxorubicin, PEDS0005_T1_doxorubicin, HepG2_doxorubicin_clone2, MCF10A_doxorubicin_clone1, MCF10A_doxorubicin_clone2)
doxo2 <- doxorubicin %>% dplyr::select(`CHLA-06-ATRT_doxorubicin`, PEDS9001_T1_doxorubicin, HepG2_doxorubicin_clone1)

# average the rows to get a representative signature for both groups
doxo1_av <- rowMeans(doxo1) |> as.data.frame()
doxo2_av <- rowMeans(doxo2) |> as.data.frame()

# change column names
colnames(doxo1_av) <- "SBS-Doxorubicin-1"
colnames(doxo2_av) <- "SBS-Doxorubicin-2"

# add back mutation types column
doxo1_av$MutationType <- rownames(doxorubicin)
doxo2_av$MutationType <- rownames(doxorubicin)

# remove rownames
rownames(doxo1_av) <- NULL
rownames(doxo2_av) <- NULL

# rearrange columns
doxo1_av <- doxo1_av %>% dplyr::select(MutationType, `SBS-Doxorubicin-1`)
doxo2_av <- doxo2_av %>% dplyr::select(MutationType, `SBS-Doxorubicin-2`)

# save tables 
write_tsv(doxo1_av, "./SBS96/doxorubicin/SBS-Doxorubicin-1.SBS96.all")
write_tsv(doxo2_av, "./SBS96/doxorubicin/SBS-Doxorubicin-2.SBS96.all")

# clean environment
rm(doxo1, doxo1_av, doxo2, doxo2_av, doxorubicin); gc()

# input these tables into SigProfilerPlotting (see python code)

## Figure 2D
# import doxorubicin SBS6 spectra 
doxorubicin <- process_mutation_matrix(c("./SBS6/doxorubicin/BT16_doxorubicin.SBS6.all",
                                         "./SBS6/doxorubicin/CHLA-06-ATRT_doxorubicin.SBS6.all",
                                         "./SBS6/doxorubicin/PEDS0005_T1_doxorubicin.SBS6.all",
                                         "./SBS6/doxorubicin/PEDS9001_T1_doxorubicin.SBS6.all",
                                         "./SBS6/doxorubicin/HepG2_doxorubicin_clone1.SBS6.all",
                                         "./SBS6/doxorubicin/HepG2_doxorubicin_clone2.SBS6.all",
                                         "./SBS6/doxorubicin/MCF10A_doxorubicin_clone1.SBS6.all",
                                         "./SBS6/doxorubicin/MCF10A_doxorubicin_clone2.SBS6.all"))

# create two data frames that represent our two groups from the cluster dendrogram
doxo1 <- doxorubicin %>% dplyr::select(BT16_doxorubicin, PEDS0005_T1_doxorubicin, HepG2_doxorubicin_clone2, MCF10A_doxorubicin_clone1, MCF10A_doxorubicin_clone2)
doxo2 <- doxorubicin %>% dplyr::select(`CHLA-06-ATRT_doxorubicin`, PEDS9001_T1_doxorubicin, HepG2_doxorubicin_clone1)

# convert to %
doxo1_percent <- freq2percent(doxo1)
doxo2_percent <- freq2percent(doxo2)

# average the rows to get a representative signature for both groups
doxo1_av <- rowMeans(doxo1_percent) |> as.data.frame()
doxo2_av <- rowMeans(doxo2_percent) |> as.data.frame()

# change column names
colnames(doxo1_av) <- "SBS-Doxorubicin-1"
colnames(doxo2_av) <- "SBS-Doxorubicin-2"

# combine into one data frame
doxo_percent <- cbind(doxo1_av, doxo2_av)

# convert each sample to percentages as well
sample_percent <- freq2percent(doxorubicin)

# add back mutation types for plotting
doxo_percent$SBS6 <- rownames(doxo_percent) 
sample_percent$SBS6 <- rownames(sample_percent)

# ensure rownames match and then remove rownames
rownames(doxo_percent) == rownames(sample_percent)
doxo_percent$SBS6 == sample_percent$SBS6
rownames(doxo_percent) <- NULL
rownames(sample_percent) <- NULL

# pivot the df
doxo_long <- pivot_longer(data = doxo_percent, cols = c("SBS-Doxorubicin-1", "SBS-Doxorubicin-2"))
sample_long <- pivot_longer(data = sample_percent, cols = colnames(sample_percent)[1:8])

# reverse factors for plotting
doxo_long$SBS6 <- factor(doxo_long$SBS6, levels = rev(unique(doxo_long$SBS6)))
sample_long$SBS6 <- factor(sample_long$SBS6, levels = rev(unique(sample_long$SBS6)))

# for sample_long, assign avg signature each sample belongs to
sample_long$sig <- ifelse(sample_long$name %in% c(colnames(sample_percent)[c(1,3,6,7,8)]), "SBS-Doxorubicin-1", "SBS-Doxorubicin-2")

# ensure signature columns for each df are factors
doxo_long$name <- factor(doxo_long$name, levels = c("SBS-Doxorubicin-2", "SBS-Doxorubicin-1"))
sample_long$sig <- factor(sample_long$sig, levels = c("SBS-Doxorubicin-2", "SBS-Doxorubicin-1"))

# plotting
p <- ggplot(doxo_long, aes(x = (value*100), y = SBS6, fill = name)) +
  geom_bar(position="dodge", stat="identity") +
  geom_point(data = sample_long, aes(x = (value*100), y = SBS6, group = sig, color = sig), inherit.aes = FALSE,
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9),
             size = 2.5) +
  ylab("SBS6") +
  xlab("% Single Nucleotide Substitutions") +
  scale_fill_manual(values = c("#f8931d", "#21aae1"), name = NULL,
                    breaks = c("SBS-Doxorubicin-1", "SBS-Doxorubicin-2")) +
  scale_color_manual(values = c("SBS-Doxorubicin-1" = "#7c490e", 
                                "SBS-Doxorubicin-2" = "#105570"), guide = "none") + 
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"), 
        axis.ticks.x = element_line(colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        axis.text = element_text(size=12)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,35), breaks = seq(0, 35, 5))

print(p)

# save plot
ggsave(plot = p, filename = "../figures/Fig2/Figure_2D.png", height = 6, width = 8)

# clean environment
rm(doxo_av, doxo_long, doxo_percent, doxo1, doxo1_av, doxo2, doxo2_av, doxorubicin, sample_long, sample_percent, p); gc()

# statistical testing performed in GraphPad Prism v11.0.2

## Figure 2E
# combine both average signatures into one table 
df <- process_mutation_matrix(c("./SBS96/doxorubicin/SBS-Doxorubicin-1.SBS96.all",
                                "./SBS96/doxorubicin/SBS-Doxorubicin-2.SBS96.all"))

# add back MutationType column
df$MutationType <- rownames(df)
df_final <- df[,c(3,1:2)]

# save
write_tsv(df_final, "./SigProfilerAssignment/input/average_doxorubicin.SBS96.all")

# clean 
rm(df, df_final); gc()

# input table into SigProfilerAssignment (see python code)

# import the activities table
data <- read.table("./SigProfilerAssignment/output/average_doxorubicin_activities.txt", header = T)

# remove columns that are all zero
data <- data %>% select(where(~ any(.x != 0)))

# get data table of %
data_table <- data
rownames(data_table) <- data_table$Samples
data_table$Samples <- NULL
data_table <- data_table |> t() |> as.data.frame()
data_table <- freq2percent(data_table)

# reformat df
data_percent <- data_table |> t() |> as.data.frame()
data_percent$Samples <- colnames(data_table)
rownames(data_percent) <- NULL
data_percent <- data_percent[,c(10, 1:9)]

# pivot wider
data_longer <- pivot_longer(data_percent, cols = colnames(data)[2:10])

# rename columns names to something more informative
colnames(data_longer)[2:3] <- c("Signature", "Activities")

# specify order of samples
samples_order <- data_longer$Samples |> unique()
data_longer$Samples <- factor(data_longer$Samples, levels = samples_order)

data_longer$Signature <- forcats::fct_reorder(data_longer$Signature, data_longer$Activities)

# plot
p <- ggplot(data_longer, aes(x = Samples, y = Activities, fill = Signature)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = "white")) +
  ylab("Percentage of mutations in each signature") +
  scale_fill_discrete(breaks = c("SBS1", 
                                 "SBS2",
                                 "SBS5", 
                                 "SBS13",
                                 "SBS18",
                                 "SBS37",
                                 "SBS40b",
                                 "SBS100",
                                 "SBS102"), 
                      palette = c("#f782ac",
                                  "#6cab57", 
                                  "#e7703d",
                                  "#6d63a2",
                                  "#cac9c9",
                                  "#efdd5c",
                                  "#882255",
                                  "#66573a",
                                  "#73c4f4")) +
  xlab("")

print(p)

ggsave(filename = "../figures/Fig2/Figure_2E.png", p, width = 4, height = 5)

# clean
rm(data, data_longer, data_percent, data_table, p, samples_order); gc()
