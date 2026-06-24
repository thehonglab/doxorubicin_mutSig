# Skinner, KT., et al, 2026. "Derivation of a doxorubicin mutational signature in cancer". Genome Research
# Contact for questions: Katie Skinner; katie.skinner@emory.edu

# Code for figures presented in Supplemental Figure 2 of manuscript

# load libraries
library(sigProcessor)

# set working directory
setwd("~/doxorubicin_mutSig/data/")

## Figure S2A
  # created in python

## Figure S2B
# load background SBS mutational matrices
background <- process_mutation_matrix(list.files("./SBS96/background/", full.names = T))

# add back MutationType column
background$MutationType <- rownames(background)
background <- background[,c(7,1:6)]
rownames(background) <- NULL

# save table
write_tsv(background, "./SigProfilerAssignment/input/background_matrix.SBS96.all")

rm(background); gc() # clean

# run SigProfilerAssignment on these background profiles (see python code)

# import the signature activities
data <- read.table("./SigProfilerAssignment/output/background_activities.txt", header = T)

# remove columns that are all zero
data <- data %>% select(where(~ any(.x != 0)))

# get data table of %
data_table <- data
rownames(data_table) <- data_table$Samples
data_table$Samples <- NULL
data_table <- data_table |> t() |> as.data.frame()
data_table <- freq2percent(data_table)

# reformat the table
data_percent <- data_table |> t() |> as.data.frame()
data_percent$Samples <- colnames(data_table)
rownames(data_percent) <- NULL
data_percent <- data_percent[,c(4, 1:3)]

# pivot wider
data_longer <- pivot_longer(data_percent, cols = colnames(data)[2:4])

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
  scale_fill_discrete(breaks = c("SBS1", "SBS5", "SBS54"), palette = c("#6cab57", "#f782ac", "#73c4f4"))

p

ggsave(filename = "../figures/FigS2/Figure_S2B.png", p, width = 4, height = 4.5)

# clean
rm(data, data_longer, data_percent, data_table, p, samples_order); gc()
