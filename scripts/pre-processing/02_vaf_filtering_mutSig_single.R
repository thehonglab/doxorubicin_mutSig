# 04-01-2026
# filtering tumour only vcf (to derive background profiles)

# load libraries
library(vcfR)
library(tidyverse)
library(argparser)

# define argument parser
vcf_input <- arg_parser(name = "vcf_input", description = "Input path to VCF")

# add input argument
vcf_input <- add_argument(parser = vcf_input, arg = "--input", help = "Input path to VCF")

# add output argument
vcf_input <- add_argument(parser = vcf_input, arg = "--output", help = "Prefix to output file")

# parse arguments
parse <- parse_args(vcf_input)

# define input argument
input <- parse$input

# define output argument
output <- parse$output

# print message to display which sample is being processed
print(paste("Processing:", input))

# read in vcf
vcf <- read.vcfR(input)

# acquire a table of variants with allele frequency information
vcf <- cbind(vcf@fix, vcf@gt) |> as.data.frame()

# define name for sample
sample <- colnames(vcf)[10]

# change the column names for untreated and exposed samples to allow processing 
colnames(vcf)[10] <- "SAMPLE"

# filter for canonical chromosomes
vcf <- vcf %>% filter(CHROM %in% paste0("chr", c(1:22, "X", "Y", "M")))

# remove variants called with targeted caller
targeted <- grepl("TARGETED", vcf$INFO)
vcf_non_targeted <- vcf[!targeted,]

# identify multi-allelic sites (i.e., multiple ALT calls, comma separated)
multi <- grepl(",", vcf_non_targeted$ALT)

# subset out the multi-allelic sites
vcf_single <- vcf_non_targeted[!multi,]

# clean environment
rm(vcf, targeted, vcf_non_targeted, multi); gc()

# parse out the allele frequencies
vcf_split <- vcf_single %>% mutate(
  # split up the keys corresponding to SAMPLE column
  FORMAT_keys = str_split(FORMAT, ":"), 
  
  # split up the variant info in SAMPLE column
  FORMAT_vals = str_split(SAMPLE, ":"), 
  
  # select the allelic depth information 
  # .x == "AF", finds where AF is in FORMAT_keys (output as TRUE/FALSE)
  # which, converts TRUE/FALSE to an index
  # .y[...], grabs the corresponding value from FORMAT_vals 
  AF = map2_chr(FORMAT_keys, FORMAT_vals, ~ .y[which(.x == "AF")])
  
) %>% dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE, AF) # select all original columns as well as the new AF columns

# filter for variants that are present at a VAF of >= 20% in exposed but not at all present in parental/untreated
vcf_filtered <- vcf_split %>% filter(AF >= 0.2)

# remove the AF columns as we don't need those anymore
vcf_final <- vcf_filtered %>% dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE)

# rename the UNTREATED and EXPOSED columns to their original names
colnames(vcf_final)[10] <- sample

# clean environment
rm(vcf_filtered, vcf_single, vcf_split, sample, input); gc()

# save as a tsv to be converted to a vcf
write.table(x = vcf_final, file = paste(output, "_var_tab.tsv", sep = ""),
            col.names = F, row.names = F, quote = F, sep = "\t")

