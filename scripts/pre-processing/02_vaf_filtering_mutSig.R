# 04-01-2026
# filtering tumour-normal vcf
# filtering rules:
  # allele frequency:
    # variant in the exposed sample should have a VAF of > 0.2
    # variants in the parental/control should be removed completely

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

# define sample names for untreated and exposed samples
untreated <- colnames(vcf)[10]
exposed <- colnames(vcf)[11]

# change the column names for untreated and exposed samples to allow processing 
colnames(vcf)[10] <- "UNTREATED"
colnames(vcf)[11] <- "EXPOSED"

# parse out the allele frequencies
vcf_split <- vcf %>% mutate(
  # split up the keys corresponding to UNTREATED and EXPOSED columns
  FORMAT_keys = str_split(FORMAT, ":"), 
  
  # split up the variant info in UNTREATED and EXPOSED columns
  FORMAT_UNTREATED_vals = str_split(UNTREATED, ":"), 
  FORMAT_EXPOSED_vals = str_split(EXPOSED, ":"),
  
  # select the allelic depth information 
  # .x == "AF", finds where AF is in FORMAT_keys (output as TRUE/FALSE)
  # which, converts TRUE/FALSE to an index
  # .y[2], grabs the corresponding value from FORMAT_vals (either UNTREATED or EXPOSED)
  UNTREATED_AF = map2_chr(FORMAT_keys, FORMAT_UNTREATED_vals, ~ .y[which(.x == "AF")]),
  EXPOSED_AF = map2_chr(FORMAT_keys, FORMAT_EXPOSED_vals, ~ .y[which(.x == "AF")])
  
) %>% dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, UNTREATED, EXPOSED, UNTREATED_AF, EXPOSED_AF) # select all original columns as well as the new AF columns

# filter for variants that are present at a VAF of >= 20% in exposed but not at all present in parental/untreated
vcf_filtered <- vcf_split %>% filter(UNTREATED_AF == 0 & EXPOSED_AF >= 0.2)

# remove the AF columns as we don't need those anymore
vcf_final <- vcf_filtered %>% dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, UNTREATED, EXPOSED)

# rename the UNTREATED and EXPOSED columns to their original names
colnames(vcf_final)[10] <- untreated
colnames(vcf_final)[11] <- exposed

# clean environment
rm(vcf, vcf_filtered, vcf_split, exposed, untreated); gc()

# save as a tsv to be converted to a vcf
write.table(x = vcf_final, file = paste(output, "_var_tab.tsv", sep = ""),
            col.names = F, row.names = F, quote = F, sep = "\t")

