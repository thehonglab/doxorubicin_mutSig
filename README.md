# Derivation of a doxorubicin mutational signatures in cancer

This repository contains all code and data associated with the manuscript "Derivation of doxorubicin mutational signatures in cancer", Skinner et al., 2026, submitted to Genome Research. 

## Repository contents

#### Data
- Mutational matrices for the specified groups of samples:
	- SBS6 (doxorubicin)
	- SBS96 (background, cisplatin, doxorubicin, patient profiles, and patient exposures)
	- DBS78 (cisplatin and doxorubicin)
	- ID30 (doxorubicin)

- SigProfilerAssignment inputs:
	- Table of SBS96 background mutational profiles
	- Table of SBS96 averaged doxorubicin profiles

- Table data (main and supplemental)

#### Scripts
- Variant filtering and mutational matrices generation (scripts/pre-processing)

- R scripts for plotting (scripts/figures):
	- Cosine similarities
	- Hierarchical clustering
	- SBS6 bar charts
	- SigProfilerAssignment activities output
	- ID30 mutational profiles

- Python scripts for plotting (scripts/figures):
	- SBS96 mutational profiles
	- DBS78 mutational profiles

- R scripts for generating all data in tables (main and supplemental) (scripts/tables)

- Additional functions used for processing and plotting of INDEL data (scripts/additional_functions)

## Processing raw data

#### Genome alignment 

We aligned FASTQ files to GrCh38p13 using Illumina DRAGEN v4.4.4.f2:

```
dragen -f -r GrCh38p13 -1 WGS/SAMPLE_1.fq.gz \
-2 WGS/SAMPLE_2.fq.gz \
--RGID SAMPLE \
--RGSM SAMPLE_1 \
--output-directory WGSOut \
--output-file-prefix SAMPLE \
--enable-map-align true \
--enable-duplicate-marking true \
--enable-map-align-output true \
--enable-cnv false \
--enable-variant-caller false \
--validate-pangenome-reference=false
```

BAM files output from the above step have been deposited to the database of Genotypes and Phenotypes (dbGaP) under the accession number: phs004740.v1.p1

#### Variant calling

Next, using the BAM files as input, we performed variant calling. We did this is two ways:

1. Tumour-only, to characterise the background mutational profiles of all parental and non-treated samples:

```
dragen -f -r GrCh38p13 -b WGSOut/CONTROL.bam \
--enable-variant-caller true \
--output-directory VCFs \
--output-file-prefix CONTROL \
--enable-map-align false
```

2. Tumour-normal, where:
- The drug exposed is the tumour 
- The parental/non-treated is the normal

This removes mutations found in the control samples (i.e., background) to isolate the mutations caused by the drug exposure

```
dragen -f -r GrCh38p13 --tumor-bam-input WGSOut/EXPOSED.bam \
-b WGSOut/CONTROL.bam \
--enable-variant-caller true \
--output-directory VCFs \
--output-file-prefix SAMPLE_withCONTROL \
--enable-map-align false
```

Hard-filtered VCFs output from the above two steps have been deposited to the database of Genotypes and Phenotypes (dbGaP) under the accession number: phs004740.v1.p1

These VCFs were used as input to our pre-processing pipeline: scripts/pre-processing

## Packages

Listed below are R and python packages needed to complete data analysis, along with installation instructions.

#### R packages:
- vcfR
```
install.packages('vcfR')
```

- tidyverse
```
install.packages("tidyverse")
```

- argparser
```
install.packages("argparse")
```

- sigProcessor
```
install.packages("remotes")
remotes::install_github("katieskinner98/sigProcessor")
```

- ggplotify
```
install.packages("ggplotify")
```

- indelsig.tools.lib
```
install.packages("devtools")
devtools::install_github("Nik-Zainal-Group/indelsig.tools.lib")
```

- data.table
```
install.packages("data.table")
```

- mSigAct
```
install.packages("remotes")
remotes::install_github("steverozen/mSigAct")

- mSigTools
```
install.packages("mSigTools")
```

- cosmicsig
```
install.packages("cosmicsig")
```
 
#### Python packages:

- SigProfilerMatrixGenerator
```
pip install SigProfilerMatrixGenerator
```

- SigProfilerPlotting
```
pip install SigProfilerPlotting
```

- SigProfilerAssignment
```
pip install SigProfilerAssignment
```

- SigProfilerTopography
```
pip install SigProfilerTopography
```
