# Derivation of a doxorubicin mutational signatures in cancer

This repository contains all code and data associated with the manuscript "Derivation of doxorubicin mutational signatures in cancer", Skinner et al., 2026, submitted to Genome Research. 

## Repository contents:

1. Data:
- Mutational matrices for the specified groups of samples:
	- SBS6 (doxorubicin)
	- SBS96 (background, cisplatin, doxorubicin, patient profiles, and patient exposures)
	- DBS78 (cisplatin and doxorubicin)
	- ID30 (doxorubicin)

- SigProfilerAssignment inputs:
	- Table of SBS96 background mutational profiles
	- Table of SBS96 averaged doxorubicin profiles

- Table data (main and supplemental)


2. Scripts:
- Variant filtering (sample pre-processing)

- R scripts for plotting:
	- Cosine similarities
	- Hierarchical clustering
	- SBS6 bar charts
	- SigProfilerAssignment activities output
	- ID30 mutational profiles

- Python scripts for plotting:
	- SBS96 mutational profiles
	- DBS78 mutational profiles

- R scripts for generating all data in tables (main and supplemental)

- Additional functions used for processing and plotting of INDEL data
