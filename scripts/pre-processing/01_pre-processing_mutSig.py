# 07-20-2025
# this script is for pre-processing of VCFs
# the first step is to filter for variants that pass QC metrics
# the second step is to remove common population variants, decoy sequences, etc... things that likely won't be the result of mutagen exposure

# import packages
import subprocess
import argparse

# first we want to define the absolute path to the blacklist bed file provided to us by the Rozen Lab
blacklist="/Users/katieskinner98/ResistanceProjectWGS/mutationalSignaturesAnalysis/blacklisted_regions/MasterFilter_GRCh38p7_srt_merge.bed"

# set up our argparse module
parser=argparse.ArgumentParser(description ='VCF pre-processing')

# add arguments
parser.add_argument('--sample_name', '-s', type=str, help='name of your sample', action="store")
parser.add_argument('--sample_dir', '-d', type=str, help='file path of the sample to be processed', action="store")
parser.add_argument('--top_dir', '-t', type=str, help='top level file path to samples to be processed (if processing multiple)', action="store")

# parse the inputted arguments
args = parser.parse_args()

# define the file path to the sample
path=args.top_dir + "/" + args.sample_dir + "/" + args.sample_name

# the first command is to unzip our VCF
subprocess.run(["bgzip", "-d", "-k", (path + ".hard-filtered.vcf.gz")])

# the second command filters for variants that pass QC metrics
input1=f"{path}.hard-filtered.vcf"
output1=f"{path}.PASS.vcf"

with open(output1, "w") as out_fh:
    subprocess.run(
        ["bcftools", "view", "-i", 'FILTER=="PASS"', input1],
        stdout=out_fh,
        check=True
    )

# the third command removes common population variants, decoy sequences, etc... things that likely wouldn't be results of mutagen exposure
input2=output1
output2=f"{path}.blacklistRemoved.hard-filtered.vcf"

with open(output2, "w") as out_fh:
    subprocess.run(
        ["bedtools", "intersect", "-a", input2, "-b", blacklist, "-header", "-v"],
        stdout=out_fh,
        check=True
    )
