# 07-21-2025
# generation of mutational matricies that we will pipe into the plot function

# import packages
import importlib
import subprocess
import sys
import os
from pathlib import Path
import argparse

# list of package names
packages = [
    'SigProfilerMatrixGenerator'
]

# check if packages are installed, if not, install them
for package in packages:
    try:
        importlib.import_module(package)
        print(f"✅ {package} is already installed.")
    except ImportError:
        print(f"📦 {package} not found. Installing...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# import packages
import SigProfilerMatrixGenerator
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

# print the version of SigProfilerMatrixGenerator
print("Version of SigProfilerMatrixGenerator is:", SigProfilerMatrixGenerator.__version__)

# set up our argparse module
parser=argparse.ArgumentParser(description ='mutational signatures matrix generation')

# add arguments
parser.add_argument('--sample_dir', '-d', type=str, help='file path of the sample to be processed', action="store")
parser.add_argument('--top_dir', '-t', type=str, help='top level file path to samples to be processed (if processing multiple)', action="store")

# parse the inputted arguments
args = parser.parse_args()

# define the file path to the sample
out_path=args.top_dir + "/" + args.sample_dir + "/" + "final_filtered/"

# acquire the sample name in the correct format (for ease, we will be using the sample_dir, but you can parse in the sample name if you'd like)
# for some samples, there will be a main directory for cell cell, then a sub-directory for clone or mouse PDX number
# so we want to acquire the sub-directory name for samples where this is the case
id=Path(args.sample_dir)
id_final=os.path.basename(id)

# now we have the final ID, we can generate our matrices
matrices = matGen.SigProfilerMatrixGeneratorFunc(id_final,
                                                 "GRCh38",
                                                 out_path,
                                                 plot=False, exome=False, bed_file=None,
                                                 chrom_based=False, tsb_stat=False, seqInfo=True,
                                                 cushion=100)
