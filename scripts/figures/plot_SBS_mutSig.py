# 07-21-2025
# plotting of mutational signatures from generated SBS mutational matricies

# import packages
import importlib
import subprocess
import sys
import os
from pathlib import Path
import argparse

# list of package names
packages = [
    'SigProfilerPlotting'
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
import sigProfilerPlotting as sigPlot

# print the version of SigProfilerMatrixGenerator
print("Version of sigProfilerPlotting is:", sigPlot.__version__)

# set up our argparse module
parser=argparse.ArgumentParser(description ='mutational signatures matrix generation')

# add arguments
parser.add_argument('--sample_dir', '-d', type=str, help='file path of the sample to be processed', action="store")
parser.add_argument('--top_dir', '-t', type=str, help='top level file path to samples to be processed (if processing multiple)', action="store")
parser.add_argument('--plot_dir', '-p', type=str, help='file path where plots should be saved', action="store")

# parse the inputted arguments
args = parser.parse_args()

# define the file path to the sample matrix
matrix_path=args.top_dir + "/" + args.sample_dir + "/final_filtered/output/SBS/"
out_path=args.plot_dir

# acquire the sample name in the correct format (for ease, we will be using the sample_dir, but you can parse in the sample name if you'd like)
# for some samples, there will be a main directory for cell cell, then a sub-directory for clone or mouse PDX number
# so we want to acquire the sub-directory name for samples where this is the case
id=Path(args.sample_dir)
id_final=os.path.basename(id)

# plotting SBS6
sigPlot.plotSBS((matrix_path + id_final + ".SBS6.all"), (out_path + "percent/SBS6/"), id_final, "6", percentage=True)

sigPlot.plotSBS((matrix_path + id_final + ".SBS6.all"), (out_path + "counts/SBS6/"), id_final, "6", percentage=False)

# plotting SBS24
sigPlot.plotSBS((matrix_path + id_final + ".SBS24.all"), (out_path + "percent/SBS24/"), id_final, "24", percentage=True)

sigPlot.plotSBS((matrix_path + id_final + ".SBS24.all"), (out_path + "counts/SBS24/"), id_final, "24", percentage=False)

# plotting SBS96
sigPlot.plotSBS((matrix_path + id_final + ".SBS96.all"), (out_path + "percent/SBS96/"), id_final, "96", percentage=True)

sigPlot.plotSBS((matrix_path + id_final + ".SBS96.all"), (out_path + "counts/SBS96/"), id_final, "96", percentage=False)
