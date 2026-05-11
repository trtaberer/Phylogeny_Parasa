#!/bin/bash
#SBATCH --job-name=buscophy
#SBATCH --output=buscophy.out
#SBATCH --error=buscophy.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=14-00:00:00

# This pipeline uses buscophy to extract busco genes from the assembled genomes
# More detail about buscophy can be found here: https://gitlab.leibniz-lib.de/smartin/buscophy

# STEP 1:
# Download and install buscophy

git clone https://gitlab.leibniz-lib.de/smartin/buscophy.git
cd buscophy

# load a module to activate minforge.
# If minforge is installed globally skip the next line.
module load miniforge/24.7.1

conda create --yes --prefix ./buscophy_env --channel bioconda --channel conda-forge snakemake mamba python==3.12

conda activate ./buscophy_env

# STEP 2:
# Run buscophy
# As this project was about a genus of moth, the --lineage flag was chosen as lepidoptera_odb10
# This can be adjusted depending on lineage of interest
# Including the --frag flag tells buscophy to also extract fragmented genes

buscophy --threads 8 --input /data/SPAdes/Renamed --lineage lepidoptera_odb10 --frag

# buscophy will run and create numerous output files within the directory buscophy

