#!/bin/bash
#SBATCH --clusters=arc
#SBATCH --job-name=Constrained
#SBATCH --output=Constrained.out
#SBATCH --error=Constrained.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00

# For comparison to LSD2, treePL can be run as a penalised likelihood approach to divergence estimation
# There are numerous steps to take for treePL
# Firstly, we need to generate 100 bootstrap trees in IQ-TREE help generate confident intervals later on

# STEP 1:
# Define our input file directories

Alignment=/data/AMAS
Tree=/data/LSD2
Model=/data/IQ-Tree/ModelSearch
Output=/data/treePL
mkdir -p $Output

# STEP 2:
# Load IQ-TREE and run
# Remember to specify bootstraps with -b. In this case we generate 100

ml IQ-TREE/2.2.2.6-gompi-2022b

iqtree2 -s $Alignment/concatenated.fasta -p $Model/MF.best_scheme.nex -m $Model/MF.best_scheme.nex -g $Tree/rooted.nwk -b 100 -bnni -pre $Output/constrained_bs -nt 32
