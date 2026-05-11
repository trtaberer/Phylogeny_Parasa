#!/bin/bash
#SBATCH --job-name=IQ-TREE_TS
#SBATCH --output=IQ-TREE_TS.out
#SBATCH --error=IQ-TREE_TS.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=14-00:00:00

# After the model search, it's time to run the tree search on the concatenated alignment.

# STEP 1:
# Define the input and output folders

Input=/data/AMAS
Model=/data/IQ-Tree/ModelSearch
Output=/data/IQ-Tree/TreeSearch
mkdir -p $Output

# STEP 2:
# Load and run IQ-TREE tree search with 1000 ultrafast bootstrap replicates (-bb 1000) and approximate likelihood-ratio test (-alrt 1000)

ml IQ-TREE/2.2.2.6-gompi-2022b

iqtree2 -s $Input/concatenated.fasta -p $Model/MF.best_scheme.nex -m $Model/MF.best_scheme.nex -bb 1000 -alrt 1000 -bnni --runs 10 -nt 32 -st DNA -pre $Output/TS


