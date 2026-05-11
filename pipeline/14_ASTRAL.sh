#!/bin/bash
#SBATCH --job-name=ASTRAL
#SBATCH --output=ASTRAL.out
#SBATCH --error=ASTRAL.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00

# Now, we can generate a consensus tree with ASTRAL-IV based on our individual gene trees

# STEP 1:
# After IQ-TREE loop has run and generated a tree for each gene, combine them all into one file:
# cat /data/Gene_trees/*.treefile > all_gene_trees.tre

# STEP 2:
# Load Java
module load Java/11.0.2

# STEP 3:
# Specify input and output folders
INPUT=/data/Gene_trees
OUTPUT=/data/ASTRAL
mkdir -p $Output

# STEP 4:
# Run ASTRAL-IV, ensuring you specify the outgroup taxon with --root

astral4 -t 32 -i $INPUT/all_gene_trees.tre --root Taeda_prasina_T1075 -o $OUTPUT/Parasa_genetree_output.tre 2> $OUTPUT/Parasa_genetree.log

# Supports multithreading with -t. Make sure this matches number of cores at top of SLURM script
# Can be run up to 64 cores; best to use at least 16 for big jobs
