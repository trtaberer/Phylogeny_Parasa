#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --array=1-2024:1
#SBATCH --job-name=Gene_tree
#SBATCH --output=Gene_tree_%j.log
#SBATCH --error=Gene_tree_%j.error

# As a second method for comparison, we want to generate individual trees from each gene alignment

# We use the alignments that were filtered with OliInSeq, trimAl and Geneious, saved in /data/filtered_alignments
# Change array from 1-[number of alignments in list]

# STEP 1:
# We will need to create a text file specifying the name of samples

Sample_List=/data/sample_listGenetreealignments.txt
# The sample list should have 1 column:
# GENE_NAME
# 1000at7088_nt_filtered_trimmed

# STEP 2:
# Use slurm array task ID to allocate sample name and directory
GENE_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $Sample_List)
GENE_NAME=$(echo $GENE_LINE | awk '{print $1}')

# STEP 3:
# Define input and output directories
Input=/data/filtered_alignments
Output=/data/Gene_trees
mkdir -p $Output

# STEP 4:
# Load IQ-Tree and define full path to gene alignment

ml IQ-TREE/2.2.2.6-gompi-2022b

GENE_FILE="$Input/$GENE_NAME.fasta"

# STEP 5:
# Run IQ-TREE with 5 likelihood searches (-runs 5)

iqtree2 -s "$GENE_FILE" -m MFP -bb 1000 -alrt 1000 --runs 5 -nt 32 -st DNA -pre "$Output/$GENE_NAME"
