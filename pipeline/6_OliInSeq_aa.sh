#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --array=1-5081:1
#SBATCH --job-name=OliInSeq_aa

# Next, we need to use OliInSeq with both the amino acid (aa) and nucleotide (nt) alignments to ensure they are corresponding
# Check the number of alignments in each of the folders buscophy/results/phylo/busco_aa_alignments_trim and buscophy/results/phylo/busco_nt_alignments_trim
# Go into folder and do
# ls -1 | wc -l

# This will show how many gene alignment files buscophy extracted from the genomes
# Note that this will depend on the number of genes in the targeted lineage of interest as well as success of identifying genes
# As ever, change the --array to this number!

# STEP 1:
# Create a sample list of aa alignment files

Sample_List=/data/sample_listOliInSeq_aa
# The sample list should have 2 columns, where spaces represent columns e.g.
# SAMPLE_NAME SAMPLE_DIRECTORY
# 9918at7088 /buscophy/results/phylo/busco_aa_alignments

# STEP 2:
# Use slurm array task ID to allocate sample name and directory
SAMPLE_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $Sample_List)
SAMPLE_NAME=$(echo $SAMPLE_LINE | awk '{print $1}')
SAMPLE_DIRECTORY=$(echo $SAMPLE_LINE | awk '{print $2}')

# STEP 3:
# Define an output path
OUTPUT=/data/buscophy/results/OliInSeq/aa_filtered
mkdir -p $Output 

# STEP 4:
# Run OliInSeq:
OliInSeq-v0.9.6 -i ${SAMPLE_DIRECTORY}/${SAMPLE_NAME}_aa.fasta -o $OUTPUT/${SAMPLE_NAME}_aa_filtered.fasta
