#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --clusters=arc
#SBATCH --time=00:05:00
#SBATCH --array=1-5081:1
#SBATCH --job-name=Trimal_nt

# We now want to remove any columns in the alignments that had more than 50% missing data
# This will again be for both aa and nt files, but using the files we just created with OliInSeq
# For this, we use trimAl

# STEP 1:
# Create a sample list of nt alignment files

Sample_List=/data/sample_listtrimAl_nt.txt
# The sample list should have 2 columns, where spaces represent columns e.g.
# SAMPLE_NAME SAMPLE_DIRECTORY
# 9918at7088 /data/buscophy/results/OliInSeq/nt_filtered

# STEP 2:
# Use slurm array task ID to allocate sample name and directory
SAMPLE_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $Sample_List)
SAMPLE_NAME=$(echo $SAMPLE_LINE | awk '{print $1}')
SAMPLE_DIRECTORY=$(echo $SAMPLE_LINE | awk '{print $2}')

# STEP 3:
# Define an output path
OUTPUT=/data/buscophy/results/Trimal/nt_filtered_trimmed
mkdir -p $Output

# STEP 4:
# Run trimal:
trimal -in ${SAMPLE_DIRECTORY}/${SAMPLE_NAME}_nt_filtered.fasta -out $OUTPUT/${SAMPLE_NAME}_nt_filtered_trimmed.fasta -gt 0.5


