#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=1-10:00:00
#SBATCH --array=1-67:1
#SBATCH --job-name=SPAdes_assembly
#SBATCH --output=SPAdes_assembly_%j.log
#SBATCH --error=SPAdes_assembly_%j.error

# Following general rules of:
# 1 node per task (i.e. 1 assembly/array)
# High amount of RAM due to genome sizes
# Omitted [SBATCH --nodes=1] as the scheduler should be able to work this out based on other requirements
# cpus-per-tak=32 as default for SPAdes - this can be increased to 64
# memory=128G as default for SPAdes - this can be increased i.e. doubled
# Time can be adjusted based on how long assembler takes for each sample; 10 hours seems safe
# Amount of time should be monitored and changed accordingly
# Amount of time is per assembly, NOT for the whole set of arrays

# Set sample name and directory 
# This pipeline requires the number of jobs to match the number of lines in file sample_listSPAdes.txt
# Change array number to number of samples being processed in the sample list

# STEP 1:
# We will need to create a text file specifying the name of samples
# We want to process filtered reads are stored following FASTP

Sample_List=/data/sample_listSPAdes.txt
# The sample list should have 3 columns, where spaces represent columns e.g.
# SAMPLE_HEAD SAMPLE_NAME SAMPLE_DIRECTORY
# T1_044 Filtered_T1_044 path/to/filtered/sample
# T1_044 Filtered_T1_044 path/to/filtered/sample

# STEP 2:
# Use slurm array task ID to allocate sample name and directory
SAMPLE_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $Sample_List)
SAMPLE_HEAD=$(echo $SAMPLE_LINE | awk '{print $1}')
SAMPLE_NAME=$(echo $SAMPLE_LINE | awk '{print $2}')
SAMPLE_DIRECTORY=$(echo $SAMPLE_LINE | awk '{print $3}')

# STEP 3:
# Define a path to the output folder.

Output=/data/SPAdes/$SAMPLE_HEAD/
mkdir -p $Output

# STEP 4:
# Load Anaconda in a virtual environment in order to run spades.py
module load Anaconda3/2024.02-1
source activate $DATA/myenv

# STEP 5:
# Run SPAdes
# Ensure memory and threads (i.e. cpus per task) are equal to or less than what is requested in SLURM list above

spades.py --isolate --pe1-1 ${SAMPLE_DIRECTORY}/${SAMPLE_NAME}_*_1.fq.gz --pe1-2 ${SAMPLE_DIRECTORY}/${SAMPLE_NAME}_*_2.fq.gz -o $Output --threads 32 --memory 128

# STEP 6:
# Run quast for quality check of assembly
# Define an output directory for quast of contigs

Output_Contigs=/data/SPAdes/${SAMPLE_HEAD}/${SAMPLE_HEAD}_quast_output_contigs
mkdir -p $Output_Contigs

# STEP 7:
# Run quast on contigs with the following

quast.py ${Output}/contigs.fasta -o $Output_Contigs --threads 4

# Rename resulting log file and copy to QUAST_output_contigs folder for easy access

cd ${Output_Contigs}
mv report.txt ${SAMPLE_HEAD}_report.txt
mkdir -p /data/SPAdes/QUAST_output_contigs
scp ${SAMPLE_HEAD}_report.txt /data/SPAdes/QUAST_output_contigs


# STEP 8:
# Repeat quast for scaffolds
# Define an output directory for quast on scaffolds

Output_Scaffolds=/data/SPAdes/${SAMPLE_HEAD}/${SAMPLE_HEAD}_quast_output_scaffolds
mkdir -p $Output_Scaffolds

# Run quast with the following

quast.py ${Output}/scaffolds.fasta -o $Output_Scaffolds --threads 4

# STEP 9:
# Rename resulting log file and copy to QUAST_output_scaffolds folder for easy access

cd ${Output_Scaffolds}
mv report.txt ${SAMPLE_HEAD}_report.txt
mkdir -p /data/SPAdes/QUAST_output_scaffolds
scp ${SAMPLE_HEAD}_report.txt /data/SPAdes/QUAST_output_scaffolds

