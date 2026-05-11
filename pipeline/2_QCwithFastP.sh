#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --array=1-67:1
#SBATCH --job-name=FastP
#SBATCH --output=FastP.%A_%a.out
#SBATCH --error=FastP.%A_%a.error

# Set sample name and directory 
# This pipeline requires the number of jobs to match the number of lines in file sample_listFastP.txt
# Change array number to number of samples being processed in the sample list

# STEP 1:
# Firstly, create a sample list with two columns, where spaces represent columns e.g.
# SAMPLE_NAME SAMPLE_DIRECTORY
# The SAMPLE_NAME is the short form name given to each sample e.g. T1_041
# The SAMPLE_DIRECTORY is the path to the folder where the raw reads have been stored on the system

Sample_List=/data/sample_listFastP.txt

# STEP 2:
# Use slurm array task ID to allocate sample name and directory

SAMPLE_NAME=$(cat $Sample_List | head -n $SLURM_ARRAY_TASK_ID | tail -1 | awk {'print $1}')
SAMPLE_DIRECTORY=$(cat $Sample_List | head -n $SLURM_ARRAY_TASK_ID | tail -1 | awk {'print $2}')

# STEP 3:
# Move into sample directory
cd $SAMPLE_DIRECTORY

# STEP 4:
# Running FastP (https://github.com/OpenGene/fastp)
# Install FastP and define its file path e.g.:

FASTP=/data/Software/fastp

# STEP 5:
# Set up for loop to conduct filtering for each read pair
for ReadPair in `ls ${SAMPLE_NAME}_*_1.fq.gz | cut -f1,2,3,4,5 -d'_'`
do

  #Use FastP to conduct automated filtering of fastq files
  #Note: based on initial test we will trim the first 10bp from start of each read
 $FASTP \
 -i ${ReadPair}_1.fq.gz \
 -o Filtered_${ReadPair}_1.fq.gz \
 -I ${ReadPair}_2.fq.gz \
 -O Filtered_${ReadPair}_2.fq.gz \
 --trim_front1 10 \
 --trim_front2 10

# STEP 6: 
Rename QC reports and move to FastP report folder

mkdir /data/FastP

mv fastp.html /data/FastP/${ReadPair}.html
mv fastp.json /data/FastP/${ReadPair}.json

done
