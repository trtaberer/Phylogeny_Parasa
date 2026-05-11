#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --array=1-67:1
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc_%j.log
#SBATCH --error=fastqc_%j.error

# Set sample name and directory 
# This pipeline requires the number of jobs to match the number of lines in file sample_listFastQC.txt
# Change array number to number of samples being processed in the sample list

# STEP 1:
# Firstly, create a sample list with two columns, where spaces represent columns e.g.
# SAMPLE_NAME SAMPLE_DIRECTORY
# The SAMPLE_NAME is the short form name given to each sample e.g. T1_041
# The SAMPLE_DIRECTORY is the path to the folder where the raw reads have been stored on the system
# Then, define where this sample list is found:

Sample_List=/data/sample_listFastQC.txt

# STEP 2:
# Use slurm array task ID to allocate sample name and directory

SAMPLE_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $Sample_List)
SAMPLE_NAME=$(echo $SAMPLE_LINE | awk '{print $1}')
SAMPLE_DIRECTORY=$(echo $SAMPLE_LINE | awk '{print $2}')

# STEP 3:
# Create a specified output directory. This is where the QC reports will end up
# The -p command of mkdir will create the output folder automatically

Output=/data/FastQC/$SAMPLE_NAME/
mkdir -p $Output

# STEP 4:
# Load and run FastQC

module load FastQC/0.11.9-Java-11

# Run FastQC using specifications from the sample list
# Both the forward (*_1.fq.gz) and reverse (*_2.fq.gz) can be included in the same line
# Specify the output as above

fastqc ${SAMPLE_DIRECTORY}/${SAMPLE_NAME}_*_1.fq.gz ${SAMPLE_DIRECTORY}/${SAMPLE_NAME}_*_2.fq.gz -o $Output

# Finally, move into the output folder, and copy the html file created into a new folder
# This will mean that all html files in the array should end up here, making them easier to download and view

cd ${Output}
mkdir -p /data/FastQC/FastQC_html
scp ${SAMPLE_NAME}*.html /data/FastQC/FastQC_html
