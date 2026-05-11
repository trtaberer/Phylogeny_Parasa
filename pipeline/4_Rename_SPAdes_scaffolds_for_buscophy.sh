#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=1-1:00:00
#SBATCH --array=1-67:1
#SBATCH --job-name=rename_SPAdes_scaffolds
#SBATCH --output=rename_SPAdes_scaffolds_%j.log
#SBATCH --error=rename_SPAdes_scaffolds_%j.error

# Processing necessary before running buscophy

# STEP 1:
# Make a text file with scaffolds, containing the following 4 columns of information:
# TAXONOMY SAMPLE_NAME SAMPLE_NAME_EDIT SAMPLE_DIRECTORY
# e.g.
# Parasa_euchlora T1_054 T1054 /path/to/SPAdes/assembly

Sample_List=/data/sample_listSPAdes_scaffolds.txt

# STEP 2:
# Use slurm array task ID to allocate sample name and directory
SAMPLE_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $Sample_List)
TAXONOMY=$(echo $SAMPLE_LINE | awk '{print $1}')
SAMPLE_NAME=$(echo $SAMPLE_LINE | awk '{print $2}')
SAMPLE_NAME_EDIT=$(echo $SAMPLE_LINE | awk '{print $3}')
SAMPLE_DIRECTORY=$(echo $SAMPLE_LINE | awk '{print $4}')

# STEP 3:
# Rename and move files into new folder

mkdir /data/SPAdes/Renamed

cd $SAMPLE_DIRECTORY/$SAMPLE_NAME
scp ${TAXONOMY}_$SAMPLE_NAME_EDIT.fas /data/SPAdes/Renamed

cd /data/SPAdes/Renamed

# STEP 4:

# Next, we need to rename the headers in the file. The following commands can be achieved by the terminal via the terminal.
# Firstly, we have to make sure there are no spaces in the headers, using grep
# grep "^>.* .*" your_file.fasta

# If there are no lines that come up, this means none of the headers contain spaces. If needed, replace spaces with underscores
# Check the first header of the fasta file with:
# head your_file.fasta

# You can also check for the abundance of underscores with
# grep "^>.*_.*" your_file.fasta
# Make sure you hit ^C after you've seen a few, otherwise they will go on for ages with large fasta files

# Next, we want to add some information before the header, so that buscophy can make sense of the file
# We want Genus, species and id number, separated with underscores
# For example, for Parasa euchlora, we want to add:
# >Parasa_euchlora_T1054_[existing header information]

# Note that headers can be identified with > symbol
# Also note that we must ensure we do NOT use additional underscores. E.g., if you are dealing with subspecies, use a - to separate the names
# >Parasa_euchlora-euchlora_T1054
# This is because of how buscophy reads the first underscores of the file. In addition, ensure there are no underscores in the id numbers
# In my case, they were originally T1_054, but we already removed the underscore in the fasta file name, and we must make sure to not include it in the header
# If you really need, you can make an additional file with the translations e.g.
# Parasa euchlora euchlora T1_054 = Parasa_euchlora_euchlora_T1054

# We use sed to change the header name of the fasta files
# We want to identify the header with >, add information to the start, ensure the remaining header information stays there, and specify the file to deal with
# It is HIGHLY suggested that you have a backup of your scaffolds before doing this, just in case something goes wrong (it could be difficult to reverse!)
# To do this, we run command:

sed -i 's/^>/>Parasa_euchlora_T1054_/' Parasa_euchlora_T1054.fas

# Verifying the change with commands like:
# head your_file.fasta
# grep "^>.*_.*" your_file.fasta

# There are other programs like awk that can do similar. It shouldn't take long to do it all.
