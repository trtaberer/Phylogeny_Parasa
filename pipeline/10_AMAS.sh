#!/bin/bash
#SBATCH --job-name=AMAS_concat
#SBATCH --output=AMAS_concat.out
#SBATCH --error=AMAS_concat.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# STEP 1:

# After filtering with OliInSeq and trimAl, further filtering is conducted with Geneious (http://www.geneious.com/)
# At this stage, you can proceed only with the nt files
# Download the filtered .fasta files and process them with Geneious on your local computer
# In Geneious (no subscription needed), remove:
	# Alignments with less than 50% samples represented
	# Alignments with a GC content over 60%
	# Alignments with a pairwise identity over 98%
	# Alignments with a sequence length of less than 400 bp

# When done, re-upload these remaining alignments to the cluster for further processing
# For this, create a folder and upload the .fasta files there:

mkdir /data/filtered_alignments

# STEP 2:

# Next, we want to create summary statistics for each alignment with AMAS
# For this, we use the alignments remaining after OliInSeq, trimAl and Geneious

# Define the folder path where these alignment files have been saved, as well as an output folder:

nt_list=/data/filtered_alignments

Output=/data/AMAS
mkdir -p $Output

# STEP 3:

# Run AMAS to create summary statistics. The output will be a single .txt file in your defined output folder

AMAS.py summary -f fasta -d dna -i $nt_list/*.fasta -o $Output/AMAS_summary.txt

# STEP 4:

# Next, for running IQ-TREE, AMAS can also be used to concatenate the alignments into one large .fasta file
# Given this command, it will create both the partition file and concatenated alignment file required for IQ-TREE

AMAS.py concat \
    -f fasta -d dna -i $nt_list/*.fasta \
    -p $Output/partitions.txt -t $Output/concatenated.fasta \
    --part-format raxml