#!/bin/bash
#SBATCH --job-name=LSD2
#SBATCH --output=LSD2.out
#SBATCH --error=LSD2.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00

# We are now able to time calibrate our phylogeny!
# Based on the two methods used, this example will utilise the phylogeny as created from the IQ-TREE run on the concatenated alignment

# STEP 1:
# LSD2 requires a rooted phylogeny in Newick format as input. As IQ-TREE outputs a .treefile, you can download this and open it in a tree viewer program like FigTree
# Here, you can root the phylogeny on your outgroup and export it as a .nwk, ready for LSD2
# When done, upload to /data/LSD2

# STEP 2:
# Load LSD2

ml LSD2/2.4.1-GCCcore-12.2.0

# STEP 3:
# Define input and output directories:

INPUT=/data/LSD2
OUTPUT=/data/LSD2/Calibration

# STEP 4:
# Create a file outgroup.txt with your outgroup. It should look like the following if using one outgroup, with the sample name as used throughout:

	#1
	#Taeda_prasina_T1075

# STEP 5:
# Create a calibrations.txt file with calibrations point throughout the phylogeny. It should look like this, specifying the mrca of the taxa with divergence values in millions of years:

	#1
	#mrca(Strigivenifera_vernata_T1087,Dalcera_abrasa_T5058) b(-72.5,-55.7)

# If there is also a calibration point at the base, this can be specified with -a "b(-[from],-[to])" in the LSD2 command itself.

# STEP 6:
# Run LSD2
	# Without confidence intervals (e.g. for BioGeoBears):

lsd2 -i $INPUT/rooted.nwk -s 3229490 -o $OUTPUT/LSD2_no_CI.timetree -g $INPUT/outgroup.txt \
     -r k -d $INPUT/calibrations.txt -a "b(-90.0,-81.0)" -z 0 \
     -l 0 -u e

	#With confidence intervals:

lsd2 -i $INPUT/rooted.nwk -s 3229490 -o $OUTPUT/LSD2_with_CI.timetree -g $INPUT/outgroup.txt  \
     -r k -d $INPUT/calibrations.txt -a "b(-90.0,-81.0)" -z 0 -f 100 \
     -l 0 -u e

# Note that -s flag is the total alignment length used to generate the rooted.nwk input tree