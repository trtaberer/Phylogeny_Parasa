#!/bin/bash
#SBATCH --job-name=treePL_step1
#SBATCH --output=step1_%A_%a.out
#SBATCH --error=step1_%A_%a.err
#SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=28G

# Load modules
ml NLopt/2.4.2-GCCcore-7.3.0
ml GCCcore/9.3.0
ml ADOL-C/2.7.2-gompi-2020a

# Save the rooted_and_ladderized_bs.ufboot file as generated in R to folder /data/treePL as needed for later

# Define filenames and directories
# Tip: There are 3 distinct stages, so it is recommended to make output folders for this:

mkdir /data/treePL/Stage1
mkdir /data/treePL/Stage2
mkdir /data/treePL/Stage3

# For STAGE 1 and 2 of treePL, the rooted.nwk tree as generated for LSD2 is our input
# For STAGE 3, the bootstrap replicate tree is used

# STAGE 1

# Create a config.txt file. An example is given in Fig. 1 of https://doi.org/10.48550/arXiv.2008.07054
# There are the following included sections:
	
	#treefile=
	#[General commands]
	#[Calibrations]
	#[Priming command]
	#[Best optimisation parameters]
	#[Cross-validation analysis]
	#[Best smoothing value]
	#[Output file of dating step]

# Following https://doi.org/10.48550/arXiv.2008.07054, STAGE 1 is to prime the run from your ML tree.
# The config files specifies the rooted.nwk tree as the input, with the [General commands], [Calibrations] and [Priming command] activated
# The other sections are # out
# The config file for this stage looks like this (ignoring the first # on each line in your actual file):

#	treefile=data/LSD2/rooted.nwk
#	[General commands]
#	numsites = 3229490
#	nthreads = 8
#	thorough
#	log_pen
#	[Calibrations]
#	mrca = ROOT Taeda_prasina_T1075 Parasa_cuernavaca_T3013
#	min = ROOT 81.0
#	max = ROOT 90.0
#	mrca = CHRYSOPOLOMIDAE Strigivenifera_venata_T1087 Dalcera_abrasa_T5058
#	min = CHRYSOPOLOMIDAE 55.7
#	max = CHRYSOPOLOMIDAE 72.5
#	[Priming command]
#	prime
#	#[Best optimisation parameters]
#	#[Cross-validation analysis]
#	#[Best smoothing value]
#	#[Output file of dating step]

CONFIG_TEMPLATE=/data/treePL/config_STAGE1.txt
RUN_DIR=/data/treePL/Stage1
RUN_ID=$SLURM_ARRAY_TASK_ID

# Make a copy of the config file with unique output names
cp $CONFIG_TEMPLATE $RUN_DIR/run_stage1_${RUN_ID}.txt

# Optionally edit the config if needed (e.g., sed to set random seed or output name)

# Run treePL
cd $RUN_DIR
treePL run_${RUN_ID}.txt

# After running 100 replicates, complete the following tasks in the terminal:

# To count the number of opt, optad and optcvad values to find the majority for STAGE 2, go to your .out files from running STAGE 1 and:

# Go to your output folder
# cd /data/treePL/Stage1

# Extract last 7 lines from all files and grep the line of interest
	# for f in step1_*.out; do tail -n 7 "$f"; done | grep '^opt =' | sort | uniq -c | sort -nr
	# for f in step1_*.out; do tail -n 7 "$f"; done | grep '^optad =' | sort | uniq -c | sort -nr
	# for f in step1_*.out; do tail -n 7 "$f"; done | grep '^optcvad =' | sort | uniq -c | sort -nr

# Make a note of the majority opt, optad and optcvad values for STAGE 2.