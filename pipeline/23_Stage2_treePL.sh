#!/bin/bash
#SBATCH --job-name=treePL_step2
#SBATCH --output=step2_%A_%a.out
#SBATCH --error=step2_%A_%a.err
#SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

# Load modules
ml NLopt/2.4.2-GCCcore-7.3.0
ml GCCcore/9.3.0
ml ADOL-C/2.7.2-gompi-2020a

# STAGE 2 is the cross-validation analysis. Reconstruct the config file like so, inputting the majority opt, optad and optcvad values as obtained from STAGE 1.
# The [Priming command] should now be commented out, and the [Best optimisation parameters] and [Cross-validation analysis] should be commented in

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
#	#[Priming command]
#	#prime
#	[Best optimisation parameters]
#	opt = 1
#	moredetail
#	optad = 2
#	moredetailad
#	optcvad = 2
#	[Cross-validation analysis]
#	randomcv
#	cviter = 10
#	cvsimaniter = 100000
#	cvstart = 10
#	cvstop = 0.0001
#	cvmultstep = 0.5
#	cvoutfile = /data/treePL/Stage2/CV.txt
#	#[Best smoothing value]
#	#[Output file of dating step]

CONFIG_TEMPLATE=/data/treePL/config_STAGE2.txt
RUN_DIR=/data/treePL/Stage2
RUN_ID=$SLURM_ARRAY_TASK_ID

# Make unique config for each job
cp $CONFIG_TEMPLATE $RUN_DIR/run_stage2_${RUN_ID}.txt

# Replace cvoutfile with a unique one
sed -i "s|cvoutfile = .*|cvoutfile = $RUN_DIR/randomcv_Parasa_${RUN_ID}.txt|" $RUN_DIR/run_stage2_${RUN_ID}.txt

cd $RUN_DIR
treePL run_stage2_${RUN_ID}.txt

# After running 100 replicates, complete the following tasks in the terminal:

#	DIR=/data/treePL/Stage2

# Make a temporary file to store best lambdas
#	TMP_FILE="best_lambdas.txt"

# Loop over files 1–100
#	for i in $(seq 1 100); do
#	FILE="${DIR}/randomcv_Parasa_${i}.txt"
#	if [ ! -f "$FILE" ]; then
#	echo "Warning: $FILE not found, skipping."
#	continue
#	fi

# Use awk to find lowest chisq and its lambda
#  	awk '
#    	$1 == "chisq:" {
#  	# Remove parentheses, get lambda
#      	gsub(/\(|\)/, "", $2);
#     	 lambda = $2;
#     	 chisq = $3;
#     	 if (NR == 1 || chisq < min) {
#     	  min = chisq;
#     	   best = lambda;
#     	 }
#    	}
#    	END {
#    	  print best;
#   	 }
#  	' "$FILE" >> "$TMP_FILE"
#	done

#	echo "Extracted best lambdas for all files."

# Count occurrences and find the plurality winner
#	echo ""
#	echo "Smoothing values ranked by count:"
#	sort "$TMP_FILE" | uniq -c | sort -nr

# Make a note of the best smoothing value for STAGE 3





