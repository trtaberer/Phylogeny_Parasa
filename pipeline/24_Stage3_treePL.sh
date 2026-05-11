#!/bin/bash
#SBATCH --job-name=treePL_step3
#SBATCH --output=step3_%A_%a.out
#SBATCH --error=step3_%A_%a.err
#SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

# Load modules
ml NLopt/2.4.2-GCCcore-7.3.0
ml GCCcore/9.3.0
ml ADOL-C/2.7.2-gompi-2020a

# In STAGE 3, the bootstrap replicates are dated using the optimisation parameters and best smoothing value
# Here, the bootstrap replicate tree is used in the config file as opposed to the rooted.nwk tree
# The [Cross-validation analysis] is now commented out, but the [Best smoothing value] and [Output file of dating step] are in.
# The config file should look as follows:

#	treefile=data/treePL/rooted_and_ladderized_bs.ufboot
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
#	#[Cross-validation analysis]
#	#randomcv
#	#cviter = 10
#	#cvsimaniter = 100000
#	#cvstart = 10
#	#cvstop = 0.0001
#	#cvmultstep = 0.5
#	#cvoutfile = /data/treePL/Stage2/CV.txt
#	[Best smoothing value]
#	100
#	[Output file of dating step]
#	/data/treePL/Stage3/treePL_dated.tre

CONFIG_TEMPLATE=/data/treePL/config_STAGE3.txt
RUN_DIR=/data/treePL/config_STAGE3
RUN_ID=$SLURM_ARRAY_TASK_ID

cp $CONFIG_TEMPLATE $RUN_DIR/run_stage3_${RUN_ID}.txt

# Make outfile unique
sed -i "s|outfile = .*|outfile = $RUN_DIR/treePL_Parasa_dated_${RUN_ID}.tre|" $RUN_DIR/run_stage3_${RUN_ID}.txt

cd $RUN_DIR
treePL run_stage3_${RUN_ID}.txt

# After running 100 replicates, complete the following tasks in the terminal:

# 	TREE_DIR=/data/treePL/config_STAGE3
# 	OUTPUT="${TREE_DIR}/treePL_Parasa_dated_ALL.tre"

#	rm -f "$OUTPUT"

#	for i in $(seq 1 100); do
#	    FILE="${TREE_DIR}/treePL_Parasa_dated_${i}.tre"
#	    if [[ -f "$FILE" ]]; then
#	        cat "$FILE" >> "$OUTPUT"
#	        # Do NOT add echo "" here — the tree file already ends with a newline!
#	    else
#	        echo "Warning: $FILE not found!"
#	    fi
#	done

#	echo "Combined trees written to: $OUTPUT"
