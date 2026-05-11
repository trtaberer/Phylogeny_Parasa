#!/bin/bash
#SBATCH --job-name=Concord_vect
#SBATCH --output=Concord_vect.out
#SBATCH --error=Concord_vect.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=01:00:00

# We can then test for concordance between the phylogenies created with both the concatenated alignment and the consensus gene tree alignments
# This particular .sh file is designed to be done in stages; use #s to blank out other STEPs as you go

# STEP 1:
# First, we need to estimate the support and quartet concordance vectors in ASTRAL
# This doesn’t take long at all, so can be done on base cluster without submitting a job

# The q file is the .tree file from the concatenated run
# The i file is the file containing the list of all the gene trees

QFILE=/data/IQ-Tree/TreeSearch
IFILE=/data/Gene_trees
OUTPUT=/data/Concordance
mkdir -p $Output

# Load and run ASTRAL and Java

ASTRAL=/data/Software/ASTRAL-master/Astral
module load Java/11.0.2

java -Xmx128G -jar $ASTRAL/astral.5.7.8.jar -q $QFILE/TS.treefile -i $IFILE/all_gene_trees.tre -t 2 -o $OUTPUT/astral_annotated.tree 2> $OUTPUT/astral_annotated.log

cd /data/Concordance

# Ensure you are in the Concordance output file before running the next steps

# STEP 2:
# Next, we want concordance vectors
# This step will take a while!
# Load IQ-TREE

ml IQ-TREE/2.2.2.6-gompi-2022b

# Specify the model, output and alignment. 
# The model is the best_model.nex file from our IQ-TREE Tree Search (TS) step
# The alignment is the concatenated alignment as made by AMAS, which was then used for the IQ-TREE MF and TS steps

MODEL=/data/IQ-Tree/TreeSearch
OUTPUT=/data/Concordance
ALIGNMENT=/data/AMAS

iqtree2 -te $OUTPUT/astral_annotated.tree -s $ALIGNMENT/concatenated.fasta -p $MODEL/TS.best_model.nex --scfl 100 --prefix scfl -nt 32

# STEP 3:
# Calculate the gene concordance vectors

INPUT1=/data/IQ-Tree/TreeSearch
INPUT2=/data/Gene_trees

iqtree2 -te $INPUT1/TS.treefile --gcf $INPUT2/all_gene_trees.tre --prefix gcf -T 32

# STEP 4:
# Finally we do a dummy analysis in IQ-TREE. The only point of this is to get the branch lengths in coalescent units 
# from the ASTRAL analysis, in a format that is output by IQ-TREE in a convenient table with IQ-TREE branch ID's.
# Note the -blfix option, which keeps the original branch lengths - this makes the scfs meaningless, but is here 
# simply to allow us to extract branch lengths in coalescent units from the ASTRAL tree in a convenient table
# We set scfl to 1, which saves time given the scfs are already meaningless - never use the sCFs from this analysis!

INPUT1=/data/Concordance
INPUT2=/data/IQ-Tree/TreeSearch
ALIGNMENT=/data/AMAS

iqtree2 -te $INPUT1/astral_annotated.tree -s $ALIGNMENT/concatenated.fasta -blfix -p $INPUT2/TS.best_model.nex --scfl 1 --prefix coalescent_bl -T 32
