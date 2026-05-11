library(GenSA)    # GenSA is better than optimx (although somewhat slower)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

# Create and appropriate folder and change into this for the analysis

setwd("/[Local_computer_pathway]/data/BioGeoBears")
getwd()

# Specify the tree file here. Remember, we want a tree with NO CONFIDENCE INTERVALS on the dates for this analysis. We made this when running LSD2
# It was saved as so earlier: /data/LSD2/Calibration/LSD2_no_CI.timetree.date.nexus

trfn = "LSD2_v2_no_CI.timetree.date.nexus"

# LSD2 (and perhaps other programs) creates trees in .nexus formats. We need to convert it to newick for BioGeoBears
# You can do this in programs like FigTree, or using the ape package:

library(ape)
tree <- read.nexus("LSD2_v2_no_CI.timetree.date.nexus")
trfn <- write.tree(tree, "LSD2_v2_no_CI.timetree.date.newick")

# Look at the raw Newick file:
moref(trfn)

# Look at your phylogeny (plots to a PDF, which avoids issues with multiple graphics in same window):
pdffn = "tree.pdf"
pdf(file=pdffn, width=9, height=12)

tr = read.tree(trfn)
tr
plot(tr)
title("Dated phylogeny of Parasa complex")
axisPhylo() # plots timescale

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)

# Create a geography file based on your taxon knowledge.

# Load geography file
geogfn = "Parasa_geog.data"
moref(geogfn)

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Maximum range size observed:
max(rowSums(dfnums_to_numeric(tipranges@df)))

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = 2

# Run DEC
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = "LSD2_v2_no_CI.timetree.date.newick"

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1

BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Enable time-stratified analysis by uncommenting and setting these files:
BioGeoBEARS_run_object$timesfn = "times.txt"
#BioGeoBEARS_run_object$distsfn = "distsfn.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "dispersal_multipliers_fn.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency_fn.txt"

# Optional files (leave commented out for now):
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Set up DEC model
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)


runslow = TRUE
resfn = "Parasa_DEC_M0_constrained_v1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}

runslow = FALSE
resfn = "Parasa_DEC_M0_constrained_v1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}

# Run DEC+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = "LSD2_v2_no_CI.timetree.date.RENAMED.newick"
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# Enable time-stratified analysis by uncommenting and setting these files:
BioGeoBEARS_run_object$timesfn = "times.txt"
#BioGeoBEARS_run_object$distsfn = "distsfn.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "dispersal_multipliers_fn.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency_fn.txt"

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Optional files (leave commented out for now):
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Parasa_DEC+J_M0_constrained_v1.Rdata"
runslow = TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
}



#######################################################
# PDF plots
#######################################################
pdffn = "Parasa_DEC_vs_DEC+J_M0_constrained_v1.pdf"
pdf(pdffn, width=10, height=12)

# --- pretty labels for plotting only (underscores -> spaces) ---
make_pretty <- function(x) gsub("_", " ", x, fixed = TRUE)

# 1) Pretty tree
tr_plot <- tr
tr_plot$tip.label <- make_pretty(tr_plot$tip.label)

# 2) Pretty tipranges (row names must match the tree labels)
tipranges_plot <- tipranges
rownames(tipranges_plot@df) <- make_pretty(rownames(tipranges_plot@df))

# (optional) safety: ensure 1:1 after prettifying
if (!setequal(tr_plot$tip.label, rownames(tipranges_plot@df))) {
  stop("Pretty labels don't match between tree and tipranges.")
}

base_colors <- c(
  "#6FB2E4",  # A
  "#5289B1",  # B
  "#8AB6A3",  # C
  "#469C76",  # D
  "#336951",  # E
  "#B7CEC2",  # F
  "#E5BB6B",  # G
  "#DCA237"   # H
)

library(grDevices)
blend_colors <- function(colors) {
  # Convert hex colors to RGB matrix (0-255)
  rgb_matrix <- t(col2rgb(colors))
  # Average RGB values
  avg_rgb <- colMeans(rgb_matrix)
  # Convert back to hex
  rgb(avg_rgb[1], avg_rgb[2], avg_rgb[3], maxColorValue=255)
}

possible_ranges_list_txt <- results_object$inputs$all_geog_states_list_usually_inferred_from_areas_maxareas
length(possible_ranges_list_txt)

colors_list_for_states <- sapply(possible_ranges_list_txt, function(state_indices) {
  # Check if state_indices is NA or contains NA
  if (all(is.na(state_indices))) {
    return("#FFFFFF")  # color for null/NA
  }
  # Map numeric indices to base colors
  selected_colors <- base_colors[state_indices + 1]  # +1 for R indexing
  if (length(selected_colors) == 1) {
    return(selected_colors)
  } else {
    # Blend colors for combinations
    blend_colors(selected_colors)
  }
})

# To check colour assignmnets:
# barplot(rep(1, length(colors_list_for_states)), col=colors_list_for_states, border=NA)



#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on Parasa M0_constrained"

max_age <- max(branching.times(tr))

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr_plot, tipranges=tipranges_plot, colors_list_for_states=colors_list_for_states)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.8, statecex=0.5, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr_plot, tipranges=tipranges_plot, colors_list_for_states=colors_list_for_states)

#######################################################
# Plot ancestral states - DECJ
#######################################################
analysis_titletxt ="BioGeoBEARS DEC+J on Parasa M0_constrained"

# Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr_plot, tipranges=tipranges_plot, colors_list_for_states=colors_list_for_states)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.8, statecex=0.5, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr_plot, tipranges=tipranges_plot, colors_list_for_states=colors_list_for_states)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it





#######################################################
#######################################################
# DIVALIKE AND DIVALIKE+J ANALYSIS
#######################################################
#######################################################
# NOTE: The BioGeoBEARS "DIVALIKE" model is not identical with 
# Ronquist (1997)'s parsimony DIVA. It is a likelihood
# interpretation of DIVA, constructed by modelling DIVA's
# processes the way DEC does, but only allowing the 
# processes DIVA allows (widespread vicariance: yes; subset
# sympatry: no; see Ronquist & Sanmartin 2011, Figure 4).
#
# DIVALIKE is a likelihood interpretation of parsimony
# DIVA, and it is "like DIVA" -- similar to, but not
# identical to, parsimony DIVA.
#
# I thus now call the model "DIVALIKE", and you should also. ;-)
#######################################################
#######################################################

#######################################################
# Run DIVALIKE
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = "LSD2_v2_no_CI.timetree.date.RENAMED.newick"
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "times.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "dispersal_multipliers_fn.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DIVALIKE model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "Parasa_DIVALIKE_M1_maxareas2_constained_v1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDIVALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKE = res
}

#######################################################
# Run DIVALIKE+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = "LSD2_v2_no_CI.timetree.date.RENAMED.newick"
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "times.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "dispersal_multipliers_fn.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DIVALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Parasa_DIVALIKE+J_M1_maxareas2_v1_constrained.Rdata"
runslow = TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDIVALIKEj = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKEj = res
}

pdffn = "Parasa_DIVALIKE_vs_DIVALIKE+J_M1_maxareas2_constrained_v1.pdf"
pdf(pdffn, width=10, height=12)

#######################################################
# Plot ancestral states - DIVALIKE
#######################################################
analysis_titletxt ="BioGeoBEARS DIVALIKE on Parasa M1_maxareas2 constrained"

# Setup
results_object = resDIVALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr_plot, tipranges=tipranges_plot, colors_list_for_states=colors_list_for_states)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.8, statecex=0.5, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr_plot, tipranges=tipranges_plot, colors_list_for_states=colors_list_for_states)

#######################################################
# Plot ancestral states - DIVALIKE+J
#######################################################
analysis_titletxt ="BioGeoBEARS DIVALIKE+J on Parasa M1_maxareas2 constrained"

# Setup
results_object = resDIVALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr_plot, tipranges=tipranges_plot, colors_list_for_states=colors_list_for_states)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.8, statecex=0.5, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr_plot, tipranges=tipranges_plot, colors_list_for_states=colors_list_for_states)

if (dev.cur() > 1) dev.off()
#dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)




#######################################################
#######################################################
# BAYAREALIKE AND BAYAREALIKE+J ANALYSIS
#######################################################
#######################################################
# NOTE: As with DIVA, the BioGeoBEARS BayArea-like model is 
# not identical with the full Bayesian model implemented 
# in the "BayArea" program of Landis et al. (2013). 
#
# Instead, this is a simplified likelihood interpretation
# of the model.  Basically, in BayArea and BioGeoBEARS-BAYAREALIKE, 
# "d" and "e" work like they do in the DEC model of Lagrange 
# (and BioGeoBEARS), and then BayArea's cladogenesis assumption
# (which is that nothing in particular happens at cladogenesis) is 
# replicated by BioGeoBEARS.
#
# This leaves out 3 important things that are in BayArea:
# 1. Distance dependence (you can add this with a distances 
#    matrix + the "x" parameter in BioGeoBEARS, however)
# 2. A correction for disallowing "e" events that drive
#    a species extinct (a null geographic range)
# 3. The neat Bayesian sampling of histories, which allows
#    analyses on large numbers of areas.
#
# The main purpose of having a "BAYAREALIKE" model is 
# to test the importance of the cladogenesis model on 
# particular datasets. Does it help or hurt the data 
# likelihood if there is no special cladogenesis process?
# 
# BAYAREALIKE is a likelihood interpretation of BayArea,
# and it is "like BayArea" -- similar to, but not
# identical to, Bayesian BayArea.
# I thus now call the model "BAYAREALIKE", and you should also. ;-)
#######################################################
#######################################################

#######################################################
# Run BAYAREALIKE
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = "LSD2_v2_no_CI.timetree.date.RENAMED.newick"
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "times.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "dispersal_multipliers_fn.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Check the inputs
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "Parasa_BAYAREALIKE_M1_maxareas2_constrained_v1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resBAYAREALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKE = res
}

#######################################################
# Run BAYAREALIKE+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = "LSD2_v2_no_CI.timetree.date.RENAMED.newick"
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "times.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "dispersal_multipliers_fn.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
# machines. I can't replicate this on my Mac machines, but it is almost certainly
# just some precision under-run issue, when optim/optimx tries some parameter value 
# just below zero.  The "min" and "max" options on each parameter are supposed to
# prevent this, but apparently optim/optimx sometimes go slightly beyond 
# these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
# slightly for each parameter:
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Parasa_BAYAREALIKE+J_M1_maxareas2_constrained_v1.Rdata"
runslow = TRUE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resBAYAREALIKEj = res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKEj = res
}



pdffn = "Parasa_BAYAREALIKE_vs_BAYAREALIKE+J_M1_maxareas2_constrained_v1.pdf"
pdf(pdffn, width=10, height=12)

base_colors <- c(
  "#6FB2E4",  # A
  "#5289B1",  # B
  "#8AB6A3",  # C
  "#469C76",  # D
  "#336951",  # E
  "#B7CEC2",  # F
  "#E5BB6B",  # G
  "#DCA237"   # H
)

library(grDevices)
blend_colors <- function(colors) {
  # Convert hex colors to RGB matrix (0-255)
  rgb_matrix <- t(col2rgb(colors))
  # Average RGB values
  avg_rgb <- colMeans(rgb_matrix)
  # Convert back to hex
  rgb(avg_rgb[1], avg_rgb[2], avg_rgb[3], maxColorValue=255)
}

# Setup
results_object = resBAYAREALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

possible_ranges_list_txt <- results_object$inputs$all_geog_states_list_usually_inferred_from_areas_maxareas
length(possible_ranges_list_txt)

colors_list_for_states <- sapply(possible_ranges_list_txt, function(state_indices) {
  # Check if state_indices is NA or contains NA
  if (all(is.na(state_indices))) {
    return("#FFFFFF")  # color for null/NA
  }
  # Map numeric indices to base colors
  selected_colors <- base_colors[state_indices + 1]  # +1 for R indexing
  if (length(selected_colors) == 1) {
    return(selected_colors)
  } else {
    # Blend colors for combinations
    blend_colors(selected_colors)
  }
})

# To check colour assignmnets:
# barplot(rep(1, length(colors_list_for_states)), col=colors_list_for_states, border=NA)

#######################################################
# Plot ancestral states - BAYAREALIKE
#######################################################
analysis_titletxt ="BioGeoBEARS BAYAREALIKE on Parasa M1_maxareas2 constrained"

# Setup
results_object = resBAYAREALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr_plot, tipranges=tipranges_plot, colors_list_for_states=colors_list_for_states)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.8, statecex=0.5, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr_plot, tipranges=tipranges_plot, colors_list_for_states=colors_list_for_states)

#######################################################
# Plot ancestral states - BAYAREALIKE+J
#######################################################
analysis_titletxt ="BioGeoBEARS BAYAREALIKE+J on Parasa M1_maxareas2 constrained"

# Setup
results_object = resBAYAREALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr_plot, tipranges=tipranges_plot, colors_list_for_states=colors_list_for_states)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.8, statecex=0.5, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr_plot, tipranges=tipranges_plot, colors_list_for_states=colors_list_for_states)

if (dev.cur() > 1) dev.off()
#dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)



#########################################################################
#########################################################################
#########################################################################
#########################################################################
# 
# CALCULATE SUMMARY STATISTICS TO COMPARE
# DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
# 
#########################################################################
#########################################################################
#########################################################################
#########################################################################

#########################################################################
#########################################################################
# REQUIRED READING:
#
# Practical advice / notes / basic principles on statistical model 
#    comparison in general, and in BioGeoBEARS:
# http://phylo.wikidot.com/advice-on-statistical-model-comparison-in-biogeobears
#########################################################################
#########################################################################

# Set up empty tables to hold the statistical results
restable = NULL
teststable = NULL

#######################################################
# Statistics -- DEC vs. DEC+J
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DEC, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DEC+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

# The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
# confer the same likelihood on the data. See: Brian O'Meara's webpage:
# http://www.brianomeara.info/tutorials/aic
# ...for an intro to LRT, AIC, and AICc

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DIVALIKE vs. DIVALIKE+J
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DIVALIKE, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DIVALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- BAYAREALIKE vs. BAYAREALIKE+J
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# BAYAREALIKE, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# BAYAREALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#########################################################################
# ASSEMBLE RESULTS TABLES: DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
#########################################################################
teststable$alt = c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
teststable$null = c("DEC", "DIVALIKE", "BAYAREALIKE")
row.names(restable) = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
restable = put_jcol_after_ecol(restable)
restable

# Look at the results!!
restable
teststable

#######################################################
# Save the results tables for later -- check for e.g.
# convergence issues
#######################################################

# Loads to "restable"
# save(restable, file="restable_v1.Rdata")
# load(file="restable_v1.Rdata")

# Loads to "teststable"
# save(teststable, file="teststable_v1.Rdata")
# load(file="teststable_v1.Rdata")

# Also save to text files
# write.table(restable, file="restable_constrained.txt", quote=FALSE, sep="\t")
# write.table(unlist_df(teststable), file="teststable_constrained.txt", quote=FALSE, sep="\t")

#######################################################
# Model weights of all six models
#######################################################
restable2 = restable

# With AICs:
AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable = cbind(restable, AICtable)
restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
conditional_format_table(restable_AIC_rellike)

# With AICcs -- factors in sample size
samplesize = length(tr$tip.label)
AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 = cbind(restable2, AICtable)
restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
restable_AICc_rellike = put_jcol_after_ecol(restable_AICc_rellike)
conditional_format_table(restable_AICc_rellike)


# Also save to text files
write.table(restable_AIC_rellike, file="restable_AIC_rellike_constrained.txt", quote=FALSE, sep="\t")
write.table(restable_AICc_rellike, file="restable_AICc_rellike_constrained.txt", quote=FALSE, sep="\t")
# 
# Save with nice conditional formatting
write.table(conditional_format_table(restable_AIC_rellike), file="restable_AIC_rellike_formatted_constrained.txt", quote=FALSE, sep="\t")
write.table(conditional_format_table(restable_AICc_rellike), file="restable_AICc_rellike_formatted_constrained.txt", quote=FALSE, sep="\t")


#######
# Just raw AIC values
# Create a list of your model results
models <- list(
  "DEC" = resDEC,
  "DEC+J" = resDECj,
  "DIVALIKE" = resDIVALIKE,
  "DIVALIKE+J" = resDIVALIKEj,
  "BAYAREALIKE" = resBAYAREALIKE,
  "BAYAREALIKE+J" = resBAYAREALIKEj
)

# Extract basic stats for each model
model_stats <- data.frame(
  Model = names(models),
  LogLik = sapply(models, get_LnL_from_BioGeoBEARS_results_object),
  Parameters = c(2, 3, 2, 3, 2, 3),  # Standard parameter counts
  stringsAsFactors = FALSE
)

# Calculate raw AIC values
model_stats$AIC <- -2 * model_stats$LogLik + 2 * model_stats$Parameters

# Calculate raw AICc values
samplesize <- length(tr$tip.label)
model_stats$AICc <- model_stats$AIC + (2 * model_stats$Parameters * (model_stats$Parameters + 1)) / (samplesize - model_stats$Parameters - 1)

# Display the table
print(model_stats)

# Save to file
write.table(model_stats, file="raw_AIC_table_constrained.txt", quote=FALSE, sep="\t", row.names=FALSE)

