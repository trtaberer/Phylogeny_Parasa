library(GenSA)    # GenSA is better than optimx (although somewhat slower)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)
library(ape)
library(qgraph)

setwd("/[Local_computer_pathway]/data/BioGeoBears/BSM")
getwd()

# Specify the tree file here. Remember, we want a tree with NO CONFIDENCE INTERVALS on the dates for this analysis!
# It was saved as so earlier: /data/LSD2/Calibration/LSD2_no_CI.timetree.date.nexus

trfn = "LSD2_v2_no_CI.timetree.date.nexus"

# LSD2 (and perhaps other programs) creates trees in .nexus formats. We need to convert it to newick for BioGeoBears
# You can do this in programs like FigTree, or using the ape package:

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
BioGeoBEARS_run_object$trfn = "LSD2_v2_no_CI.timetree.date.RENAMED.newick"

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
# Time-stratified Biogeographic Stochastic Mapping (BSM)
#######################################################
model_name = "DECJ_M3_timestrat"
res = resDECj

pdffn = paste0("Parasa_", model_name, "_v1.pdf")
pdf(pdffn, height=6, width=6)

analysis_titletxt = paste0(model_name, " on Parasa")

# Setup
results_object = res
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#######################################################
# Stochastic mapping on DEC M3b stratified with islands coming up
#######################################################
clado_events_tables = NULL
ana_events_tables = NULL
lnum = 0

#######################################################
# Get the inputs for Biogeographical Stochastic Mapping
# Note: this can be slow for large state spaces and trees, since 
# the independent likelihoods for each branch are being pre-calculated
# E.g., for 10 areas, this requires calculation of a 1024x1024 matrix
# for each branch.  On a tree with ~800 tips and thus ~1600 branches, this was about 1.6 gigs
# for storage of "BSM_inputs_file.Rdata".
# Update: 2015-09-23 -- now, if you used multicore functionality for the ML analysis,
# the same settings will be used for get_inputs_for_stochastic_mapping().
#######################################################
BSM_inputs_fn = "BSM_inputs_file.Rdata"
runInputsSlow = TRUE
if (runInputsSlow)
{
  # debug:
  # cluster_already_open=FALSE; rootedge=FALSE; statenum_bottom_root_branch_1based=NULL; printlevel=1; min_branchlength=0.000001
  stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
  save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
} else {
  # Loads to "stochastic_mapping_inputs_list"
  load(BSM_inputs_fn)
} # END if (runInputsSlow)

# Check inputs (doesn't work the same on unconstr)
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
set.seed(seed=as.numeric(Sys.time()))

runBSMslow = TRUE
if (runBSMslow == TRUE)
{
  # Saves to: RES_clado_events_tables.Rdata
  # Saves to: RES_ana_events_tables.Rdata
  # Bug check:
  # stochastic_mapping_inputs_list=stochastic_mapping_inputs_list; maxnum_maps_to_try=100; nummaps_goal=50; maxtries_per_branch=40000; save_after_every_try=TRUE; savedir=getwd(); seedval=12345; wait_before_save=0.01; master_nodenum_toPrint=0
  
  BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=1000, nummaps_goal=1000, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01, master_nodenum_toPrint=0)
  
  RES_clado_events_tables = BSM_output$RES_clado_events_tables
  RES_ana_events_tables = BSM_output$RES_ana_events_tables
} else {
  # Load previously saved...
  
  # Loads to: RES_clado_events_tables
  load(file="RES_clado_events_tables.Rdata")
  # Loads to: RES_ana_events_tables
  load(file="RES_ana_events_tables.Rdata")
  BSM_output = NULL
  BSM_output$RES_clado_events_tables = RES_clado_events_tables
  BSM_output$RES_ana_events_tables = RES_ana_events_tables
} # END if (runBSMslow == TRUE)

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables
head(clado_events_tables[[1]])
head(ana_events_tables[[1]])
length(clado_events_tables)
length(ana_events_tables)

include_null_range = TRUE
areanames = names(tipranges@df)
areas = areanames
max_range_size = 2

# Note: If you did something to change the states_list from the default given the number of areas, you would
# have to manually make that change here as well! (e.g., areas_allowed matrix, or manual reduction of the states_list)
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)

colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)




############################################
# Setup for painting a single stochastic map
############################################
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = TRUE
clado_events_table = clado_events_tables[[1]]
ana_events_table = ana_events_tables[[1]]

# cols_to_get = names(clado_events_table[,-ncol(clado_events_table)])
# colnums = match(cols_to_get, names(ana_events_table))
# ana_events_table_cols_to_add = ana_events_table[,colnums]
# anagenetic_events_txt_below_node = rep("none", nrow(ana_events_table_cols_to_add))
# ana_events_table_cols_to_add = cbind(ana_events_table_cols_to_add, anagenetic_events_txt_below_node)
# rows_to_get_TF = ana_events_table_cols_to_add$node <= length(tr$tip.label)
# master_table_cladogenetic_events = rbind(ana_events_table_cols_to_add[rows_to_get_TF,], clado_events_table)

############################################
# Open a PDF
############################################
pdffn = paste0(model_name, "_single_stochastic_map_n1.pdf")
pdf(file=pdffn, height=6, width=6)

# Convert the BSM into a modified res object
master_table_cladogenetic_events = clado_events_tables[[1]]
resmod = stochastic_map_states_into_res(res=res, master_table_cladogenetic_events=master_table_cladogenetic_events, stratified=stratified)

plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)

# Paint on the branch states
paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=5, lty=par("lty"), root.edge=TRUE, stratified=stratified)

plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)

############################################
# Close PDF
############################################
dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
# Plot all 50 stochastic maps to PDF
#######################################################
# Setup
include_null_range = include_null_range
areanames = areanames
areas = areanames
max_range_size = max_range_size
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = stratified

# Loop through the maps and plot to PDF
pdffn = paste0(model_name, "_", length(clado_events_tables), "BSMs_v1.pdf")
pdf(file=pdffn, height=6, width=6)

nummaps_goal = 1000
for (i in 1:nummaps_goal)
{
  clado_events_table = clado_events_tables[[i]]
  analysis_titletxt = paste0(model_name, " - Stochastic Map #", i, "/", nummaps_goal)
  plot_BSM(results_object=res, clado_events_table=clado_events_table, stratified=stratified, analysis_titletxt=analysis_titletxt, addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, show.tip.label=TRUE, include_null_range=include_null_range)
} # END for (i in 1:nummaps_goal)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
# Summarize stochastic map tables
#######################################################
length(clado_events_tables)
length(ana_events_tables)

head(clado_events_tables[[1]][,-20])
tail(clado_events_tables[[1]][,-20])

head(ana_events_tables[[1]])
tail(ana_events_tables[[1]])

areanames = names(tipranges@df)
actual_names = areanames
actual_names

# Get the dmat and times (if any)
dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)
dmat_times

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables

# Simulate the source areas
BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res, clado_events_tables, ana_events_tables, areanames)
clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
ana_events_tables = BSMs_w_sourceAreas$ana_events_tables

# Count all anagenetic and cladogenetic events
counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)

summary_counts_BSMs = counts_list$summary_counts_BSMs
print(conditional_format_table(summary_counts_BSMs))

# Histogram of event counts
hist_event_counts(counts_list, pdffn=paste0(model_name, "_histograms_of_event_counts.pdf"))

#######################################################
# Print counts to files
#######################################################
tmpnames = names(counts_list)
cat("\n\nWriting tables* of counts to tab-delimited text files:\n(* = Tables have dimension=2 (rows and columns). Cubes (dimension 3) and lists (dimension 1) will not be printed to text files.) \n\n")
for (i in 1:length(tmpnames))
{
  cmdtxt = paste0("item = counts_list$", tmpnames[i])
  eval(parse(text=cmdtxt))
  
  # Skip cubes
  if (length(dim(item)) != 2)
  {
    next()
  }
  
  outfn = paste0(tmpnames[i], ".txt")
  if (length(item) == 0)
  {
    cat(outfn, " -- NOT written, *NO* events recorded of this type", sep="")
    cat("\n")
  } else {
    cat(outfn)
    cat("\n")
    write.table(conditional_format_table(item), file=outfn, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
  } # END if (length(item) == 0)
} # END for (i in 1:length(tmpnames))
cat("...done.\n")

#######################################################
# Check that ML ancestral state/range probabilities and
# the mean of the BSMs approximately line up
#######################################################
library(MultinomialCI)    # For 95% CIs on BSM counts
check_ML_vs_BSM(res, clado_events_tables, model_name, tr=NULL, plot_each_node=FALSE, linreg_plot=TRUE, MultinomialCI=TRUE)


install.packages("qgraph")

library(qgraph)

# (A) Overall network
pdf("qgraph_overall_dispersal.pdf", width=7, height=7)
qgraph(
  disp_overall_rel,
  layout = "circle",      # or "spring" if you want force-directed
  edge.labels = TRUE,
  minimum = 1,            # as specified
  cut = 4,
  labels = colnames(disp_overall_rel),
  directed = TRUE
)
dev.off()

# (B) Per time slice networks
if (!is.null(disp_by_slice_rel)) {
  dir.create("qgraph_time_slices", showWarnings = FALSE)
  for (nm in names(disp_by_slice_rel)) {
    M <- disp_by_slice_rel[[nm]]
    diag(M) <- 0
    fn <- file.path("qgraph_time_slices", paste0("qgraph_disp_", nm, ".pdf"))
    pdf(fn, width=7, height=7)
    qgraph(
      M,
      layout = "circle",
      edge.labels = TRUE,
      minimum = 0.25,     # as specified
      cut = 1,
      labels = colnames(M),
      directed = TRUE
    )
    dev.off()
  }
}






####summary stats
# === BEGIN: create single CSV summary from counts_list ===

# Defensive checks
if (!exists("counts_list")) stop("counts_list not found in environment. Run count_ana_clado_events() first.")
# Determine how many maps were run
n_maps <- max(length(if (exists("clado_events_tables")) clado_events_tables else NULL),
              length(if (exists("ana_events_tables")) ana_events_tables else NULL), na.rm = TRUE)
if (is.infinite(n_maps) || is.na(n_maps)) n_maps <- NA_integer_

# Helper to safely get an element from counts_list or NULL
get_if_exists <- function(name) {
  if (name %in% names(counts_list)) return(counts_list[[name]])
  return(NULL)
}

# Grab mean & sd matrices for "all dispersals"
mat_mean_all <- get_if_exists("all_dispersals_counts_fromto_means")
mat_sd_all   <- get_if_exists("all_dispersals_counts_fromto_sds")
cube_all     <- get_if_exists("all_dispersals_counts_cube")  # 3D: from x to x map

# Also anagenetic-only and founder-only (optional)
mat_mean_ana <- get_if_exists("ana_dispersals_counts_fromto_means")
mat_sd_ana   <- get_if_exists("ana_dispersals_counts_fromto_sds")
cube_ana     <- get_if_exists("anagenetic_dispersals_counts_cube")

mat_mean_founder <- get_if_exists("founder_counts_fromto_means")
mat_sd_founder   <- get_if_exists("founder_counts_fromto_sds")
cube_founder     <- get_if_exists("founder_counts_cube")

# Sanity: ensure mat_mean_all exists
if (is.null(mat_mean_all)) stop("counts_list$all_dispersals_counts_fromto_means not found. Cannot proceed.")

# Row/col names
from_names <- rownames(mat_mean_all)
to_names   <- colnames(mat_mean_all)
if (is.null(from_names) || is.null(to_names)) {
  # fallback: use numeric indices
  from_names <- as.character(seq_len(nrow(mat_mean_all)))
  to_names   <- as.character(seq_len(ncol(mat_mean_all)))
}

# Build data.frame of all directed pairs (exclude self->self)
pairs <- expand.grid(from = from_names, to = to_names, stringsAsFactors = FALSE)
pairs <- pairs[pairs$from != pairs$to, , drop = FALSE]

# Pull mean & sd for all
get_mat_val <- function(mat, r, c) {
  if (is.null(mat)) return(NA_real_)
  # ensure indexing by names works
  return(as.numeric(mat[r, c]))
}

pairs$mean_all <- mapply(function(r,c) get_mat_val(mat_mean_all, r, c), pairs$from, pairs$to)
pairs$sd_all   <- mapply(function(r,c) get_mat_val(mat_sd_all, r, c), pairs$from, pairs$to)

# Scaled version (as in paper)
pairs$mean_div100_all <- pairs$mean_all / 100

# maps_with_event_all and 95% quantiles from cube (if cube exists)
pairs$maps_with_event_all <- NA_integer_
pairs$q025_all <- NA_real_
pairs$q975_all <- NA_real_

if (!is.null(cube_all)) {
  dims <- dim(cube_all)  # expect c(n_from, n_to, n_maps)
  # Check dims match matrix dims
  # We'll assume cube ordering aligns with mat_mean_all row/colnames
  for (i in seq_len(nrow(pairs))) {
    rname <- pairs$from[i]
    cname <- pairs$to[i]
    # convert names to indices
    r_idx <- which(from_names == rname)
    c_idx <- which(to_names == cname)
    if (length(r_idx)==1 && length(c_idx)==1) {
      vec <- cube_all[r_idx, c_idx, ]
      # If cube stores totals across maps (rare), we'll still compute quantiles on its third dimension
      pairs$maps_with_event_all[i] <- sum(vec > 0, na.rm = TRUE)
      pairs$q025_all[i] <- as.numeric(quantile(vec, probs = 0.025, na.rm = TRUE))
      pairs$q975_all[i] <- as.numeric(quantile(vec, probs = 0.975, na.rm = TRUE))
    } else {
      pairs$maps_with_event_all[i] <- NA_integer_
      pairs$q025_all[i] <- NA_real_
      pairs$q975_all[i] <- NA_real_
    }
  }
} else {
  # If no cube, try infer maps_with_event from counts_list$all_dispersals_totals_list if present
  tot_list_all <- get_if_exists("all_dispersals_totals_list")
  if (!is.null(tot_list_all)) {
    # tot_list_all may be a list of matrices per map; count maps where >0
    for (i in seq_len(nrow(pairs))) {
      rname <- pairs$from[i]; cname <- pairs$to[i]
      cnt <- 0
      for (m in seq_along(tot_list_all)) {
        M <- tot_list_all[[m]]
        if (!is.null(M) && rname %in% rownames(M) && cname %in% colnames(M)) {
          v <- M[rname, cname]
          if (!is.na(v) && v > 0) cnt <- cnt + 1
        }
      }
      pairs$maps_with_event_all[i] <- cnt
    }
  }
}

# Add anagenetic columns if available
if (!is.null(mat_mean_ana)) {
  pairs$mean_ana <- mapply(function(r,c) get_mat_val(mat_mean_ana, r, c), pairs$from, pairs$to)
  pairs$sd_ana   <- mapply(function(r,c) get_mat_val(mat_sd_ana, r, c), pairs$from, pairs$to)
} else {
  pairs$mean_ana <- NA_real_; pairs$sd_ana <- NA_real_
}

# Add founder columns if available
if (!is.null(mat_mean_founder)) {
  pairs$mean_founder <- mapply(function(r,c) get_mat_val(mat_mean_founder, r, c), pairs$from, pairs$to)
  pairs$sd_founder   <- mapply(function(r,c) get_mat_val(mat_sd_founder, r, c), pairs$from, pairs$to)
} else {
  pairs$mean_founder <- NA_real_; pairs$sd_founder <- NA_real_
}

# Add n_maps column
pairs$n_maps <- ifelse(is.na(n_maps), NA_integer_, n_maps)

# Reorder columns for readability
out_df <- pairs[, c("from","to","n_maps",
                    "mean_all","sd_all","mean_div100_all","maps_with_event_all","q025_all","q975_all",
                    "mean_ana","sd_ana","mean_founder","sd_founder")]

# Sort by descending mean_all for convenience
out_df <- out_df[order(-out_df$mean_all, out_df$from, out_df$to), ]

# Write CSV
outfn <- "BSM_pairwise_summary_all_vs_ana_founder.csv"
write.csv(out_df, file = outfn, row.names = FALSE)
cat("Wrote CSV to", outfn, "\n")
# Print top 20 rows for quick check
print(head(out_df, 20))

# === END ===


# === Export raw pairwise counts with rows = simulations and columns = "AtoB" style ===

# REQUIREMENTS: counts_list$all_dispersals_counts_cube must exist in your environment.

if (!exists("counts_list")) stop("counts_list not found. Run count_ana_clado_events() first.")
if (!("all_dispersals_counts_cube" %in% names(counts_list))) stop("counts_list$all_dispersals_counts_cube not found.")

cube_all <- counts_list$all_dispersals_counts_cube
dims <- dim(cube_all)
if (length(dims) != 3) stop("Expected a 3D cube (from x to x map); got dims: ", paste(dims, collapse = "x"))
n_from <- dims[1]; n_to <- dims[2]; n_maps <- dims[3]

# Region names (from & to). If your cube already has names, we'll use them.
from_names <- rownames(cube_all)
to_names   <- colnames(cube_all)

# OPTIONAL: If you want to force A..H labels (only do if you know the order matches),
# uncomment and adapt the next lines:
if (n_from == 8 && n_to == 8) {
   forced_labels <- c("A","B","C","D","E","F","G","H")
   from_names <- forced_labels
   to_names   <- forced_labels
 }

# Fallback to generic names if missing
if (is.null(from_names)) from_names <- paste0("from", seq_len(n_from))
if (is.null(to_names))   to_names   <- paste0("to", seq_len(n_to))

# Build ordered list of directed pairs (exclude self->self)
pairs_idx <- expand.grid(i = seq_len(n_from), j = seq_len(n_to), stringsAsFactors = FALSE)
pairs_idx <- pairs_idx[pairs_idx$i != pairs_idx$j, , drop = FALSE]

# Order pairs in a reproducible way: by from (A..H) then to (A..H)
# We'll use the names given (which may be A..H or full region names). If you want a specific A..H ordering, set forced_labels above.
pairs_idx$from_name <- from_names[pairs_idx$i]
pairs_idx$to_name   <- to_names[pairs_idx$j]

# Create column names like "AtoB" but safe (no spaces)
make_colname <- function(f, t) {
  # Replace spaces or punctuation if present
  f_clean <- gsub("[^0-9A-Za-z]", "", f)
  t_clean <- gsub("[^0-9A-Za-z]", "", t)
  paste0(f_clean, "to", t_clean)
}
pairs_idx$colname <- mapply(make_colname, pairs_idx$from_name, pairs_idx$to_name, USE.NAMES = FALSE)

# Optionally reorder columns alphabetically by colname:
pairs_idx <- pairs_idx[order(pairs_idx$colname), ]

# Build output data.frame: rows = maps, cols = the pair columns
out_mat <- matrix(0L, nrow = n_maps, ncol = nrow(pairs_idx))
colnames(out_mat) <- pairs_idx$colname
rownames(out_mat) <- paste0("map_", seq_len(n_maps))

for (k in seq_len(n_maps)) {
  M_k <- cube_all[,,k]        # from x to matrix for map k
  # Fill out columns in the same order as pairs_idx
  # For each directed pair (i,j), extract M_k[i,j]
  vals <- numeric(nrow(pairs_idx))
  for (r in seq_len(nrow(pairs_idx))) {
    i <- pairs_idx$i[r]; j <- pairs_idx$j[r]
    vals[r] <- as.numeric(M_k[i, j])
  }
  out_mat[k, ] <- vals
}

# Convert to data.frame with "map" column first
out_df <- as.data.frame(out_mat, stringsAsFactors = FALSE, check.names = FALSE)
out_df <- cbind(map = paste0("map_", seq_len(n_maps)), out_df)

# Save CSV
outfn <- "BSM_raw_pairwise_counts_maps_as_rows.csv"
write.csv(out_df, file = outfn, row.names = FALSE, quote = FALSE)
cat("Wrote CSV:", outfn, "\n")
cat("Dimensions: rows (maps) =", nrow(out_df), "columns (pairs) =", ncol(out_df)-1, "\n")

# Quick preview
print("First 6 rows and first 12 columns:")
print(utils::head(out_df[, seq_len(min(ncol(out_df), 13))], 6))

# Load counts cube
cube <- counts_list$all_dispersals_counts_cube
n_maps <- dim(cube)[3]

# Compute mean number of events per map for each pair (i->j)
mean_matrix <- apply(cube, c(1,2), mean, na.rm = TRUE)

# Zero the diagonal
diag(mean_matrix) <- 0

# Compute total emigration (row sums) and immigration (column sums)
emigration_totals  <- rowSums(mean_matrix, na.rm = TRUE)
immigration_totals <- colSums(mean_matrix, na.rm = TRUE)

# Combine into a display table like the supplement
regions <- colnames(mean_matrix)
left_block <- round(mean_matrix, 2)
right_block <- data.frame(Emigration = round(emigration_totals, 2),
                          Immigration = round(immigration_totals, 2))

# Combine horizontally
combined_table <- cbind(left_block, "", right_block)
write.csv(combined_table, "BSM_summary_pairwise_with_emi_imm.csv", row.names = TRUE)
cat("Saved: BSM_summary_pairwise_with_emi_imm.csv\n")



