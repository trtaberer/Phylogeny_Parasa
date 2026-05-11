library(ape)

# After generating our bootstrap replicates, the trees need to be rooted on the outgroup and ladderised in R

# STEP 1:
# Download files from /data/treePL and set this as your working directory

setwd("[Local_computer_pathway]/data/treePL")

# STEP 2:
# Read the bootstrap trees

trees <- read.tree("constrained_bs.ufboot")

# STEP 3:
# Re-root on the outgroup and ladderize each tree

processed_trees <- lapply(trees, function(tr) {
  if ("Taeda_prasina_T1075" %in% tr$tip.label) {
    tr_rooted <- root(tr, outgroup = "Taeda_prasina_T1075", resolve.root = TRUE)
    ladderize(tr_rooted, right = FALSE)  # decreasing node order (root to tips)
  } else {
    stop("Taeda_prasina_T1075 not found in a tree! Check tree tips.")
  }
})

# STEP 4:
# Write all processed trees to a new file
write.tree(processed_trees, file = "rooted_and_ladderized_bs.ufboot")

cat("All trees re-rooted on Taeda_prasina_T1075 and ladderized (decreasing node order).\n")
