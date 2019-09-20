#' This script is set up to run the BEPA package with Cara's warfare project
#' 
#' This code and the functions in BEPA definitely need checking over!

require(phylopath)
require(dplyr)
require(plyr)
require(ape)
source("BEPA.R")
#require(foreach)
#require(doParallel)
#require(parallel)

#n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
#registerDoParallel(cores = (n_cores-2))

data <- read.csv("data_for_warfare_path.csv", stringsAsFactors = F)
colnames(data)

# Need to define the column that identifes the tips on the tree to use
taxa_names <- "Tree_Taxa_Name"

# Defining the row names in the data.frame as the taxa_names
row.names(data) <- data[,taxa_names]

# Specifying the variables to be used in path models
variable_list <- colnames(data)[-which(colnames(data) == taxa_names)]

# Check this looks ok  
variable_list

# This is unncessary when all columns are included in eiher variable_list, but tidies things up when more columns are included
data_use <- data[,c(taxa_names, variable_list)]
data_use <- na.omit(data_use)

#' ----------------------------------------------
# Reading in tree
tree <- read.nexus("tree_for_warfare_path.tree")

# Checking tree tip overlap against dataset
table(row.names(data) %in% tree$tip.label)
# Should all be true

#' Pruning trees
tree_taxa_to_drop <- tree$tip.label[!tree$tip.label %in% as.character(data[,taxa_names])]

# Pruning Trees
#pruned_trees <-  lapply(trees, function(trees) drop.tip(trees, tree_taxa_to_drop))
pruned_tree <- drop.tip(tree, tree_taxa_to_drop)

# ----------------------------------------------------
#' Getting all possible path combinations
#' Note: This is slow, and becomes exponentially slower the more variables are included.

#' Relationships specified in the optional exlusions arguement of get_all_models will not be included in the output.
#' Currently this is only set up to work with bivariate relationships.

#exclusions <- c(NPP_cont ~ qtstor,
#                NPP_cont ~ leader,
#                NPP_cont ~ RS,
#                NPP_cont ~ pathogen_scaled,
#                NPP_pred ~ qtstor,
#                NPP_pred ~ leader,
#                NPP_pred ~ RS,
#                NPP_pred ~ pathogen_scaled,
#                pathogen_scaled ~ qtstor,
#                pathogen_scaled ~ leader,
#                pathogen_scaled ~ RS)

all_models_raw <- get_all_models(variable_list = variable_list)
#all_models_raw <- get_all_models(variable_list = variable_list, exclusions = exclusions)

#' ----------------------------------------------
#' Converting all model from all_models_raw to matricies
all_models <- strings_to_model_sets(model_strings = all_models_raw)

#' ----------------------------------------------
#' Comparing all models

result <- phylo_path(all_models, 
                     data = data_use, 
                     tree = pruned_tree, 
                     model = 'lambda',
                     method = "logistic_MPLE")

# Summary of results
result_summary <- summary(result)

# Plotting the results summary
plot(result_summary[1:10,])


#' -----------------------------------
#' Plotting models
#' 
#' Plotting best model
best_model <- best(result)

pdf(file = "Best Model.pdf", width = 10, height = 6, compress = F, useDingbats = F)

  plot(best_model)
  
dev.off()

# Average across best models
average_model_full <- average(result, avg_method = "full")

pdf(file = "Average of Models.pdf", width = 10, height = 6, compress = F, useDingbats = F)

  plot(average_model_full)
  
dev.off()

#' More stuff here:
#' https://cran.r-project.org/web/packages/phylopath/vignettes/intro_to_phylopath.html


