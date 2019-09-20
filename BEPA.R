#' Brute-force Exploratory Path Analysis (BEPA)
#' 
#' Package still in development, please do not circulate.
#' 
#' Joseph Watts - me@josephwatts.org
#' 
#' The BEPA package contains functions to help perform exploratory path analysis in R.
#' As the name suggests, this is not done efficiently or elegantly.
#' The package can be used to identify all possible path model structures, given a set of variables.
#' All of these models can then be compared. 
#' Time taken increases exponentially with the number of variables.
#' 
#' This package is written specifically for use with the phylopath package, which itself is based on the phyloglm and phylolm packages.
#' The functions should also be useable for exploratory path analyses more generally, but this has not been tested.
#' 
#' ------------------------------------------
#' Working notes:
#' 1. This code is generally ugly and inefficient. Still needs checking and tidying.
#' 
#' 2. The exclusion function could be made more general
#' Would probably be good to be able to exclude a wider range of models in the future.
#' For example, might want to exclude multiple predictors of a single variable, or require a variable to be predicted by others.
#' 
#' 3. Need to look into parallel processing
#' Currently there is a commented out piece of code that would parallelise the most time consuming 
#' part of the script, but it wont work on Windows (mclapply doesnt run in parallel on Windows.)
#' 
#' 4. Need to look into how the run this with multiple trees.
#' Maybe worth looking into creating a function to run multiple trees.
#' 
#' 5. The outputs of these function still need careful checking.
#' Need to make sure that the get_all_models function is not missing any models that should be included.
#' 
#' ------------------------------------------
require(phylopath)
require(dplyr)
require(plyr)
#require(foreach)
#require(doParallel)
#require(parallel)

#n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
#registerDoParallel(cores = (n_cores-2))

#' ------------------------------------------
#' get_all_models function
#' 
#' The get_all_models function generates a list of all possible acylic models, given a set of variables.
#' This function returns a data.frame in which there is a column for each of the variables specified.
#' The rows in this data.frame are 
#' Using the "exclusions" argument, it is possible to exclude all models in which on variable is predicted by another.

get_all_models <- function(variable_list, exclusions){
  
  #' ----------------
  # Some basic checkes of the variables supplied
  
  if(sum(duplicated(variable_list))>0){
    stop("Each variable name must be unique")
    
  }else if(sum(grepl(" ", variable_list)) > 0){
    stop("Variable names cannot contain spaces")
    
  } else if(sum(grepl("~", variable_list)) > 0){
    stop("Variable names cannot contain the '~' symbol")
    
  }
  
  #' ----------------
  #' Defining sub-functions
  
  tidy_strings <- function(x){
    x_items <- strsplit(x, " ") %>% unlist()
    x_items[!x_items %in% c("~", "+")]
  }
  
  #' ----------------
  #' Number of variables
  n_vars = length(variable_list)
  
  #' ----------------
  #' Defining all possible components of a path model
  combos <- vector()
  
  for(i in 1:(n_vars - 1)){
    i_combos <- combn(variable_list, m = i, simplify = F)
    combos <- c(combos, i_combos)
  }
  
  every_formulae <- list()
  
  #' This loop goes through each variable in the variable_list and defines every possible combination of predictor variables
  for(i in 1:n_vars){
    
    # Defining the response variable
    i_var <- variable_list[i]
    
    # Creating a blank list, to contain every possible combination of predictor variables for the response variable
    i_formulae <- list()
    
    # This loop goes through the list of combos and selects only those in which the response is not included (unless the response is the only variable).
    for(j in 1:length(combos)){
      
      #' Need to include the possibility that the variable is not predicted by any other variables.
      if(length(combos[[j]]) == 1 & sum(i_var == combos[[j]]) == 1){
        j_formula <- i_var 
        i_formulae <- c(i_formulae, j_formula)
        
        #' This selects as predictors every combination of variables that do not include the response variable.   
      }else if(! i_var %in% combos[[j]]){
        j_formula <- paste0(i_var, " ~ ", paste0(unlist(combos[j]), collapse =  " + "))
        i_formulae <- c(i_formulae, j_formula)
      }
      
    }
    
    every_formulae <- c(every_formulae, i_formulae)
    
  }
  
  # ------------
  #' Dropping any excluded components of the path model
  
  # If exclusions have been specified, then drop them from the every_formulae list
  if(!missing(exclusions)){
    
    every_formulae_characters <- lapply(every_formulae, tidy_strings)
    
    exclusions_characters <- lapply(as.character(exclusions), tidy_strings)
    
    # This provides an index of whether the variable is to be included or not  
    every_formulae_include <- rep(T,length(every_formulae_characters))
    
    for(i in 1:length(every_formulae_characters)){
      formula_characters_i <- every_formulae_characters[[i]]
      
      if(length(formula_characters_i) > 1){
        i_response <- formula_characters_i[1]
        i_predictors <- formula_characters_i[2:length(formula_characters_i)]
        
        for(e in 1:length(exclusions)){
          e_response <- exclusions_characters[[e]][1]
          e_predictors <- exclusions_characters[[e]][2:length(exclusions_characters[[e]])]
          
          if(i_response == e_response & e_predictors %in% i_predictors){
            every_formulae_include[i] <- FALSE
            #' Note, this will continue running through the list of exclusions, even if one has been met. 
            #' Should probably jump to the end of the j and i iterations at this point
          }
          
        }
        
      }
      
    }
    
    every_formulae <- every_formulae[every_formulae_include]
    
  }
  
  # ------------
  # Function to test whether the formulae is circular
  
  model_filter <- function(formulae_list){
    f <- formulae_list[ ! is.na(formulae_list)]
    
    tidy_strings <- function(x){
      x_items <- strsplit(x, " ") %>% unlist()
      x_items[!x_items %in% c("~", "+")]
    }
    
    f_tidy <- lapply(f, tidy_strings)
    
    var_matrix <- data.frame(matrix(nrow = n_vars,
                                    ncol = n_vars))
    colnames(var_matrix) <- variable_list
    row.names(var_matrix) <- variable_list
    
    var_matrix[,] <- 0
    
    for(f_i in 1:length(f_tidy)){
      if(length(f_tidy[[f_i]]) > 1){
        var_col <- f_tidy[[f_i]][1]
        var_rows <- f_tidy[[f_i]][2:length(f_tidy[[f_i]])]
        var_matrix[var_rows, var_col] <- 1 
      }
      
    }
    
    # By default the function will output False
    output <- F
    
    # If the var_matric is entirely full of 0s, then this model will also be excluded.
    # NOTE: THIS SHOULD ONLY BE APPLIED AT THE END
    #if(sum(rowSums(var_matrix)) > 0){
    
    tryCatch({
      
      #' This is a bit of a hack, but uses the basiSet function in the ggm package to identify whether the var_matrix 
      #' is acyclic (tests whether there is circularity in the predictor variable structure)
      #' This function will fail when the matrix is not acyclic, and the tryCatch function is used to prevent this
      #' from stopping the function.
      #' NOTE: If I wanted to develop this further, I could potentially find the source code the for basiSet function
      #' and adapt it for the current purposes. 
      f_test <- ggm::basiSet(as.matrix(var_matrix))
      
      if(length(f_test) > 1){
        output <- T
      }
      
    }, error=function(e){})
    
    #    }
    
    # This outputs TRUE if the model is acyclic and FALSE if the model contains circularity.
    output
    
  }
  
  # ------------
  # Setting up the first variable combos
  
  first_var_formulae <- every_formulae[grepl(paste0("^", variable_list[1]), every_formulae)] %>% 
    unlist()
  all_formulae_combos <- data.frame(matrix(nrow = length(first_var_formulae), 
                                           ncol = n_vars))
  all_formulae_combos[,1] <- first_var_formulae
  colnames(all_formulae_combos) <- variable_list
  
  # Using a loop to identify all other usable path model structures
  for(i in 2:n_vars){
    
    # Defining the response variable for this iteration
    i_var <- variable_list[i]
    
    # Select all formulae that predict i_var
    i_formulae <- every_formulae[grepl(paste0("^", i_var), every_formulae)] %>% unlist()
    
    # Number of rows at the start of this iteration
    n_rows_start <- nrow(all_formulae_combos)
    
    # Duplicate each existing formulae by the number of formulae in which i_var is predicted
    all_formulae_combos <- all_formulae_combos[rep(1:n_rows_start, each = length(i_formulae)),]
    
    
    all_formulae_combos[,i_var] <- rep(i_formulae, times = n_rows_start)
    
    # This is the part of this loop that takes a while. Probably need to work out an efficient way to make this parallel.
    combos_to_include <- apply(all_formulae_combos, 1, model_filter)
    
    # This line below should parallelise the function, but won't work on Windows machines
    #combos_to_include <- mclapply(alply(all_formulae_combos, .margins = 1), model_filter, mc.cores = n_cores - 2)
    
    all_formulae_combos <- all_formulae_combos[combos_to_include, ]
    
  }
  
  
  # There should be one row in which no variables are predicted by any other variables. This needs to be dropped. 
  r_index <- NULL
  
  for(r in 1:nrow(all_formulae_combos)){
    r_index[r] <- ! sum(as.character(all_formulae_combos[r, ]) == colnames(all_formulae_combos)) == ncol(all_formulae_combos)
  }
  
  all_formulae_combos <- all_formulae_combos[r_index,]
  
  all_formulae_combos
  
}

#' ------------------------------------------
#' strings_to_model_sets function
#' 
#' This function converts the strings generated in the get_all_models function into 

# Requires the columns to be named after each variable (should be)
strings_to_model_sets <- function(list_model_strings){
  
  variables <- colnames(list_model_strings)
  
  # Function to covert strings to fomula
  string2formula <- function(x){
    x[!x == names(x)] %>% 
      as.list() %>%
      lapply(as.formula)
  }
  
  # Converting strings to formulae
  list_model_formulae <- apply(list_model_strings, MARGIN = 1, string2formula)
  
  # List to become the list of model matricies
  all_model_matricies <- list()
  
  for(i in 1:length(list_model_formulae)){
    all_model_matricies[i] <- define_model_set(list_model_formulae[[i]])
    # If there are variables missing from the matrix, add them in here.
    if(nrow(all_model_matricies[i][[1]]) < length(variables)){
      missing_variables <- variables[!variables %in% row.names(all_model_matricies[i][[1]])]
      
      for(j in 1:length(missing_variables)){
        all_model_matricies[i][[1]] <- rbind(all_model_matricies[i][[1]], 0)
        row.names(all_model_matricies[i][[1]])[nrow(all_model_matricies[i][[1]])] <- missing_variables[j]
        
        all_model_matricies[i][[1]] <- cbind(all_model_matricies[i][[1]], 0)
        colnames(all_model_matricies[i][[1]])[ncol(all_model_matricies[i][[1]])] <- missing_variables[j]      
        
      }
    }
  }
  
  names(all_model_matricies) <- paste0("m", 1:length(all_model_matricies))
  
  all_model_matricies
  
}

