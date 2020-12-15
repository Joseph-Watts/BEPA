#' Brute-force Exploratory Path Analysis (BEPA)
#' 
#' Joseph Watts - https://github.com/Joseph-Watts
#' 
#' The BEPA package contains functions to help perform exploratory path analysis in R.
#' 
#' The package can be used to identify all possible path model structures, given a set of 
#' variables. All of these models can then be compared. Time taken increases exponentially 
#' with the number of variables.
#' 
#' This package is written specifically for use with the phylopath package, which itself 
#' is based on the phyloglm and phylolm packages.
#' The functions should also be useable for exploratory path analyses more generally, but 
#' this has not been tested.
#' 
#' Currently this has only been tested on Windows 10. 
#' 
#' ------------------------------------------
#' Working notes:
#' 1. The exclusion function could be made more general
#' Would probably be good to be able to exclude a wider range of models in the future.
#' For example, might want to exclude multiple predictors of a single variable, or 
#' require a variable to be predicted by others.
#'  
#' 2. When a variable is not predicted by any other variables the get_all_models just 
#' puts the variable name in the variable column. This makes it easier to identify,
#' but cannot be directly converted into a formula. It might be worth representing such
#' cases as predicted by themselves (e.g. "a ~ a") as this is how they are represented
#' when converting to DAGs. Not sure.
#' 
#' 3. Might be worth putting warning messages in when too many variables are specified.
#' Also, when there are too many cores specified. As core count increases so does total
#' memory allocation required. Memory can be a bottle neck when a large (>5) number 
#' of variables are included.
#' 
#' 8. To get around RAM constraints, it might be possible to run the function in sections,
#' save the outputs to a drive, then at the end combine them together? 
#' 
#' Possible future directions:
#' 1. Rather than having to the function anew every time, maybe it is worth saving an index
#' and then just reading it back in and then filtering it? 
#' 
#' ------------------------------------------
require(phylopath)
require(tidyverse)
require(gsubfn)
require(doParallel)
require(parallel)

#' ----------------
#' Defining sub-functions

#' ----
#' Alternative to gsubfn with arguments ordered for ease of use with apply
alt_gsubfn <- function(x, 
                       pattern,
                       replacement){
  
  gsubfn(pattern = pattern,
         replacement = replacement,
         x = x) %>%
    return()
  
}

#' ----
#' Function to reformat and tidy formulae text strings
tidy_strings <- function(x){
  x_items <- strsplit(x, " ") %>% unlist()
  x_items[!x_items %in% c("~", "+")] %>%
    return()
}

#' ----
#' Function to test whether formulae contain circularity
model_filter <- function(formulae_list,
                         m_variables){
  
  m_n_vars <- length(m_variables)
  
  f <- formulae_list[ ! is.na(formulae_list)]
  
  tidy_strings <- function(x){
    x_items <- strsplit(x, " ") %>% unlist()
    x_items[!x_items %in% c("~", "+")]
  }
  
  f_tidy <- lapply(f, tidy_strings)
  
  var_matrix <- data.frame(matrix(nrow = m_n_vars,
                                  ncol = m_n_vars))
  colnames(var_matrix) <- m_variables
  row.names(var_matrix) <- m_variables
  
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
  
  tryCatch({
    
    #' This is a bit of a hack, but uses the basiSet function in the ggm package 
    #' to identify whether the var_matrix is acyclic (tests whether there is 
    #' circularity in the predictor variable structure)
    #' This function will fail when the matrix is not acyclic, and the tryCatch 
    #' function is used to prevent this from stopping the function.
    #' NOTE: If I wanted to develop this further, I could potentially find the 
    #' source code the for basiSet function and adapt it for the current purposes. 
    
    f_test <- ggm::basiSet(as.matrix(var_matrix))
    
    if(length(f_test) > 1){
      output <- T
    }
    
  }, error=function(e){})
  
  #    }
  
  #' This outputs TRUE if the model is acyclic and FALSE if the model contains 
  #' circularity.
  return(output)
  
}


x_not_in_y <- function(x, pattern){
  
  grepl(pattern = pattern,
        x = x) %>%
    sum() == 0 %>%
    return()
  
}

#' Function to covert strings to fomula
string2formula <- function(x){
  
  x <- as.list(x) %>% unlist()
  
  #' If there are any isolated variables, they need to be redefined as a circular formula.
  #' This is how they need to be formatted for the phylopath::define_model_set function
  if(any(x == names(x))){
    
    isolated_variables <- x[x == names(x)]
    
    for(i in 1: length(isolated_variables)){
      
      x[isolated_variables[i]] <- paste(isolated_variables[i], "~", isolated_variables[i])
      
    }
    
  }
  
  # Format as a list
  as.list(x) %>%
    
    # Convert to formula
    lapply(as.formula)
  
}


#' ------------------------------------------
#' get_all_models function
#' 
#' The get_all_models function generates a list of all possible acyclic models, 
#' given a set of variables.
#' This function returns a data.frame in which there is a column for each of the 
#' variables specified. The entire rows in this data.frame represents a model.
#' Using the "exclusions" argument it is possible to exclude all models in which 
#' one variable is predicted by another.


get_all_models <- function(variable_list, 
                           exclusions = NULL, 
                           parallel = F,
                           max_variables_in_models = NULL,
                           required_variables = NULL, #' Note: this argument only functions when max_variables_in_model 
                           #' is less than length(variable_list). Might be worth giving a warning when this arguments
                           #' is specified but not actually doing anything.
                           n_cores = NULL){
  
  testing = F
  
  if(testing){
    
    variable_list = c("Aa11", 
                      "Bb22", 
                      #"Cc33", 
                      "Dd44", 
                      "Ee55")
    
    exclusions = c("Aa11 ~ Bb22", 
                   "Dd44 ~ Bb22") 
    #exclusions = NULL
    parallel = F
    #max_variables_in_models = 4
    max_variables_in_models = NULL
    required_variables = NULL 
    #required_variables = "Ee55"
    #n_cores = 4
    n_cores = NULL
    
  }
  
  
  # If n_cores is set to more than 1 but parallel is not True, give a warning
  if((!is.null(n_cores)) & parallel == F){
    
    
    if(n_cores > 1){
      
      warning("n_cores can only be greater than 1 when parallel is set to T: running on a single thread")
      
    }
    
  }
  
  
  #' Checking that n_cores argument is valid if parallel is true
  if(parallel){
    
    total_system_threads <- detectCores(all.tests = FALSE, logical = TRUE)
    
    if(is.null(n_cores)){
      
      n_cores <- total_system_threads - 4
      
      #' If the system has less than 5 logical cores, and the number of cores is not specified, 
      #' turn off parallel
      if(n_cores < 2){
        
        parallel = F
        
        print("parallel set to TRUE, n_cores not defined, & system has less than 5 core: turning parallel off")
        
      }else{
        
        print(paste0("parallel set to TRUE, n_cores not defined: n_cores set to ", n_cores))
        
      }
      
    }else if(n_cores > total_system_threads){
      
      n_cores <- total_system_threads
      
      warning(paste0("n_cores argument set to a greater number of threads than available in system: setting n_cores to ", 
                     n_cores))
      
    }else if(n_cores < 1){
      
      n_cores <- 1
      
      warning(paste0("n_cores argument set to less than 1: setting n_cores to ", 
                     n_cores))
      
    }
    
  }
  
  #' If the parallel argument is true, set up a cluster
  if(parallel){

    #' Make cluster
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)  
    
    #' Load the packages within the cluster
    clusterEvalQ(cl, library("tidyverse"))
    clusterEvalQ(cl, library("gsubfn"))
    
  }
  
  
  #' ----------------
  #' Some basic checks of the variables supplied
  
  if(sum(duplicated(variable_list)) > 0){
    
    stop("Each variable name must be unique")
    
  }else if(sum(grepl(" ", variable_list)) > 0){
    
    stop("Variable names cannot contain spaces")
    
  } else if(sum(grepl("~", variable_list)) > 0){
    
    stop("Variable names cannot contain the '~' symbol")
    
  }
  
  #' Converting variable names into single letters
  #' This saves systems memory and speeds things up a little
  variables <- letters[1:length(variable_list)]
  
  if(!is.null(exclusions)){
    
    #' Loop to replace variable names in exclusions
    for(i in 1:length(variables)){
      
      if(length(grep(variable_list[i], variable_list)) > 1){
        stop("Variable names embeded within other names not allowed") 
        #' this is required for the variable character length reductions
      }
      
      exclusions <- gsub(variable_list[i], variables[i], exclusions)
      
    }
    
  }
  
  if(!is.null(required_variables)){
    
    #' Loop to replace variable names in exclusions
    for(i in 1:length(variables)){
      
      if(length(grep(variable_list[i], variable_list)) > 1){
        stop("Variable names embeded within other names not allowed") 
        #' this is required for the variable character length reductions
      }
      
      required_variables <- gsub(variable_list[i], variables[i], required_variables)
      
    }
    
  }
  
  
  
  #' ----------------
  #' Number of variables
  n_vars = length(variables)
  
  #' ----------------
  #' Defining all possible components of a path model
  combos <- vector()
  
  if(is.null(max_variables_in_models)){
    
    max_predictors <- n_vars - 1
    
  } else {
    
    max_predictors <- max_variables_in_models - 1
    
  }
  
  
  for(i in 1:max_predictors){
    
    i_combos <- combn(variables, m = i, simplify = F)
    
    combos <- c(combos, i_combos)
    
  }
  
  
  every_formulae <- list()
  
  #' This loop goes through each variable in variables and defines every 
  #' possible combination of predictor variables
  for(i in 1:n_vars){
    
    #' Defining the response variable
    i_var <- variables[i]
    
    #' Creating a blank list, to contain every possible combination of predictor 
    #' variables for the response variable
    i_formulae <- list()
    
    #' This loop goes through the list of combos and selects only those in which 
    #' the response is not included (unless the response is the only variable).
    for(j in 1:length(combos)){
      
      #' Need to include the possibility that the variable is not predicted by any 
      #' other variables.
      if(length(combos[[j]]) == 1 & sum(i_var == combos[[j]]) == 1){
        j_formula <- i_var 
        i_formulae <- c(i_formulae, j_formula)
        
        #' This selects as predictors every combination of variables that do not 
        #' include the response variable.   
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
  if(!is.null(exclusions)){
    
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
            #' Note, this will continue running through the list of exclusions, 
            #' even if one has been met. 
            #' It would be more efficient to jump to the end of the j and i loops at this point?
          }
          
        }
        
      }
      
    }
    
    every_formulae <- every_formulae[every_formulae_include]
    
  }
  
  #' ------------
  #' Setting up the first variable combos
  
  first_var_formulae <- every_formulae[grepl(paste0("^", variables[1]), every_formulae)] %>% 
    unlist()
  
  all_formulae_combos <- data.frame(matrix(nrow = length(first_var_formulae), 
                                           ncol = n_vars))
  
  all_formulae_combos[,1] <- first_var_formulae
  
  colnames(all_formulae_combos) <- variables
  
  #' Using a loop to identify all other usable path model structures
  #' This is where most of the time is used 
  for(i in 2:n_vars){
    
    # Defining the response variable for this iteration
    i_var <- variables[i]
    
    # Select all formulae that predict i_var
    i_formulae <- every_formulae[grepl(paste0("^", i_var), every_formulae)] %>% unlist()
    
    # Number of rows at the start of this iteration
    n_rows_start <- nrow(all_formulae_combos)
    
    # Duplicate each existing formulae by the number of formulae in which i_var is predicted
    all_formulae_combos <- all_formulae_combos[rep(1:n_rows_start, each = length(i_formulae)),]
    
    all_formulae_combos[,i_var] <- rep(i_formulae, times = n_rows_start)
    
    #' This is the part of this loop that takes a while. 
    
    if(parallel){
      
      clusterExport(cl= cl, 
                    varlist = ls(), 
                    envir = environment())
      
      combos_to_include <- parApply(cl = cl,
                                    X = all_formulae_combos, 
                                    MARGIN = 1, 
                                    FUN = model_filter,
                                    m_variables = variables)
      
    }else{
      
      combos_to_include <- apply(X = all_formulae_combos, 
                                 MARGIN = 1, 
                                 FUN = model_filter,
                                 m_variables = variables)
      
    }
    
    all_formulae_combos <- all_formulae_combos[combos_to_include, ]
    
  }
  
  #' There should be one row in which no variables are predicted by any other variables. 
  #' This is dropped here.
  r_index <- NULL
  
  for(r in 1:nrow(all_formulae_combos)){
    r_index[r] <- ! sum(as.character(all_formulae_combos[r, ]) == colnames(all_formulae_combos)) == ncol(all_formulae_combos)
  }
  
  all_formulae_combos <- all_formulae_combos[r_index,]
  
  #' Quick tidy, these names will be replaced later
  rownames(all_formulae_combos) <- c()
  
  #' Dropping combos that do not include required variables
  if(!is.null(required_variables)){
    
    for(i in 1:length(required_variables)){
      
      i_col_index <- which(colnames(all_formulae_combos) == required_variables[i])
      
      #' Index of all rows in which the required variable is not predicted by any other variables
      i_not_response <- which(all_formulae_combos[, i_col_index] == required_variables[i])
      
      #' Index of all rows in which i is not a predictor variable
      i_not_predictor <- apply(X = all_formulae_combos[ , -i_col_index], 
                               MARGIN = 1,
                               FUN = x_not_in_y,
                               required_variables[i]) %>% 
        which()
      
      test <- all_formulae_combos[- i_not_response, ]
      test <- test[- i_not_predictor, ]
      
      #' Drop rows if the required variable is neither a response variable or a predictor variable
      i_drop <- intersect(i_not_response, i_not_predictor)
      
      all_formulae_combos <- all_formulae_combos[- i_drop, ]
      
    }
    
  }
  
  #' Renaming the model columns
  row.names(all_formulae_combos) <- paste0("m", 1:nrow(all_formulae_combos))
  
  #' Replacing the single letter variable names with the full variable names again
  if(parallel){
    
    all_formulae_combos <- parApply(cl = cl,
                                    X = all_formulae_combos,
                                    MARGIN = 2,
                                    FUN = alt_gsubfn,
                                    "\\S+",
                                    setNames(as.list(variable_list), variables)) %>%
      as.data.frame(stringsAsFactors = F)
    
    #' Stopping cluster
    stopCluster(cl)
    
  } else {
    
    all_formulae_combos <- apply(X = all_formulae_combos,
                                 MARGIN = 2,
                                 FUN = alt_gsubfn,
                                 "\\S+",
                                 setNames(as.list(variable_list), variables)) %>%
      as.data.frame(stringsAsFactors = F)
    
  }
  
  colnames(all_formulae_combos) <- variable_list
  
  #' Doneskies
  return(all_formulae_combos)
  
}



#' ------------------------------------------
#' strings_to_model_sets function
#' 
#' This function converts the strings generated in the get_all_models 
#' function into matrices representing directed acyclic graphs (DAGs).
#' Takes the whole data.frame as the input.
#' 
#' Requires the columns to be named after each variable, which should 
#' be the case if the get_all_models function is used


strings_to_model_sets <- function(model_strings,
                                  parallel = F,
                                  n_cores = NULL){
  
  testing <- F
  
  if(testing){
    
    model_strings <- get_all_models(variable_list = c("Aa11", 
                                                      "Bb22", 
                                                      #"Cc33", 
                                                      "Dd44", 
                                                      "Ee55"))
    
    parallel = T
    n_cores = NULL
    
  }
  
  #' -----
  #' Quick checks of arguments
  
  # If n_cores is set to more than 1 but parallel is not True, give a warning
  if((!is.null(n_cores)) & parallel == F){
    
    
    if(n_cores > 1){
      
      warning("n_cores can only be greater than 1 when parallel is set to T: running on a single thread")
      
    }
    
  }
  
  
  #' Checking that n_cores argument is valid if parallel is true
  if(parallel){
    
    total_system_threads <- detectCores(all.tests = FALSE, logical = TRUE)
    
    if(is.null(n_cores)){
      
      n_cores <- total_system_threads - 4
      
      #' If the system has less than 5 logical cores, and the number of cores is not specified, 
      #' turn off parallel
      if(n_cores < 2){
        
        parallel = F
        
        print("parallel set to TRUE, n_cores not defined, & system has less than 5 core: turning parallel off")
        
      }else{
        
        print(paste0("parallel set to TRUE, n_cores not defined: n_cores set to ", n_cores))
        
      }
      
    }else if(n_cores > total_system_threads){
      
      n_cores <- total_system_threads
      
      warning(paste0("n_cores argument set to a greater number of threads than available in system: setting n_cores to ", 
                     n_cores))
      
    }else if(n_cores < 1){
      
      n_cores <- 1
      
      warning(paste0("n_cores argument set to less than 1: setting n_cores to ", 
                     n_cores))
      
    }
    
  }
  
  #' If the parallel argument is true, set up a cluster
  if(parallel){

    #' Make cluster
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)  
    
    #' Load the packages within the cluster
    clusterEvalQ(cl, library("tidyverse"))
  }
  
  
  #' ------------
  #' Defining variables
  #' Taking the variable names from the data frame 
  variables <- colnames(model_strings)
  
  
  #' ------------
  #' Converting strings to formulae
  list_model_formulae <- list()
  
  # List to become the list of model matrices
  all_model_matrices  <- list()
  
  if(parallel){
    
    list_model_formulae <- parApply(cl = cl,
                                    X = model_strings,
                                    MARGIN = 1,
                                    FUN = string2formula)
    
    
    all_model_matrices <- parSapply(cl = cl,
                                    X = list_model_formulae,
                                    FUN = define_model_set)
    
    #' stopping cluster
    stopCluster(cl)
    
  }else{
    
    list_model_formulae <- apply(X = model_strings,
                                 MARGIN = 1,
                                 FUN = string2formula)
    
    all_model_matrices <- sapply(X = list_model_formulae,
                                 FUN = define_model_set)

  }
  
  names(all_model_matrices ) <- row.names(model_strings)
  
  return(all_model_matrices)
  
}
