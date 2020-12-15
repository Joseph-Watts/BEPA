#' Running without parallel processing

start_time <- Sys.time()
get_all_parF_1 <- get_all_models(variable_list = c("Aa11", 
                                 "Bb22", 
                                 #"Cc33", 
                                 "Dd44", 
                                 "Ee55"))
pF_1_time <- Sys.time() - start_time
pF_1_time

get_all_parF_2 <- get_all_models(variable_list = c("Aa11", 
                                                   "Bb22", 
                                                   #"Cc33", 
                                                   "Dd44", 
                                                   "Ee55"),
                                 exclusions = c("Aa11 ~ Bb22", 
                                                "Dd44 ~ Bb22") ,
                                 parallel = F,
                                 max_variables_in_models = 3,
                                 required_variables = "Ee55",
                                 n_cores = NULL)


#' ----
#' Multicore functions


start_time <- Sys.time()
get_all_4core_1 <- get_all_models(variable_list = c("Aa11", 
                                                     "Bb22", 
                                                     #"Cc33", 
                                                     "Dd44", 
                                                     "Ee55"),
                                   parallel = T,
                                   n_cores = 4)
c4_1_time <- Sys.time() - start_time
c4_1_time


get_all_4core_2 <- get_all_models(variable_list = c("Aa11", 
                                                   "Bb22", 
                                                   #"Cc33", 
                                                   "Dd44", 
                                                   "Ee55"),
                                 exclusions = c("Aa11 ~ Bb22", 
                                                "Dd44 ~ Bb22") ,
                                 parallel = T,
                                 max_variables_in_models = 3,
                                 required_variables = "Ee55",
                                 n_cores = 4)

start_time <- Sys.time()
get_all_8core_1 <- get_all_models(variable_list = c("Aa11", 
                                                   "Bb22", 
                                                   #"Cc33", 
                                                   "Dd44", 
                                                   "Ee55"),
                                 parallel = T,
                                 n_cores = 8)
c8_1_time <- Sys.time() - start_time
c8_1_time






#' strings_to_model_sets
testing_1 <- strings_to_model_sets(get_all_parF_1,
                      parallel = T,
                      n_cores = 4)


testing_2 <- strings_to_model_sets(get_all_8core_1)


