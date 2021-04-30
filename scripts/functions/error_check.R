

### error checking functions and workflow for lotus


#####     Functions for error checking     #####

##' Get species names
##' 
##' Gets the names of all the species present in the pseudoabsences
##' file provided to the function so that we know which species we attempted to
##' run models for
##' This will only work for my directory structure and needs to have the right names inputted 
##' 

get_spp_names <- function(taxa, pseudoabs_type){
  
  require(tidyverse)
  
  # get names from the pseudoabsence data fed to the model
  if(taxa == 'moth'){
    
    load(paste0("/home/users/thoval/DECIDE/data/species_data/moths/moth_", pseudoabs_type, ".rdata")) ## moths
    name <- gsub(pattern = ' ', 
                 replacement = '_',
                 x = names(ab1))
    
    name_df <- data.frame(species = name, 
                          name_index = seq(0,length(name)-1))
    
  } else if(taxa == 'butterfly'){
    
    load(paste0("/home/users/thoval/DECIDE/data/species_data/butterflies/butterfly_", pseudoabs_type, ".rdata")) ## butterflies
    name <- gsub(pattern = ' ', 
                 replacement = '_',
                 x = names(res_out))
    
    name_df <- data.frame(species = name, 
                          name_index = seq(0,length(name)-1))
    
  } else {stop('!!! whoooaaaahhhh there boy, du calme! Name  not right !!!')}
  
  return(list(name = name, 
              name_df = name_df))
  
}


##' Check for failed models
##' 
##' Based on the output of the get_spp_names() function, find the models that failed
##' to produce any output.
##' Need to give 'species_names' argument the output of get_spp_names() function
##' returns a data frame with the failed models.

get_failed_models <- function(taxa, species_names, pseudoabs_type, models = c('lr', 'gam', 'rf', 'me')) {
  
  require(purrr)
  
  spp_out <- lapply(c(1:length(species_names$name)), FUN = function(name_index){
    
    er_mod <- lapply(c(1:length(models)), FUN = function(m){
      
      check_models <- list.files(paste0('/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/outputs/', 
                                        taxa, '/SDM_Bootstrap_', 
                                        models[m], '_', pseudoabs_type), 
                                 pattern = paste0(species_names$name[name_index]),
                                 full.names = TRUE)
      
      if(length(check_models)<=1){
        
        errored_models <- data.frame(taxa = taxa, 
                                     species = species_names$name_df$species[name_index], 
                                     name_index_mods = species_names$name_df$name_index[name_index], 
                                     model = models[m])
        
      } else { errored_models <- NULL }
      
      return(errored_models)
      
    })
    
    failed_models <- do.call('rbind', er_mod)
    
    return(failed_models)
    
  })
  
  failed_spp_mod <- do.call('rbind', spp_out)
  
  return(failed_spp_mod)
  
}


##' get error logs
##' 
##' Get the error logs associated with each of the failed models
##' The 'failed_models' argument takes the output of get_failed_models() function
##' Returns a list

get_error_logs <- function(taxa, failed_models, pseudoabs_type, models = c('lr', 'gam', 'rf', 'me')) {
  
  # taxa=taxa
  # failed_models = failed_mods
  # pseudoabs_type = pa_name
  # models = c('lr', 'gam', 'rf', 'me')
  
  require(purrr)
  require(tidyverse)
  
  spp_out <- lapply(c(1:length(failed_models$name_index_mods)), FUN = function(n){
    
    # get the models that failed for a species species
    failed_mods_list <- failed_models$model[n]
    
    er_mod <- lapply(c(1:length(failed_mods_list)), FUN = function(m){
      
      # I want to find the name_index of the species
      # then search for that index in the error outputs
      error_location <- list.files(paste0('/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/outputs/', 
                                          taxa, '/log_files/SDM_Bootstrap_', 
                                          failed_mods_list[m]), 
                                   pattern = paste0('-', failed_models$name_index_mods[n], '.out'),
                                   full.names = TRUE) %>% 
        sort(decreasing = TRUE)
      
      if(length(error_location) > 1){
        
        warning(paste("!!!   multiple matching arguments in error log folder for species",
                      failed_models$species[n], "and model",
                      failed_models$model[n], ". Index", failed_models$name_index_mods[n], 
                      "check to make sure that error log is correct and consider deleting old logs for each model run   !!!"))
        
      }
      
      # read the error log into R, skipping the NULL lines
      error_log <- readLines(error_location[1], skipNul = T)
      
      # find the lines that have "error or Error"
      # give context by including the 10 lines before the first error 
      # and to the last line in the file
      if(length(grep(pattern="error|Error", error_log))>0){
        
        er_lines <- grep(pattern="error|Error", error_log)
        
        cropped_log <- error_log[(min(er_lines)-10):length(error_log)]
        
      } else if(length(grep(pattern="error|Error", error_log))==0){
        
        warning(paste('!!!   No error message in log file for species', 
                      failed_models$species[n], 'and model',
                      failed_models$model[n]))
        cropped_log <- error_log 
        
      }
      
      return(cropped_log)
      
    })
    
    names(er_mod) <- paste(failed_models$species[n], failed_mods_list, sep = '_')
    
    return(er_mod)
    
  })
  
  # name the entries of the list
  names(spp_out) <- failed_models$species
  
  return(spp_out)
  
}




#####      Workflow for checking

# taxa
taxa = 'butterfly'

# pseudoabsence type (i.e. the model run)
pa_name = 'PA_thinned_10000nAbs'

# get the species names
sp_name <- get_spp_names(taxa, pa_name)

# check each species for failed models
failed_mods <- get_failed_models(taxa, sp_name, pa_name)

# get the associated error logs for each species
error_logs <- get_error_logs(taxa, failed_mods, pa_name)

error_logs
