

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

get_failed_models <- function(taxa, species_names, pseudoabs_type, index_file, models = c('lr', 'gam', 'rf', 'me')) {
  
  require(purrr)
  
  # taxa = taxa
  # species_names = sp_name 
  # pseudoabs_type = pa_name
  # index_file = file_for_lotus
  
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
                                     name_index_mods = index_file$X[paste0(gsub(' ', replacement = '_', x=index_file$species), index_file$model) %in% ## use the file for lotus to find the lotus number of the file of interest
                                                                      paste0(species_names$name_df$species[name_index], models[m])]-1, # -1 to match the slurm output (zero-indexed)
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
##' Returns a list of all the error files

get_error_logs <- function(taxa, failed_models, pseudoabs_type, models = c('lr', 'gam', 'rf', 'me'), slurm_location) {
  
  # taxa=taxa
  # failed_models = failed_mods
  # pseudoabs_type = pa_name
  # models = c('lr', 'gam', 'rf', 'me')
  
  require(purrr)
  require(tidyverse)
  
  spp_out <- lapply(c(1:length(failed_models$name_index_mods)), FUN = function(n){
    
    # I want to find the name_index of the species
    # then search for that index in the error outputs
    error_location <- list.files(slurm_location, 
                                 pattern = paste0('-', failed_models$name_index_mods[n], '.err'),
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
      
      cropped_log <- error_log
      
    } else if(length(grep(pattern="error|Error", error_log))==0){
      
      warning(paste('!!!   No error message in log file for species', 
                    failed_models$species[n], 'and model',
                    failed_models$model[n]))
      cropped_log <- error_log 
      
    }
  
    names(error_log) <- paste(failed_models$species[n], failed_models$model[n], sep = '_')
    
    return(error_log)
    
  })
  
  # lapply statement to get the errors on a single line
  err <- lapply(1:length(spp_out), FUN = function(x){
    er_lines <- grep(pattern="error|Error", spp_out[[x]])
    return(spp_out[[x]][er_lines])
  })
  
  failed_models$error <- do.call(c, err)
  
  # name the entries of the list
  names(spp_out) <- failed_models$species
  
  return(list(error_logs = spp_out, dataframe = failed_models))
  
}





#####      Workflow for checking

# taxa
taxa = 'butterfly'

# pseudoabsence type (i.e. the model run)
pa_name = 'PA_thinned_10000nAbs'

# set working directory to rslurm folder 
# to be able to access the file_for_lotus
fls_loc <- paste0('/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/scripts/',taxa, '/_rslurm_lrgamrfme_', pa_name)
setwd(fls_loc)

# read in file for lotus
file_for_lotus <- read.csv("file_for_lotus.csv")

# get the species names
sp_name <- get_spp_names(taxa, pa_name)

# check each species for failed models
failed_mods <- get_failed_models(taxa, sp_name, pa_name, index_file = file_for_lotus)
failed_mods

# get the associated error logs for each species
error_logs <- get_error_logs(taxa, failed_mods, pa_name, slurm_location = fls_loc)
error_logs # manually look at what the errors were and decide which ones need to be rerun.




#####     resubmit scripts that failed 

# get the lines in file_for_lotus that failed
file_for_lotus_resub <- file_for_lotus[file_for_lotus$X %in% failed_mods[1:12,]$name_index_mods,] ## change the indexing numbers for each error set
write_csv(file_for_lotus_resub, 'file_for_lotus_resub.csv')

#####    Automated lotus script
library(rslurm)
library(tidyverse)

# taxa
taxa=taxa

## choose model of interest
model = c('lr', 'gam', 'rf', 'me')

# pseudoabsence name
pa_name = pa_name

# number of bootstraps
k = 10

# number of knots for gam
knots_gam = 4


#####    arguments for slurm_apply     #####

# queue for lotus
queue_name = 'long-serial'

# time requirement
time = '47:59:59'

# memory requirement in MB
mem_req = 40000


# get the parameters for function
pars <- data.frame(name_index = seq(1, dim(file_for_lotus_resub)[1]))


## the function
slurm_sdm_boot <- function(name_index) {
  
  require(raster)
  require(dismo)
  require(tidyverse)
  require(mgcv)
  require(randomForest)
  require(glmnet)
  source("/home/users/thoval/DECIDE/scripts/Edited_Rob_Functions.R")
  
  # read the file_for_lotus.csv
  file_for_lotus <- read.csv('../file_for_lotus_resub.csv')
  
  
  #####     setting the parameters from the file for lotus     #####
  
  # taxa
  taxa = file_for_lotus$taxa[name_index] #'moth' # moth, butterfly
  
  ## choose model of interest
  model = file_for_lotus$model[name_index] # 'rf' # one of c('lr', 'gam', 'rf', 'me')
  warning(paste('!!!   ',model))
  
  # pseudoabsence name
  pa_name = file_for_lotus$pa_name[name_index] # 'PA_thinned_10000nAbs'
  
  # number of bootstraps
  k = file_for_lotus$k[name_index]
  
  # number of knots for gam
  knots_gam = file_for_lotus$knots_gam[name_index]
  
  
  
  # where to save at the end
  outPath = paste0("/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/outputs/", taxa, "/SDM_Bootstrap_", model, "_", pa_name,"/")
  
  # # load file with parameters
  # file_for_lotus = read.csv(paste0("/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/scripts/", taxa, "/_rslurm_", model, "_", pa_name, "/file_for_lotus.csv"))
  
  # load environmental data
  env_dat <- raster::stack("/home/users/thoval/DECIDE/data/environmental_data/envdata_fixedcoasts_nocorrs_100m_GB.gri")
  
  
  # load pseudoabsences based on taxa and name
  if(taxa == 'moth'){
    
    load(paste0("/home/users/thoval/DECIDE/data/species_data/moths/", taxa, '_', pa_name, ".rdata"))
    
  } else if(taxa == 'butterfly'){
    
    load(paste0("/home/users/thoval/DECIDE/data/species_data/butterflies/", taxa, '_', pa_name, ".rdata"))
    ab1 <- res_out
    
  }
  
  ## Find the species of interest
  species = file_for_lotus$species[name_index]
  warning(paste('!!!   species   !!!  ', species, '  !!!   species   !!!'))
  
  # run the model
  sdm <- fsdm(species = species, model = model, 
              climDat = env_dat, spData = ab1, knots_gam = knots_gam,
              k = k, # number of bootstraps
              write = F, outPath = outPath)
  
  # choose the type and index for predict function
  if(model == 'lr'|model == 'gam'){
    type <- "response"
    index <- NULL
  } else if(model == 'rf'){
    type <- "prob"
    index <- 2
  }
  
  ## bootstrapped models
  print(paste0('#####   predicting from bootstrapped models   #####')) 
  
  ## predict from each of the bootstrapped models
  ## different workflow for lrReg and other methods
  mod_preds <- get_predictions(model_outs = sdm,
                               model = model, 
                               env_data = env_dat)
  
  
  ## save files ##
  print("#####     Saving files     #####")
  species_name <- gsub(pattern = ' ', replacement = '_', species) # get species name without space
  
  # save prediction raster
  print("#####     Saving prediction raster     #####")
  writeRaster(x = mod_preds$mean_predictions, 
              filename = paste0(outPath, model, "_SDMs_", species_name, "_meanpred.grd"),
              format = 'raster', overwrite = T)
  
  # # save quantile max min
  # print("#####     Saving quantile max min raster     #####")
  # writeRaster(x = mod_preds$quant_minmax, 
  #             filename = paste0(outPath, model, "_SDMs_", species_name, "_quantilemaxmin.grd"),
  #             format = 'raster', overwrite = T)
  
  # save standard deviation raster
  print("#####     Saving standard deviation raster     #####")
  writeRaster(x = mod_preds$quant_range, 
              filename = paste0(outPath, model, "_SDMs_", species_name, "_quantilerange.grd"),
              format = 'raster', overwrite = T)
  
  # write AUC to file for easy-access
  print("#####     Writing AUC to file     #####")
  write.csv(x = data.frame(raw_AUC = sdm$AUC,
                           meanAUC = sdm$meanAUC),
            file = paste0(outPath, model, "_SDMs_", species_name, "_AUC_values.csv"))
  
  # write data to file too
  print("#####     Writing data to file     #####")
  write.csv(x = sdm$Data,
            file = paste0(outPath, model, "_SDMs_", species_name, "_Data.csv"))
  
  # save subset model output
  print("#####     Saving model output     #####")
  
  # remove data from model output
  sdm$Data <- NULL
  
  # outout of model to store
  model_output <- list(species = species_name,
                       model = model,
                       sdm_output = sdm,
                       number_validations = k)
  
  save(model_output, file = paste0(outPath, model, "_SDMs_", species_name, 
                                   ".rdata"))
  
  return(model_output)
  
}

#### slurm apply call
sdm_slurm <- slurm_apply(slurm_sdm_boot,
                         params = pars,
                         jobname = paste0('resub_',paste0(model, collapse = ''), '_', pa_name),
                         nodes = dim(pars)[1],
                         cpus_per_node = 1,
                         slurm_options = list(partition = queue_name,
                                              time = as.character(time),
                                              mem = mem_req,
                                              error = paste0(paste(model, collapse = '_'),'_error-%j-%a.err')),
                         rscript_path = '',
                         submit = T)

get_job_status(sdm_slurm)