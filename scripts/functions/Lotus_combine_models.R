## project directory
setwd("~/DECIDE/DECIDE_WP1")

###   Combine predictions from different models on lotus
library(tidyverse)
library(raster)

#####  Get parameters
# use the names from the pseudoabsences

# taxa for slurm output and parameter loading
# - still need to change taxa within function 
taxa = 'butterfly'
pseudoabs_type = '10000nAbs' ## same, still need to change within function when changing

if(taxa == 'moth'){
  
  ## for moths
  load('scripts/lotus/moth/pseudoabsences/moth_PA_unthinned_10000nAbs.rdata')
  pars <- data.frame(name_index = seq(1, length(names(ab1))))
  
} else if(taxa == 'butterfly'){
  
  ## for butterflies
  load('scripts/lotus/butterfly/pseudoabsence_scripts/butterfly_PA_unthinned_10000nAbs.rdata')
  pars <- data.frame(name_index = seq(1, length(names(res_out))))
  
}



calculate_ensemble <- function(name_index) {
  
  require(tidyverse)
  require(raster)
  
  ### 1. parameters for function
  taxa = 'butterfly' # moth
  pseudoabs_type = '10000nAbs' # which model run name to go through
  auc_cutoff = 0.75 ## just a suggestion - might need some thought (although AUC values so stupidly high might not be a problem until we're using a different score metric)
  models = c('lr', 'gam', 'rf', 'me') #, 'lrReg') ## lrReg hasn't worked for any species yet
  
  
  ## no need to change these unless changing directories
  main_directory = '/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/outputs'
  output_directory = paste0(main_directory, '/',taxa, '/combined_model_outputs/')
  
  
  ### 2. setup storage for model loops
  # error outputs
  error_out <- list()
  
  # model outputs
  mod_pred_out <- list()
  qrange_out <- list()
  auc_out <- list()
  errored_models <- list()
  
  
  ### 3. get names and species of interest 
  if(taxa == 'moth'){
    
    load("/home/users/thoval/DECIDE/data/species_data/moths/moth_PA_unthinned_10000nAbs.rdata") ## moths
    names <- gsub(pattern = ' ', 
                  replacement = '_',
                  x = names(ab1))
    
  } else if(taxa == 'butterfly'){
    
    load("/home/users/thoval/DECIDE/data/species_data/butterflies/butterfly_PA_unthinned_10000nAbs.rdata") ## butterflies
    names <- gsub(pattern = ' ', 
                  replacement = '_',
                  x =  names(res_out))
    
  } else {stop('whoooaaaahhhh there boy, du calme! Name  not right')}
  
  
  print(names[name_index])
  
  
  ### 4. Go through each model in turn
  for(m in 1:length(models)){
    
    ### first, check that the model for that species exists
    ### if it doesn't, skip to next model
    {check_models <- list.files(paste0(main_directory, '/', taxa, '/SDM_Bootstrap_', models[m], '_', pseudoabs_type), 
                                pattern = paste0(names[name_index]),
                                full.names = TRUE)
    
    if(length(check_models)<=1){
      
      print(paste('!!!   model', models[m], 'failed for species', names[name_index], '  !!!'))
      
      errored_models[[m]] <- data.frame(taxa = taxa, 
                                        species = names[name_index], 
                                        model = models[m])
      
      next
    }}
    
    
    ### load the mean predictions
    mp <- list.files(paste0(main_directory, '/', taxa, '/SDM_Bootstrap_', models[m], '_', pseudoabs_type), 
                     pattern = paste0(names[name_index], "_meanpred.grd"),
                     full.names = TRUE)
    
    mod_preds <- raster(mp)
    names(mod_preds) <- paste0(names[name_index], '_', models[m],'_mean_pred')
    mod_pred_out[[m]] <- mod_preds
    
    ### load the quantile range
    qr <- list.files(paste0(main_directory, '/', taxa, '/SDM_Bootstrap_', models[m], '_', pseudoabs_type), 
                     pattern = paste0(names[name_index], "_quantilerange.grd"),
                     full.names = TRUE)
    
    qrange <- raster(qr)
    names(qrange) <- paste0(names[name_index], '_', models[m], '_quantile_range')
    qrange_out[[m]] <- qrange 
    
    ### load the auc values
    aucval <- list.files(paste0(main_directory, '/', taxa, '/SDM_Bootstrap_', models[m], '_', pseudoabs_type), 
                         pattern = paste0(names[name_index], "_AUC_values.csv"),
                         full.names = TRUE)
    
    auc_val <- read.csv(aucval)
    auc_val$model_id <- paste0(names[name_index], '_', models[m])
    auc_out[[m]] <- auc_val
    
    
  }
  
  if(!length(mod_pred_out)){ # stop run if no models worked for a given species
    
    stop(paste("No models ran for species", names[name_index]))
    
  }
  
  # get things ready for combining
  means <- do.call('stack', mod_pred_out)
  qranges <- do.call('stack', qrange_out)
  
  # change auc into vector format for weighted average
  aucs <- do.call('c', lapply(auc_out, FUN=function(x) unique(x$meanAUC))) # removes NULL objects
  
  # check that number of entries match between raster and auc
  # they should already be in the correct order
  if(length(aucs) != nlayers(means)) {stop('Number of AUC weights does not match number of models')}
  
  
  ####   Store the output auc scores so we know which models were used - although should get that from the AUC table
  auc_out <- do.call(rbind, auc_out)[,3:4] %>% 
    unique() %>% 
    filter(meanAUC >= auc_cutoff)
  
  write.csv(auc_out, 
            file = paste0(output_directory, names[name_index], '_aucOuts.csv'))
  
  # check the auc values and drop any bad models
  if(any(aucs < auc_cutoff)){
    
    means <- dropLayer(means, which(aucs < 0.75)) # which() provides the index of which ones meet the statement
    qranges <- dropLayer(qranges, which(aucs < 0.75))
    
    aucs <- aucs[aucs>=auc_cutoff] # remove aucs < cutoff
    
  }
  
  
  ### 5. Create average model prediction 
  print("#####     calculating weighted ensemble raster     #####")
  wt_mean <- raster::weighted.mean(x = means, w = aucs)
  wt_qr <- raster::weighted.mean(x = qranges, w = aucs)
  
  
  ## save the outputs
  print("#####     Saving prediction raster     #####")
  writeRaster(x = wt_mean, 
              filename = paste0(output_directory, species_name, "_weightedmean.grd"),
              format = 'raster', overwrite = T)
  
  writeRaster(x = wt_qr, 
              filename = paste0(output_directory, species_name, "_weightedvariation.grd"),
              format = 'raster', overwrite = T)
  
  
  # store models that failed
  error_out <- do.call('rbind', errored_models)
  write.csv(error_out, 
            file = paste0(output_directory, names[name_index], '_failed_models.csv'))
  
}




#####     The script to create the job submission files
setwd(paste0('scripts/lotus/', taxa, '/combine_models_scripts'))

sdm_slurm <- slurm_apply(calculate_ensemble,
                         params = pars,
                         jobname = paste0(taxa, "_weighted_ensemble"),
                         nodes = length(unique(pars$name_index)),
                         cpus_per_node = 1,
                         slurm_options = list(partition = "short-serial",
                                              time = "23:59:59",
                                              mem = "30000",
                                              error = paste0('/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/outputs/', taxa, '/combined_model_outputs/error/combinmod_', pseudoabs_type,'-%j-%a.out')),
                         rscript_path = '',
                         submit = F)



