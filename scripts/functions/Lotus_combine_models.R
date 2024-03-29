## project directory
setwd("~/DECIDE/DECIDE_WP1")
rm(list = ls())

###   Combine predictions from different models on lotus
library(tidyverse)
library(raster)
library(rslurm)

#####  Get parameters

# taxa for slurm output and parameter loading
taxa = 'moth' # moth, butterfly, nightflying_moth
pseudoabs_type = 'PA_thinned_10000nAbs'
auc_cutoff = 0.75 ## just a suggestion - might need some thought (although AUC values so stupidly high might not be a problem until we're using a different score metric)
sd_comb='mean_score'

# pseudoabsence date
pa_date_suffix = '2021_12_08'

# models = c('lr', 'gam', 'rf', 'me') #, 'lrReg') ## lrReg hasn't worked for any species yet

# use the species names from the pseudoabsences
load(paste0("outputs/pseudoabsences/", taxa, "_", pseudoabs_type, "_", pa_date_suffix, ".rdata"))

if(taxa == 'butterfly') {
  pars <- data.frame(name_index = seq(1, length(names(butt_out))))
  species <- names(butt_out)
} else if(taxa == 'moth'){
  pars <- data.frame(name_index = seq(1, length(names(moth_out))))
  species <- names(moth_out)
} else if(taxa == 'nightflying_moth') {
  pars <- data.frame(name_index = seq(1, length(names(nfmoth_out))))
  species <- names(nfmoth_out)
}

###   function
calculate_ensemble <- function(name_index) {
  
  require(tidyverse)
  require(raster)
  
  # read the file_for_lotus.csv
  file_for_lotus <- read.csv('file_for_lotus.csv')
  
  ### 1. parameters for function
  # species
  names = file_for_lotus$species[name_index]
  
  # taxa
  taxa = file_for_lotus$taxa[name_index] 
  
  # pseudoabsence type
  pseudoabs_type = file_for_lotus$pseudoabs_type[name_index] # which model run name to go through
  
  # auc cutoff for dropping models
  auc_cutoff = file_for_lotus$auc_cutoff[name_index] ## just a suggestion - might need some thought (although AUC values so stupidly high might not be a problem until we're using a different score metric)
  
  # how to combine standard deviation
  sd_comb=file_for_lotus$sd_comb[name_index]
  
  # models to check
  models = c('lr', 'gam', 'rf', 'me') #, 'lrReg') ## lrReg hasn't worked for any species yet
  
  # pseudoabsence date
  pa_date_suffix = file_for_lotus$pa_date_suffix[name_index]
  
  ## no need to change these unless changing directories
  main_directory = '/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/outputs'
  output_directory = paste0(main_directory, '/',taxa, '/combined_model_outputs/', pseudoabs_type, '_', pa_date_suffix, '/')
  dir.create(output_directory)
  
  
  ### 2. setup storage for model loops
  # error outputs
  error_out <- list()
  
  # model outputs
  mod_pred_out <- list()
  qrange_out <- list()
  qminmax_out <- list()
  auc_out <- list()
  errored_models <- list()
  
  
  print(names)
  
  
  ### 4. Go through each model in turn
  for(m in 1:length(models)){
    
    ### first, check that the model for that species exists
    ### if it doesn't, skip to next model
    {check_models <- list.files(paste0(main_directory, '/', taxa, '/SDM_Bootstrap_', models[m], '_', pseudoabs_type, '_', pa_date_suffix), 
                                pattern = paste0(names),
                                full.names = TRUE)
    
    if(length(check_models)<=1){
      
      print(paste('!!!   model', models[m], 'failed for species', names, '  !!!'))
      
      errored_models[[m]] <- data.frame(taxa = taxa, 
                                        species = names, 
                                        model = models[m])
      
      next
    }}
    
    
    ### load the mean predictions
    mp <- list.files(paste0(main_directory, '/', taxa, '/SDM_Bootstrap_', models[m], '_', pseudoabs_type, '_', pa_date_suffix), 
                     pattern = paste0(names, "_meanpred.grd"),
                     full.names = TRUE)
    
    mod_preds <- raster(mp)
    names(mod_preds) <- paste0(names, '_', models[m],'_mean_pred')
    mod_pred_out[[m]] <- mod_preds
    
    ### load the quantile range
    qr <- list.files(paste0(main_directory, '/', taxa, '/SDM_Bootstrap_', models[m], '_', pseudoabs_type, '_', pa_date_suffix), 
                     pattern = paste0(names, "_quantilerange.grd"),
                     full.names = TRUE)
    
    qrange <- raster(qr)
    names(qrange) <- paste0(names, '_', models[m], '_quantile_range')
    qrange_out[[m]] <- qrange
    
    # ### load the quantile min max
    # # removed because moving to standard deviation
    # qmm <- list.files(paste0(main_directory, '/', taxa, '/SDM_Bootstrap_', models[m], '_', pseudoabs_type), 
    #                   pattern = paste0(names, "_quantilemaxmin.grd"),
    #                   full.names = TRUE)
    # 
    # qminmax <- raster::stack(qmm)
    # names(qminmax) <- c(paste0(names, '_', models[m], '_quantile_min'),
    #                     paste0(names, '_', models[m], '_quantile_max'))
    # qminmax_out[[m]] <- qminmax
    
    
    ### load the auc values
    aucval <- list.files(paste0(main_directory, '/', taxa, '/SDM_Bootstrap_', models[m], '_', pseudoabs_type, '_', pa_date_suffix), 
                         pattern = paste0(names, "_AUC_values.csv"),
                         full.names = TRUE)
    
    auc_val <- read.csv(aucval)
    auc_val$model_id <- paste0(names, '_', models[m])
    auc_out[[m]] <- auc_val
    
    
  }
  
  if(!length(mod_pred_out)){ # stop run if no models worked for a given species
    
    stop(paste("No models ran for species", names))
    
  }
  
  ###   Combine means of models that ran
  
  if(any(sapply(mod_pred_out, is.null))){
    
    # means and quantile range for each species
    means <- do.call('stack', mod_pred_out[-which(sapply(mod_pred_out, is.null))])
    qranges <- do.call('stack', qrange_out[-which(sapply(qrange_out, is.null))])
    
    # # get max and min of corresponding raster layers
    # min_stack <- stack(lapply(qminmax_out[-which(sapply(qminmax_out, is.null))], FUN = function(x) x[[1]]))
    # max_stack <- stack(lapply(qminmax_out[-which(sapply(qminmax_out, is.null))], FUN = function(x) x[[2]]))
    
  } else if(!any(sapply(mod_pred_out, is.null))){
    
    # means and quantile range for each species
    means <- do.call('stack', mod_pred_out)
    qranges <- do.call('stack', qrange_out)
    
    # # get max and min of corresponding raster layers
    # min_stack <- stack(lapply(qminmax_out, FUN = function(x) x[[1]]))
    # max_stack <- stack(lapply(qminmax_out, FUN = function(x) x[[2]]))
    
  }
  
  
  # change auc into vector format for weighted average
  aucs <- do.call('c', lapply(auc_out, FUN=function(x) unique(x$meanAUC))) # removes NULL objects
  
  # check that number of entries match between raster and auc
  # they should already be in the correct order
  if(length(aucs) != nlayers(means)) {stop('Number of AUC weights does not match number of models')}
  
  
  ### check the auc values and drop any bad models
  # if all models are below auc_cutoff then...
  if(all(aucs < auc_cutoff)) {
    
    ####   Store the output auc scores so we know which models were used - although should get that from the AUC table
    auc_df <- do.call(rbind, auc_out)[,3:4] %>% 
      unique() %>% 
      mutate(used = 'used_bad_models')
    
    # need to create auc_weights function for the next if() statement to work when all models < cutoff
    auc_weights <- aucs
    
    write.csv(auc_df, 
              file = paste0(output_directory, names, "_", sd_comb, '_aucOuts.csv'))
    
  } else if(any(aucs < auc_cutoff) | all(aucs > auc_cutoff)) { # this else if statement works with all cases where at least some of the models have AUC > than the cutoff.
    
    means <- dropLayer(means, which(aucs < auc_cutoff)) # which() provides the index of which ones meet the statement
    qranges <- dropLayer(qranges, which(aucs < auc_cutoff))
    # min_stack <- dropLayer(min_stack, which(aucs < auc_cutoff))
    # max_stack <- dropLayer(max_stack, which(aucs < auc_cutoff))
    
    # remove the auc values that fall below cutoff to make sure weighted.mean functions work
    auc_weights <- aucs[aucs>=auc_cutoff] # remove aucs < cutoff
    
    ####   Store the output auc scores so we know which models were used - although should get that from the AUC table
    auc_df <- do.call(rbind, auc_out)[,3:4] %>% 
      unique() %>% 
      mutate(used = ifelse(meanAUC >= auc_cutoff, 'used', 'dropped'))
    
    if(any(auc_df$used == 'dropped')) {print(paste("!!  Some models were dropped because they were lower than", 
                                                   auc_cutoff, ". Models dropped:", auc_df$model_id[auc_df$used == 'dropped']))}
    
    write.csv(auc_df, 
              file = paste0(output_directory, names, "_", sd_comb, '_aucOuts.csv'))
    
  } 
  
  
  ### 5. Create average model prediction 
  if(length(auc_weights) > 1){ # if there is more than one model then calculate weighted average
    
    print("#####     calculating weighted ensemble raster     #####")
    wt_mean <- raster::weighted.mean(x = means, w = auc_weights)
    
    # whether to take the mean standard deviation (weighted by auc)
    # or take the maximum standard deviation
    if(sd_comb == 'mean_score'){
      wt_qr <- raster::weighted.mean(x = qranges, w = auc_weights)
    } else if(sd_comb =='max_score'){
      wt_qr <- calc(qranges, fun = max)
    }
    
  } else { # otherwise just save the means as the wt_mean and ranges
    
    print("#####     Storing single model as a raster     #####")
    wt_mean <- means
    wt_qr <- qranges
    
  }
  
  # # Find the maximum and minimum quantiles across all models that ran
  # if(nlayers(min_stack)==1){
  #   
  #   # if only one model run then just store the min and max stacks
  #   # as the quantile output - makes code run faster.
  #   min_quant <- min_stack
  #   max_quant <- max_stack
  #   
  # } else if(nlayers(min_stack>1)){
  #   
  #   min_quant <- calc(min_stack, min)
  #   max_quant <- calc(max_stack, max)
  #   
  # }
  # 
  # minmax_stack <- stack(min_quant, max_quant)
  # names(minmax_stack) <- c('min_quant', 'max_quant')
  # mod_quant_rnge <- minmax_stack[[2]] - minmax_stack[[1]]
  
  
  ## save the outputs
  print("#####     Saving prediction raster     #####")
  writeRaster(x = wt_mean, 
              filename = paste0(output_directory, names, "_", pseudoabs_type, '_', pa_date_suffix, "_", sd_comb, "_weightedmeanensemble.grd"),
              format = 'raster', overwrite = T)
  
  writeRaster(x = wt_qr, 
              filename = paste0(output_directory, names, "_", pseudoabs_type, '_', pa_date_suffix, "_", sd_comb, "_weightedvariationensemble.grd"),
              format = 'raster', overwrite = T)
  
  # writeRaster(x = mod_quant_rnge, 
  #             filename = paste0(output_directory, names, "_", pseudoabs_type, "_rangeensemblequantiles.grd"),
  #             format = 'raster', overwrite = T)
  
  
  # store models that failed
  error_out <- do.call('rbind', errored_models)
  write.csv(error_out, 
            file = paste0(output_directory, names, "_", pseudoabs_type, '_', pa_date_suffix, "_", sd_comb, '_failed_models.csv'))
  
}




#####     The script to create the job submission files
setwd(paste0('scripts/lotus/', taxa, '/combine_models_scripts'))

sdm_slurm <- slurm_apply(calculate_ensemble,
                         params = pars,
                         jobname = paste0(taxa, "_", pseudoabs_type, "_", pa_date_suffix, "_", sd_comb, "_weighted_ensemble"),
                         nodes = length(unique(pars$name_index)),
                         cpus_per_node = 1,
                         slurm_options = list(partition = "short-serial",
                                              time = "02:00:00",
                                              mem = "30000",
                                              error = paste0('error_combinmod_', pseudoabs_type, "_", sd_comb,'-%j-%a.out')),
                         rscript_path = '',
                         sh_template = "../sdm_scripts/jasmin_submit_sh.txt",
                         submit = F)


## create the file for lotus
file_for_lotus <- data.frame(species = gsub(pattern = ' ', replacement = '_', x = unique(species)),
                             taxa = rep(taxa, length = length(unique(species))),
                             pseudoabs_type = rep(pseudoabs_type, length = length(unique(species))),
                             pa_date_suffix = pa_date_suffix,
                             auc_cutoff = rep(auc_cutoff, length = length(unique(species))),
                             sd_comb = rep(sd_comb, length = length(unique(species))))
head(file_for_lotus)

write.csv(file_for_lotus, 
          file = paste0('_rslurm_', taxa, "_", pseudoabs_type, "_", pa_date_suffix, "_", sd_comb, "_weighted_ensemble", "/file_for_lotus.csv"))
