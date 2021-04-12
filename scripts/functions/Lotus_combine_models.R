
###   Combine predictions from different models on lotus


#####     1. setup
taxa = 'moth'
main_directory = '/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/outputs'
output_directory = paste0(main_directory, '/',taxa, '/combined_model_outputs')
models = c('lr', 'lrReg', 'gam', 'rf', 'me')
pseudoabs_type = '10000nAbs'# which model run name to go through

auc_cutoff = 0.75 ## just a suggestion


# get species names
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
  
} else {print('whoooaaaahhhh there boy, du calme!')}


#####    2. start loop by going through all the names

# # sdm outputs for each species
# species_stack <- list()

# error outputs
error_out <- list()

# model outputs
mod_pred_out <- list()
qrange_out <- list()
auc_out <- list()

for(n in 1:length(names)){
  
  print(names[n])
  
  # initiate model list within for loop so that it gets replaced when starting a new species
  # otherwise we might get some weird overlaps
  model_stack <- list()
  errored_models <- list()
  
  # go through each model in turn
  for(m in 1:length(models)){
    
    ### first, check that the model for that species exists
    ### if it doesn't, skip to next species
    {check_models <- list.files(paste0(main_directory, '/', taxa, '/SDM_Bootstrap_', models[m], '_', pseudoabs_type), 
                                pattern = paste0(names[n]),
                                full.names = TRUE)
    
    if(length(check_models)<=1){
      
      print(paste('!!!   model', models[m], 'failed for species', names[n], '  !!!'))
      
      errored_models[[m]] <- data.frame(taxa = taxa, 
                                        species = names[n], 
                                        model = models[m])
      
      next
    }}
    
    
    ### get the mean predictions
    mp <- list.files(paste0(main_directory, '/', taxa, '/SDM_Bootstrap_', models[m], '_', pseudoabs_type), 
                     pattern = paste0(names[n], "_meanpred.grd"),
                     full.names = TRUE)
    
    mod_preds <- raster(mp)
    names(mod_preds) <- paste0(names[n], '_', models[m],'_mean_pred')
    mod_pred_out[[m]] <- mod_preds
    
    ### get the quantile range
    qr <- list.files(paste0(main_directory, '/', taxa, '/SDM_Bootstrap_', models[m], '_', pseudoabs_type), 
                     pattern = paste0(names[n], "_quantilerange.grd"),
                     full.names = TRUE)
    
    qrange <- raster(qr)
    names(qrange) <- paste0(names[n], '_', models[m], '_quantile_range')
    qrange_out[[m]] <- qrange 
    
    ### get the auc values
    aucval <- list.files(paste0(main_directory, '/', taxa, '/SDM_Bootstrap_', models[m], '_', pseudoabs_type), 
                         pattern = paste0(names[n], "_AUC_values.csv"),
                         full.names = TRUE)
    
    auc_val <- read.csv(aucval)
    auc_val$model_id <- paste0(names[n], '_', models[m])
    auc_out[[m]] <- auc_val
    
    
  }
  
  # store models that errored
  error_out[[n]] <- do.call('rbind', errored_models)
  
  # get things ready for combining
  means <- do.call('stack', mod_pred_out)
  qranges <- do.call('stack', qrange_out)
  
  # get auc into vector format
  aucs <- do.call('c',sapply(auc_out, FUN=function(x) unique(x$meanAUC))) # removes NULL objects  
  
  # check that number of entries match
  # they should already be in the correct order
  ####   Store the output so we know which models were used - although should get that from the AUC table
  
  
  # check the auc values and drop any bad models
  if(any(aucs < 0.75)){
    
    means <- dropLayer(means, which(aucs < 0.75))
    qranges <- dropLayer(qranges, which(aucs < 0.75))
    
    
  }
  
  #### create average model prediction 
  wt_mean <- raster::weighted.mean(x=means, w = aucs)
  wt_qr <- raster::weighted.mean(x=qranges, w = aucs)
  
  ## save the outputs
  
}



error_df <- do.call('rbind', error_out)




