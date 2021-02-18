

### workflow for lotus models



slurm_sdm_boot <- function(index) {
  # presAbs, envDat, models, species, k, knots, outPath
  
  require(raster)
  require(dismo)
  require(tidyverse)
  require(mgcv)
  require(randomForest)
  source("/home/users/thoval/DECIDE/scripts/Edited_Rob_Functions.R")
  
  load("/home/users/thoval/DECIDE/data/species_data/moths/moth_PA_unthinned_matched.rdata")
  env_dat <- raster::stack("/home/users/thoval/DECIDE/data/environmental_data/edat_nocorrs.gri")
  k=10
  knots=4
  outPath = "/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/outputs/SDM_Bootstrap_lr/"
  file_for_lotus = read.csv("/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/scripts/_rslurm_SDM_Bootstrap_lr/file_for_lotus.csv")
  
  species = file_for_lotus$species[index]
  
  # commented out the other models because it took too long to include
  # all models in the same job
  models = c('lr')#, 'gam', 'rf', 'me')
  
  all_mods <- list()
  
  for(i in models){
    
    # run model i
    sdm <- fsdm(species = species, model = i, 
                climDat = env_dat, spData = ab1, knots = knots,
                k = k, # number of bootstraps
                prediction = T, # with prediction
                inters = F, 
                write = F, outPath = outPath)
    
    
    # run standard error for lr and gam
    if(i == 'lr'| i ==  'gam'){
      
      se_out <- predict(env_dat, sdm$Model, type = "response", se.fit = TRUE, index = 2)
      
    } else {
      se_out <- NULL 
    }
    
    # choose the type and index for predict function
    if(i == 'lr'|i == 'gam'){
      type <- "response"
      index <- NULL
    } else if(i == 'rf'){
      type <- "prob"
      index <- 2
    }
    
    ## bootstrapped models
    print(paste0('#####   predicting from bootstrapped models   #####')) 
    
    ## predict from each of the bootstrapped models
    boots <- lapply(sdm$Bootstrapped_models, FUN = function(x) predict(env_dat, x, type=type, index=index))
    boot_out <- do.call("stack", boots)
    
    ## quantiles
    print(paste0('#####   getting quantiles   #####'))
    mean_preds <- calc(boot_out, fun = mean, na.rm = T)
    quant_preds <- calc(boot_out, fun = function(x) {quantile(x, probs = c(0.05, 0.95), na.rm = T)})
    rnge <- quant_preds[[2]]-quant_preds[[1]]
    
    all_mods[[i]] <- list(model = i,
                          sdm_output = sdm,
                          se_prediction = se_out,
                          mean_predictions = mean_preds, # mean across all bootstraps
                          quantile_predictions = quant_preds, # the lower and upper bound quantiles
                          quantile_range = rnge) # the upper minus the lower quantile
    
  }
  
  save(all_mods, file = paste0(outPath, "SDMs_", species, 
                               ".rdata"))
  
  return(all_mods)
  
}

