

### workflow for lotus models
getwd()

## function for lotus
slurm_sdm_boot <- function(index) {
  
  require(raster)
  require(dismo)
  require(tidyverse)
  require(mgcv)
  require(randomForest)
  source("/home/users/thoval/DECIDE/scripts/Edited_Rob_Functions.R")
  
  # load environmental data
  env_dat <- raster::stack("/home/users/thoval/DECIDE/data/environmental_data/edat_nocorrs.gri")
  
  
  
  #####     CHANGE FOR EACH RUN     #####
  
  # load pseudoabsences
  load("/home/users/thoval/DECIDE/data/species_data/moths/moth_PA_unthinned_matched.rdata")
  
  # number of bootstraps
  k = 10
  
  # number of knots for gam
  knots = 4
  
  # where to save at the end
  outPath = "/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/outputs/SDM_Bootstrap_lr/"
  
  # for the species of interest
  file_for_lotus = read.csv("/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/scripts/_rslurm_SDM_Bootstrap_lr/file_for_lotus.csv")
  
  ## Find the species of interest
  species = file_for_lotus$species[index]
  
  ## choose model of interest
  model = 'lr'# one of c('lr', 'gam', 'rf', 'me')
  
  # run the model of interest
  sdm <- fsdm(species = species, model = model, 
              climDat = env_dat, spData = ab1, knots = knots,
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
  boots <- lapply(sdm$Bootstrapped_models, FUN = function(x) predict(env_dat, x, type=type, index=index))
  boot_out <- do.call("stack", boots)
  
  ## quantiles
  print(paste0('#####   getting quantiles   #####'))
  mean_preds <- calc(boot_out, fun = mean, na.rm = T)
  quant_preds <- calc(boot_out, fun = function(x) {quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)})
  rnge <- quant_preds[[2]]-quant_preds[[1]]
  
  model_output <- list(model = model,
                       sdm_output = sdm,
                       mean_predictions = mean_preds, # mean probability of presence across all bootstraps
                       quantile_predictions = quant_preds, # the lower and upper bound quantiles
                       quantile_range = rnge) # the upper minus the lower quantile
  
  ## save file
  save(model_output, file = paste0(outPath, "SDMs_", species, 
                                   ".rdata"))
  
  return(model_output)
  
}

