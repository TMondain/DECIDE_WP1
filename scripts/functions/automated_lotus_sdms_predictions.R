setwd("~/DECIDE/DECIDE_WP1") # to return to project directory

#####    Automated lotus script
library(rslurm)
library(tidyverse)

# taxa
taxa = 'moth' # moth, butterfly

# pseudoabsence name
pa_name = 'PA_thinned_10000nAbs'

# name of model to run
model = 'me'# one of c('lr', 'gam', 'rf', 'me')

# queue for lotus
queue_name = 'long-serial'

# time requirement
time = '47:59:59'

# memory requirement in MB
mem_req = 30000


# load data for parameters
if(taxa == 'moth'){
  
  ## for moths
  dfm_df <- read_csv("data/edited_insect_data/moth/DayFlyingMoths.csv")

} else if(taxa == 'butterfly'){
  
  ## for butterflies
  dfm_df <- read_csv("data/edited_insect_data/butterfly/BNM_NRMS_No_Duplicates_east_north.csv")
  
}

# get the parameters for function
pars <- data.frame(name_index = seq(1, length(names(dfm_df))))


## the function
slurm_sdm_boot <- function(name_index) {
  
  require(raster)
  require(dismo)
  require(tidyverse)
  require(mgcv)
  require(randomForest)
  require(glmnet)
  source("/home/users/thoval/DECIDE/scripts/Edited_Rob_Functions.R")
  
 
  
  
  #####     CHANGE FOR EACH RUN     #####
  
  # taxa
  taxa = 'moth' # moth, butterfly
  
  ## choose model of interest
  model = 'me' # one of c('lr', 'gam', 'rf', 'me')
  
  # pseudoabsence name
  pa_name = 'PA_thinned_10000nAbs'
 
  # number of bootstraps
  k = 10
  
  # number of knots for gam
  knots_gam = 4
  
  
  # where to save at the end
  outPath = paste0("/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/outputs/", taxa, "/SDM_Bootstrap_", model, "_", pa_name,"/")
  
  # for the species of interest 
  file_for_lotus = read.csv("/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/scripts/",taxa, "/_rslurm_", model, "_", pa_name, "/file_for_lotus.csv")
  
  # load environmental data
  env_dat <- raster::stack("/home/users/thoval/DECIDE/data/environmental_data/edat_nocorrs_nosea.gri")
  
  # load pseudoabsences based on taxa and name
  if(taxa == 'moth'){
    
    load(paste0("/home/users/thoval/DECIDE/data/species_data/moths/", pa_name, ".rdata"))
    
  } else if(taxa == 'butterfly'){
    
    load(paste0("/home/users/thoval/DECIDE/data/species_data/butterflies/", pa_name, ".rdata"))
    ab1 <- res_out
    
  }
  
  ## Find the species of interest
  species = file_for_lotus$species[name_index]

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
  if(model != 'lrReg') {
    
    boots_out <- raster::stack(lapply(sdm$Bootstrapped_models, FUN = function(x) predict(env_dat, x, type=type, index=index)))
    
    ## quantiles
    print(paste0('#####   getting quantiles   #####'))
    mean_preds <- calc(boots_out, fun = mean, na.rm = T)
    quant_preds <- calc(boots_out, fun = function(x) {quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)})
    rnge <- quant_preds[[2]]-quant_preds[[1]]
    
    
  } else if(model == 'lrReg') { 
    
    print(paste0('#####   predicting for lrReg bootstrapped models   #####'))
    
    ## convert variables to matrix
    covsMat <- as.matrix(rasterToPoints(env_dat)) # convert variables to matrix
    
    ## predict from lrReg model
    boots <- stack(lapply(sdm$Bootstrapped_models, FUN = function(x) {
      
      pred <- predict(x, covsMat[, 3:ncol(covsMat)], type = "response") # predict from matrix
      
      pred <- as.matrix(cbind(covsMat[, 1:2], pred)) # combine predictions with east - norths
      
      if (any(is.na(pred[, 3]))) pred <- pred[-which(is.na(pred[,3])), ] # get rid of NAs
      
      pred_rast <- rasterize(pred[, 1:2], env_dat[[1]], field = pred[, 3]) # turn into a raster of probabilities
      
      ## return predictions
      return(pred_rast)
      
    }))
    
    ## check for intercept-only models
    uniqueVals <-lapply(1:k,
                        function(x) { length(cellStats(boots[[x]], unique)) })
    
    drop <- which(uniqueVals <= 2) ## i.e. the mean and NA, because with intercept-only models the only values are 0.5 and NA
    
    # assign intercept only models an AUC of NULL - important for weighted average later
    if(any(drop)){
      
      sdm$AUC[drop] <- NA
      
      print(paste("Dropping", length(drop), "intercept-only model(s). Intercept-only models are given an AUC value of NA so they can be identified.
                        Where 1:(k-1) models are intercept only, only the non-intercept models are included in the final average. Where all models are intercept-only,
                        their predictions are returned but should not be used."))
      
    }
    
    
    if (length(drop) == k | length(drop) == 0) {
      
      ## quantiles
      print(paste0('#####   getting quantiles lrReg   #####'))
      mean_preds <- mean(boots) # where all models are intercept-only, takes the mean to avoid errors later but AUC scores are NA which means they are dropped for final ensembles
      quant_preds <- calc(boots, fun = function(x) {quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)})
      rnge <- quant_preds[[2]]-quant_preds[[1]]
      
      
    } else {
      
      print(paste0('#####   getting quantiles lrReg   #####'))
      mean_preds <- mean(boots[[-drop]])
      quant_preds <- calc(boots[[-drop]], fun = function(x) {quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)})
      rnge <- quant_preds[[2]]-quant_preds[[1]]
      
    }
    
    ## if some models were intercept-only then recalculate k (number of models used)
    k <- k - length(drop)
    
  }
  
  ## save files ##
  print("#####     Saving files     #####")
  species_name <- gsub(pattern = ' ', replacement = '_', species) # get species name without space
  
  # save prediction raster
  print("#####     Saving prediction raster     #####")
  writeRaster(x = mean_preds, 
              filename = paste0(outPath, model, "_SDMs_", species_name, "_meanpred.grd"),
              format = 'raster', overwrite = T)
  
  # save quantile max min
  print("#####     Saving quantile max min raster     #####")
  writeRaster(x = quant_preds, 
              filename = paste0(outPath, model, "_SDMs_", species_name, "_quantilemaxmin.grd"),
              format = 'raster', overwrite = T)
  
  # save quantile range raster
  print("#####     Saving quantile range raster     #####")
  writeRaster(x = rnge, 
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

# set working directory to correct taxa/location
setwd(paste0('scripts/lotus/', taxa, '/sdm_scripts/'))

#### slurm apply call
sdm_slurm <- slurm_apply(slurm_sdm_boot,
                         params = pars,
                         jobname = paste0(model, '_', pa_name),
                         nodes = length(unique(dfm_df$sp_n)),
                         cpus_per_node = 1,
                         slurm_options = list(partition = queue_name,
                                              time = as.character(time),
                                              mem = mem_req,
                                              error = paste0('/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/outputs/', taxa, '/log_files/SDM_Bootstrap_', model,'/SDM_Bootstrap_', model,'-%j-%a.out')),
                         rscript_path = '',
                         submit = F)

file_for_lotus = data.frame(species = unique(dfm_df$sp_n))
head(file_for_lotus)

write.csv(file_for_lotus, 
          file = paste0('_rslurm_', model, '_', pa_name, "/file_for_lotus.csv"))
            






