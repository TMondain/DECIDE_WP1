setwd("~/DECIDE/DECIDE_WP1") # to return to project directory
rm(list = ls())

#####    Automated lotus script
library(rslurm)
library(tidyverse)


#####     CHANGE FOR EACH RUN     #####

# taxa
taxa = 'butterfly' # moth, butterfly, nightflying_moth

# date records created 
recs_date = "2021_12_06"

## choose model of interest
model = c('lr', 'gam', 'rf', 'me')

# pseudoabsence name
pa_name = 'PA_thinned_10000nAbs'

# pseudoabsence date
pa_date_suffix = '2021_12_08'

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


# load data for parameters
if(taxa == 'moth'){
  
  ## for moths
  dfm_df <- read_csv(paste0("data/edited_insect_data/moth/DayFlyingMoths_EastNorths_no_duplicates_", recs_date, ".csv"))
  
} else if(taxa == 'butterfly'){
  
  ## for butterflies
  dfm_df <- read_csv(paste0("data/edited_insect_data/butterfly/butterfly_EastNorths_no_duplicates_", recs_date, ".csv"))
  
} else if(taxa == 'nightflying_moth'){
  
  dfm_df <- read_csv(paste0("data/edited_insect_data/moth/NightFlyingMoths_EastNorths_no_duplicates_", recs_date, ".csv"))
  
}

# get the parameters for function
pars <- data.frame(name_index = seq(1, length(unique(dfm_df$sp_n))*length(model)))


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
  file_for_lotus <- read.csv('file_for_lotus.csv')
  
  
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
  
  # date of pseudoabsence generation
  pa_date_suffix = file_for_lotus$pa_date_suffix[name_index]
  
  
  # where to save at the end
  outPath = paste0("/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/outputs/", taxa, "/SDM_Bootstrap_", model, "_", pa_name, "_", pa_date_suffix, "/")
  dir.create(outPath, recursive = TRUE)
  
  # # load file with parameters
  # file_for_lotus = read.csv(paste0("/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/scripts/", taxa, "/_rslurm_", model, "_", pa_name, "/file_for_lotus.csv"))
  
  # load environmental data
  env_dat <- raster::stack("/home/users/thoval/DECIDE/data/environmental_data/envdata_fixedcoasts_nocorrs_100m_GB.gri")
  
  
  # load pseudoabsences based on taxa and name
  load(paste0("/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/outputs/pseudoabsences/", taxa, '_', pa_name, "_", pa_date_suffix, ".rdata"))
  if(taxa == 'butterfly') {
    ab1 <- butt_out
  } else if(taxa == 'moth'){
    ab1 <- moth_out
  } else if(taxa == 'nightflying_moth') {
    ab1 <- nfmoth_out
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

dir.create(paste0('scripts/lotus/', taxa, '/sdm_scripts/'), recursive = TRUE)

# set working directory to correct taxa/location
setwd(paste0('scripts/lotus/', taxa, '/sdm_scripts/'))

#### slurm apply call
sdm_slurm <- slurm_apply(slurm_sdm_boot,
                         params = pars,
                         jobname = paste0(taxa, paste0(model, collapse = ''), '_', pa_name, '_', pa_date_suffix),
                         nodes = dim(pars)[1],
                         cpus_per_node = 1,
                         slurm_options = list(partition = queue_name,
                                              time = as.character(time),
                                              mem = mem_req,
                                              error = paste0(paste(model, collapse = '_'),'_error-%j-%a.err')),
                         rscript_path = '',
                         sh_template = "jasmin_submit_sh.txt",
                         submit = F)


file_for_lotus <- expand.grid(species = unique(dfm_df$sp_n),
                              taxa = taxa,
                              model = model,
                              pa_name = pa_name,
                              k = k,
                              pa_date_suffix = pa_date_suffix,
                              knots_gam = knots_gam)


# # old file for lotus for running each model separately - deprecated
# file_for_lotus = data.frame(species = unique(dfm_df$sp_n),
#                             taxa = rep(taxa, length = length(unique(dfm_df$sp_n))),
#                             model = rep(model, length = length(unique(dfm_df$sp_n))),
#                             pa_name = rep(pa_name, length = length(unique(dfm_df$sp_n))),
#                             k = rep(k, length = length(unique(dfm_df$sp_n))),
#                             knots_gam = rep(knots_gam, length = length(unique(dfm_df$sp_n))))


head(file_for_lotus)

write.csv(file_for_lotus, 
          file = paste0('_rslurm_', taxa, paste0(model, collapse = ''), '_', pa_name, '_', pa_date_suffix, "/file_for_lotus.csv"))







