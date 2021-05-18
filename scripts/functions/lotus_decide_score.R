setwd("~/DECIDE/DECIDE_WP1")

## Lotus create decide score
library(rslurm)

# going to run it as a slurm_apply call to run butterflies and moths

# parameter files
pars = data.frame(index = 1:2)

# pseudoabsence name (i.e. the model)
pseudoabs_name = 'PA_thinned_10000nAbs'

# decide score method
method = 'var_sqroute_preds'

# combine decide score method
comb_method = 'weight_mean'

# main directory
main_dir = '/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/outputs/'

decide_score <- function(index){
  
  require(raster)
  
  # read in file
  file_for_lotus <- read.csv('file_for_lotus.csv')
  
  # route directory
  main_dir = file_for_lotus$main_dir[index]
  
  # get the parameters set
  taxa = file_for_lotus$taxa[index]
  pseudoabs_name = file_for_lotus$pseudoabs_name[index]
  method = file_for_lotus$method[index]
  comb_method = file_for_lotus$comb_method[index]
  
  ## load in species
  # get names
  fls <- paste0(main_dir, taxa,'/combined_model_outputs/',pseudoabs_name)
  
  # predictions
  preds <- raster::stack(list.files(fls, 
                                    pattern = '_weightedmeanensemble.grd', 
                                    full.names = T))
  
  # variation
  var <- raster::stack(list.files(fls, 
                                  pattern = '_rangeensemblequantiles.grd', 
                                  full.names = T))
  
  
  ## combine them together to create decide score
  if(method == 'var_sqroute_preds'){
    dec_spp <- preds*sqrt(var)
  } else { stop(paste('!!!   script only coded for methods var_sqroute_preds   !!!')) }
  
  # get mean across all species
  if(comb_method == 'mean'){
    decide <- calc(dec_spp, mean)
  } else if(comb_method == 'weight_mean'){
    decide <- raster::weighted.mean(dec_spp, preds)
  } else { stop(paste('!!!   script only coded for combining methods mean and weight_mean   !!!')) }
  
  
  ## write single raster 
  writeRaster(x = decide, 
              filename = paste0(main_dir, 'decide_scores/outputs/', taxa, '_',  pseudoabs_name, '_decide_score_', method, '_', comb_method, '.grd'),
              format = 'raster', overwrite = T)
  
}


#####     The script to create the job submission files
setwd(paste0('scripts/lotus/decide_scores/'))

sdm_slurm <- slurm_apply(decide_score,
                         params = pars,
                         jobname = paste0("decide_scores_", pseudoabs_name, "_", method, "_", comb_method),
                         nodes = length(unique(pars$index)),
                         cpus_per_node = 1,
                         slurm_options = list(partition = "short-serial",
                                              time = "23:59:59",
                                              mem = "30000"),
                                              # error = paste0(main_dir,'decide_scores/', '_rslurm_decide_scores_', pseudoabs_name, '_', method, '/', method,'_error-%j-%a.out')),
                         rscript_path = '',
                         submit = F)


# create file for lotus
file_for_lotus <- data.frame(taxa = c('moth', 'butterfly'),
                             pseudoabs_name = rep(pseudoabs_name, 2),
                             method = rep(method, 2), 
                             main_dir = rep(main_dir, 2),
                             comb_method = rep(comb_method, 2))

write.csv(file_for_lotus, 
          file = paste0('_rslurm_decide_scores_', pseudoabs_name, '_', method,'_', comb_method, '/file_for_lotus.csv'))
