

## Lotus create decide score
library(rslurm)

# going to run it as a slurm_apply call to run butterflies and moths

# parameter files
pars = data.frame(index = 1:2)

# pseudoabsence name (i.e. the model)
pseudoabs_name = 'PA_thinned_10000nAbs'

file_for_lotus <- data.frame(taxa = c('moth', 'butterfly'),
                             pseudoabs_name = rep('PA_thinned_10000nAbs', 2))

decide_score <- function(index){
  
  # read in file
  file_for_lotus <- read.csv('file_for_lotus.csv')
  
  # get the parameters set
  taxa = file_for_lotus$taxa[index]
  preudoabs_name = file_for_lotus$preudoabs_name[index]
  
  ## load in species
  # get names
  fls <- paste0('/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/outputs/', taxa,'/combined_model_outputs/',pseudoabs_name)
  
  # predictions
  preds <- raster::stack(list.files(fls, 
                                    pattern = '_weightedmeanensemble.grd', 
                                    full.names = T))
  
  # variation
  var <- raster::stack(list.files(fls, 
                                  pattern = '_rangeensemblequantiles.grd', 
                                  full.names = T))
  
  
  ## combine them together to create decide score
  dec_spp <- preds*sqrt(var)
  
  # get mean across all species
  decide <- calc(dec_spp, mean)
  
  ## write single raster 
  
  
  
}