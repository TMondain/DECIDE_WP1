setwd("~/DECIDE/DECIDE_WP1")
rm(list = ls())

## Lotus create decide score
library(rslurm)

# going to run it as a slurm_apply call to run butterflies and moths

# pseudoabsence name (i.e. the model)
pseudoabs_name = 'PA_thinned_10000nAbs'

# date for models
pa_date_suffix = "2021_12_08"

# date for records
recs_date_suffix = '2021_12_06'

# Combine SD method (to be  tagged onto pseudo abs name)
sd_comb = 'mean_score' # or 'mean_score' or ''

# decide score method
method = 'var_only' #c('sqroot_var_preds', 'equal_weighting', 'sqroot_preds_var', 'var_only') # 'var_sqroot_preds', 'equal_weighting', 'preds_sqroot_var', 'var_only'

# combine decide score method
comb_method = 'mean' # c('weight_mean', 'mean')

# main directory
main_dir = '/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/outputs/'

# how to weight decide score by records
weight_record = 'time_since' # NULL if no weighting by records is desired

# taxa
taxa = 'nightflying_moth' # c('moth', 'butterfly', 'nightflying_moth')

# parameter file
pars = data.frame(index = 1:(length(method) * length(comb_method) * length(taxa)))

decide_score <- function(index){
  
  require(raster)
  require(tidyverse)
  require(lubridate)
  require(sp)
  
  # read in file
  file_for_lotus <- read.csv('file_for_lotus.csv')
  
  # route directory
  main_dir = file_for_lotus$main_dir[index]
  
  # get the parameters set
  taxa = file_for_lotus$taxa[index]
  pseudoabs_name = file_for_lotus$pseudoabs_name[index]
  sd_comb = file_for_lotus$sd_comb[index]
  method = file_for_lotus$method[index]
  comb_method = file_for_lotus$comb_method[index]
  weight_record = file_for_lotus$weight_record[index]
  recs_date_suffix = file_for_lotus$recs_date_suffix[index]
  pa_date_suffix = file_for_lotus$pa_date_suffix[index]
  
  ## load in species
  # get names
  fls <- paste0(main_dir, taxa,'/combined_model_outputs/', pseudoabs_name, '_', pa_date_suffix)
  
  print('!! stacking rasters')
  
  # predictions
  preds <- raster::stack(list.files(fls, 
                                    pattern = paste0(sd_comb, '_weightedmeanensemble.grd'), 
                                    full.names = T))
  
  # variation
  var <- raster::stack(list.files(fls, 
                                  pattern = paste0(sd_comb, '_weightedvariationensemble.grd'), 
                                  full.names = T))
  
  print('!! calculating decide score')
  
  ## combine them together to create decide score
  if(method == 'sqroot_var_preds'){
    dec_spp <- preds*sqrt(var)
  } else if(method =='equal_weighting'){
    dec_spp <- preds*var
  } else if(method == 'sqroot_preds_var'){
    dec_spp <- sqrt(preds)*var
  } else if(method == 'var_only'){
    dec_spp <- var
  } else { stop(paste('!!!   script only coded for methods:
                     var_sqroot_preds, equal_weighting, preds_sqroot_var  !!!')) }
  
  # get mean across all species
  if(comb_method == 'mean'){
    decide <- calc(dec_spp, mean)
  } else if(comb_method == 'weight_mean'){
    decide <- raster::weighted.mean(dec_spp, preds)
  } else { stop(paste('!!!   script only coded for combining methods mean and weight_mean   !!!')) }
  
  
  ## weight the DECIDE score by records
  if(!is.null(weight_record)){
    
    print('!! loading records')
    
    # load records
    if(taxa == 'moth'){
      
      # this being hard-coded is bad!!!
      dfm_full <- read.csv(paste0("/home/users/thoval/DECIDE/data/species_data/DayFlyingMoths_EastNorths_no_duplicates_", recs_date_suffix, ".csv"))
      
    } else if(taxa == 'butterfly'){
      
      # this being hard-coded is bad!!!
      dfm_full <- read.csv(paste0("/home/users/thoval/DECIDE/data/species_data/butterfly_EastNorths_with_duplicates_", recs_date_suffix, ".csv"))
      
    } else if(taxa == 'nightflying_moth') {
      
      dfm_full <- read.csv(paste0("/home/users/thoval/DECIDE/data/species_data/NightFlyingMoths_EastNorths_no_duplicates_", recs_date_suffix, ".csv"))
      
    } else { stop('!! Only coded for butterflies and moths !!') }
    
    
    # simplify dataframe
    dfm <- dfm_full %>% 
      dplyr::select(lon, lat, date, year = yr, species = sp_n, common_name=com_n)
    
    # # convert dataframe to sf object
    # st_dfm <- st_as_sf(dfm, coords = c('lon', 'lat'), crs = 27700)
    
    
    ## weight the records
    
    ### Function to create the weighting layer - stupid doing this for each lotus run - needless extra time
    count_records <- function(records_df, # sf object or data frame/data table of record counts
                              template_raster, # raster that you want the counts to be aggregated by
                              Coords = 27700, 
                              weight_by_time = TRUE)
    {
      
      # Get the counts of records per cell and store as data frame
      if(class(records_df)[1] == 'sf'){
        
        record_counts <- records_df %>% 
          mutate(lon = unlist(map(records_df$geometry, 1)),
                 lat = unlist(map(records_df$geometry, 2)))  %>% 
          as.data.frame() %>% 
          mutate(geometry = NULL)
        
      } else if(class(records_df)[1] == 'data.table'|class(records_df)[1] == 'data.frame'){
        
        record_counts <- records_df
        
      } else { stop('!! records_df must be class sf or data.frame') }
      
      if(weight_by_time){
        
        print('!! summarising records for weight by time') 
        
        record_weights <- record_counts %>% 
          group_by(lon, lat) %>% 
          dplyr::summarise(last_rec = max(year), # what is the most recent year?
                           last_date = ymd(max(date)), # what is the most recent date?
                           yrs_since_last_rec = (year(Sys.Date())-last_rec), # number of years since the last record in a grid cell
                           days_since_last_rec = as.numeric(Sys.Date()) - as.numeric(last_date), # number of days since the last record in a grid cell
                           score = 1/yrs_since_last_rec) %>% # get the big values small and vice-versa, so that larger numbers are 'bad' and smaller numbers are 'good', so that it matches the number of records layer that's also outputted by this function
          ungroup()
        
      } else if(!weight_by_time){
        
        record_weights <- record_counts
        
      } else {stop('Weight by time mst be TRUE/FALSE')}
      
      # convert to a spatial points data frame to match with raster
      xy <- record_weights[,c("lon","lat")]
      spdf.recs <- SpatialPointsDataFrame(coords = xy, data = record_weights,
                                          proj4string = CRS(paste0("+init=epsg:", Coords))) 
      
      ### create counts raster to store number of records in each cell ###
      # this is used to make an 'effort' layer
      n_recs <- template_raster
      
      # make a raster of zeroes for input
      n_recs[!is.na(n_recs)] <- 0
      
      if(weight_by_time){
        
        # get cell numbers for the spdf.recs data frame
        cells  <- (cellFromXY(n_recs,spdf.recs))
        
        # fill those cells with the score, have done some tests and 98% sure they are in the same order but ideally would need a match statement here
        n_recs[cells] <- spdf.recs$score
        
      } else if(!weight_by_time){
        
        # get the cell index for each point and make a table:
        counts  <- table(cellFromXY(n_recs,spdf.recs))
        
        # fill in the raster with the counts from the cell index:
        n_recs[as.numeric(names(counts))] <- counts
        
      }
      
      # mask cells falling outside of raster values
      n_recs <- mask(n_recs, template_raster)
      
      return(n_recs)
      
    }
    
    
    # create the weighting layer
    # weight by time since last record or raw number of records (sightings unique to each day)
    if(weight_record == 'time_since'){
      
      print('!! counting records using time since method')
      
      rec_counts <- count_records(records_df = dfm, 
                                  template_raster = decide,
                                  weight_by_time = TRUE)
      
    } else if(weight_record == 'number_record'){
      
      print('!! counting records total numbers')
      
      rec_counts <- count_records(records_df = dfm, 
                                  template_raster = decide,
                                  weight_by_time = FALSE)
      
    } else { warning('!! Currently only coded for time since and number_records') }
    
    
    ## smooth the DECIDE score
    
    ### Function to weight the DECIDE score
    
    smooth_recording <- function(weighted_layer, # layer to be weighted
                                 effort_raster,
                                 sm_matrix = matrix(c(0,    0.09, 0.11, 0.09,    0,
                                                      0.09, 0.21, 0.33, 0.21, 0.09,
                                                      0.11, 0.33,    1, 0.33, 0.11,
                                                      0.09, 0.21, 0.33, 0.21, 0.09,
                                                      0,    0.09, 0.11, 0.09,    0), 
                                                    nrow = 5, ncol = 5),
                                 fun = sum,
                                 recording_impact = 10) # the larger the value the less impact it has
    {
      
      smoothed_effort <- focal(x = effort_raster, 
                               w = sm_matrix,
                               fun = fun,
                               pad = TRUE,
                               na.rm = TRUE,
                               padValue = NA,
                               NAonly = FALSE)
      
      # Convert recording to weighting
      weighting <- smoothed_effort
      weighting <- 1/(1+(weighting*recording_impact)) # Outputs a layer of 1s where there are no records and progressively goes to 0 where there are more records
      # have added a 1 to stop the possibility of dividing by 0
      
      # set minimum value to original layer
      M <- minValue(weighted_layer)
      
      # for now, just multiplying the weighted_layer by the weighting
      adjusted_score <- weighted_layer * weighting
      adjusted_score[adjusted_score < M] <- M # change the lowest records to be equal to the original minimum value in the area 
      
      
      return(list(weighting_layer = weighting,
                  weighted_score = adjusted_score))
      
    }
    
    print('!! smoothing the effort layer')
    
    # do the smoothing
    sr <- smooth_recording(weighted_layer = decide,
                           effort_raster = rec_counts,
                           recording_impact = 5)
    
    decide <- sr$weighted_score
    
  }
  
  print('!! saving')
  
  ## write single raster 
  writeRaster(x = decide, 
              filename = paste0(main_dir, 'decide_scores/outputs/', taxa, '_',  pseudoabs_name, '_', pa_date_suffix, '_', sd_comb, '_decide_score_', method, '_', comb_method, '_', weight_record, '.grd'),
              format = 'raster', overwrite = T)
  
}


#####     The script to create the job submission files
setwd(paste0('scripts/lotus/decide_scores/'))

sdm_slurm <- slurm_apply(decide_score,
                         params = pars,
                         jobname = paste0("decide_scores_", paste(taxa, collapse = '_'), '_', pseudoabs_name, '_', pa_date_suffix, '_', sd_comb, "_", paste(method, collapse = '_'), "_", comb_method, '_', weight_record),
                         nodes = length(unique(pars$index)),
                         cpus_per_node = 1,
                         slurm_options = list(partition = "long-serial",
                                              time = "167:59:59",
                                              mem = "20000",
                                              output = paste0(paste(method, collapse = '_'),'_out-%j-%a.out'),
                                              error = paste0(paste(method, collapse = '_'),'_error-%j-%a.err')),
                         rscript_path = '',
                         sh_template = "jasmin_submit_sh.txt",
                         submit = F)


# create file for lotus
file_for_lotus <- expand.grid(taxa = taxa,
                              pseudoabs_name = pseudoabs_name,
                              pa_date_suffix = pa_date_suffix,
                              recs_date_suffix = recs_date_suffix,
                              sd_comb = sd_comb,
                              method = method, 
                              main_dir = main_dir,
                              comb_method = comb_method, 
                              weight_record = weight_record)

write.csv(file_for_lotus, 
          file = paste0('_rslurm_decide_scores_', paste(taxa, collapse = '_'), '_', pseudoabs_name, '_', pa_date_suffix, '_', sd_comb, '_', paste(method, collapse = '_'), '_', comb_method, '_', weight_record, '/file_for_lotus.csv'))
