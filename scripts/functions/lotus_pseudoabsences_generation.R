# Create pseudoabsences for butterflies and moths

slurm_pseudoabs_fun <- function(species_name, nAbs = 10000, matchPres = FALSE, minYear = 2000, maxYear = 2021, recThresh = 5, env_dat, remove_duplicates = TRUE, species_data, function_loc = "/home/users/thoval/DECIDE/scripts/Edited_Rob_Functions.R"){
  
  require(raster)
  require(tidyverse)
  
  print(species_name)
  
  # load functions
  source(as.character(function_loc))
  
  # load raster
  endat <- raster::stack(as.character(env_dat))
  
  # load species data
  sp_dat <- read_csv(as.character(species_data))
  
  # get a species and year column
  sp_dat <- sp_dat %>% 
    mutate(year = yr,
           species = sp_n,
           thinned_id = paste(species, TO_GRIDREF, lon, lat))  
  
  dim(sp_dat)
  
  # remove duplicates
  if(remove_duplicates){
    print('!!!!!     DUPLICATES REMOVED    !!!!!')
    sp_dat <- sp_dat[!duplicated(sp_dat$thinned_id),]
    # dim(sp_dat)
  }
  
  
  # create the pseudoabsences and store as a named list
  ab1 <- list(cpa(spdat = sp_dat, species = species_name, 
                  matchPres = matchPres, nAbs = nAbs,
                  minYear = minYear, maxYear = maxYear, recThresh = recThresh, 
                  screenRaster = endat))
  
  names(ab1) <- species_name
  
  print('!!!!!     Finished pseudoabsences    !!!!!')
  
  
  return(ab1)
  
}


library(rslurm)
library(tidyverse)

setwd('/gws/nopw/j04/ceh_generic/thoval/DECIDE/SDMs/scripts/pseudoabsences')

# create the parameters file

# butterfly data
butt_df <- read_csv("/home/users/thoval/DECIDE/data/species_data/butterfly_EastNorths_no_duplicates_2021_12_06.csv")
butt_sp <- unique(butt_df$sp_n)

# moth data
moth_df <- read_csv("/home/users/thoval/DECIDE/data/species_data/DayFlyingMoths_EastNorths_no_duplicates_2021_12_06.csv")
moth_sp <- unique(moth_df$sp_n)

# night-flying moth data
nfmoth_df <- read_csv("/home/users/thoval/DECIDE/data/species_data/NightFlyingMoths_EastNorths_no_duplicates_2021_12_06.csv")
nfmoth_sp <- unique(nfmoth_df$sp_n)


# parameters
pars <- data.frame(species_name = c(butt_sp, moth_sp, nfmoth_sp), 
                   env_dat = '/home/users/thoval/DECIDE/data/environmental_data/envdata_fixedcoasts_nocorrs_100m_GB.grd',
                   species_data = c(rep("/home/users/thoval/DECIDE/data/species_data/butterfly_EastNorths_no_duplicates_2021_12_06.csv", length(butt_sp)),
                                    rep("/home/users/thoval/DECIDE/data/species_data/DayFlyingMoths_EastNorths_no_duplicates_2021_12_06.csv", length(moth_sp)),
                                    rep("/home/users/thoval/DECIDE/data/species_data/NightFlyingMoths_EastNorths_no_duplicates_2021_12_06.csv", length(nfmoth_sp))),
                   nAbs = 10000, 
                   matchPres = FALSE, 
                   minYear = 2000, 
                   maxYear = 2021, 
                   recThresh = 5, 
                   remove_duplicates = TRUE,
                   function_loc = "/home/users/thoval/DECIDE/scripts/Edited_Rob_Functions.R",
                   stringsAsFactors = FALSE)

dim(pars)

## slurm apply call
pseudo_slurm <- slurm_apply(slurm_pseudoabs_fun,
                            params = pars,
                            jobname = 'pseudoabsences',
                            nodes = nrow(pars),
                            cpus_per_node = 1,
                            slurm_options = list(partition = "short-serial-4hr", 
                                                 time = "03:59:59", 
                                                 mem = "20000",
                                                 error = 'pseudo_abs_err-%j-%a.out',
                                                 account = "short4hr"),
                            submit = TRUE)







#####     get results set working directory to the folder with slurm cll in it     #####

library(rslurm)
library(tidyverse)

# create the parameters file

# butterfly data
butt_df <- read_csv("/home/users/thoval/DECIDE/data/species_data/butterfly_EastNorths_no_duplicates.csv")
butt_sp <- unique(butt_df$sp_n)

# moth data
moth_df <- read_csv("/home/users/thoval/DECIDE/data/species_data/DayFlyingMoths_EastNorths_no_duplicates.csv")
moth_sp <- unique(moth_df$sp_n)

# night-flying moth data
nfmoth_df <- read_csv("/home/users/thoval/DECIDE/data/species_data/NightFlyingMoths_EastNorths_no_duplicates_2021_12_06.csv")
nfmoth_sp <- unique(nfmoth_df$sp_n)

# parameters
pars <- data.frame(species_name = c(butt_sp, moth_sp, nfmoth_sp), 
                   env_dat = '/home/users/thoval/DECIDE/data/environmental_data/envdata_fixedcoasts_nocorrs_100m_GB.grd',
                   species_data = c(rep("/home/users/thoval/DECIDE/data/species_data/butterfly_EastNorths_no_duplicates_2021_12_06.csv", length(butt_sp)),
                                    rep("/home/users/thoval/DECIDE/data/species_data/DayFlyingMoths_EastNorths_no_duplicates_2021_12_06.csv", length(moth_sp)),
                                    rep("/home/users/thoval/DECIDE/data/species_data/NightFlyingMoths_EastNorths_no_duplicates_2021_12_06.csv", length(nfmoth_sp))),
                   nAbs = 10000, 
                   matchPres = FALSE, 
                   minYear = 2000, 
                   maxYear = 2021, 
                   recThresh = 5, 
                   remove_duplicates = TRUE,
                   function_loc = "/home/users/thoval/DECIDE/scripts/Edited_Rob_Functions.R",
                   stringsAsFactors = FALSE)

# slurm call with submit=FALSE
pseudo_slurm <- slurm_apply(slurm_pseudoabs_fun,
                            params = pars,
                            jobname = 'pseudoabsences',
                            nodes = nrow(pars),
                            cpus_per_node = 1,
                            slurm_options = list(partition = "short-serial-4hr", 
                                                 time = "03:59:59", 
                                                 mem = "6000",
                                                 error = 'pseudo_abs_err-%j-%a.out',
                                                 account = "short4hr"),
                            submit = FALSE)

# decide what to call these pseudoabsences
pa_name <- 'PA_thinned_10000nAbs'

# get the outputs
res <- get_slurm_out(pseudo_slurm, outtype = 'raw')
names(res[[1]])

# combine the slurm out
res_out <- lapply(res, FUN = function(x){
  out <- x[[1]]
})

# get the species name - important in case one fails
res_nam <- lapply(res, FUN = function(x) names(x))

# give the correct names
names(res_out) <- res_nam
names(res_out)

# separate the moths and butterflies
butt_out <- res_out[names(res_out) %in% butt_sp]
moth_out <- res_out[names(res_out) %in% moth_sp]
nfmoth_out <- res_out[names(res_out) %in% nfmoth_sp]

# suffix to give the data
suffix = gsub('-', '_', Sys.Date())

# write outputs
save(butt_out, file = paste0('butterfly_', pa_name, '_', suffix, '.rdata'))
save(moth_out, file = paste0('moth_', pa_name, '_', suffix, '.rdata'))
save(moth_out, file = paste0('nightflying_moth_', pa_name, '_', suffix, '.rdata'))

# find the missing species
missing_spp_butt <- butt_sp[!butt_sp %in% names(butt_out)]
missing_spp_butt

missing_spp_moth <- moth_sp[!moth_sp %in% names(moth_out)]
missing_spp_moth

missing_spp_nf_moth <- nfmoth_sp[!nfmoth_sp %in% names(nfmoth_out)]
missing_spp_nf_moth



