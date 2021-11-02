
library(raster)
library(dismo)


########################################
#####     HadUK weather variables


# Load all the weather data into R
# have files for 2010 to 2016 inclusive of monthly min and max temp and rainfall
# setwd("C:/Users/thoval/OneDrive - NERC/Documents/Analyses/")

path <- "../../../Documents/DECIDE/DECIDE_WP1/data/raw_data/environmental/HadUK_dat/"

vars <- c("rainfall", "tasmax", "tasmin")

out_var <- list()

for(i in vars) {
  
  # list all the files with a given structure
  files <- list.files(pattern = i, path = path, full.names = T)
  
  # initialise list
  st_l <- list()
  
  for(j in 1:length(files)){
    
    # get the files for each year
    # need to use stack so that each month is a different layer
    env_d <- raster::stack(files[j])
    
    # rename layers 1:12 - months of the year
    names(env_d) <- paste("M", seq(1:12), sep = "_")  
    
    # store as an object to use as an index to average across later
    ind <- names(env_d)
    
    # store the object in the list to use later
    st_l[[j]] <- env_d
    
    
  }
  
  print(i)
  
  # now want to average the same layer across multiple raster stacks
  # basically because I want to get the average min, max temp and rainfall 
  # for each month over 6 years, 2010 to 2015 inclusive
  st_o <- stackApply(stack(st_l), # creates a stack of all objects in the list
                     fun = mean, # takes the mean
                     indices = ind, # provide an index to get the names
                     na.rm = T)
  out_var[[i]] <- st_o
  
}

# plot rainfall in jan and july
plot(out_var$rainfall[[1]])
plot(out_var$rainfall[[7]])
# looks sensible

# # project them to the "right" coords system
# basically so I can refer to lat lon for plotting
rain <- projectRaster(out_var$rainfall, crs = "+proj=longlat +datum=WGS84 +no_defs", res = 0.008333333)
tasmax <- projectRaster(out_var$tasmax, crs = "+proj=longlat +datum=WGS84 +no_defs", res = 0.008333333)
tasmin <- projectRaster(out_var$tasmin, crs = "+proj=longlat +datum=WGS84 +no_defs", res = 0.008333333)

# crop the files to same extent as LCM
ext <- extent(matrix(c(-9,49, 2,60), ncol = 2))
rain <- crop(rain, ext)
tasmax <- crop(tasmax, ext)
tasmin <- crop(tasmin, ext)
plot(tasmin)

# convert to the 19 bioclimatic variables
had_bv <- dismo::biovars(prec = rain, tmin = tasmin, tmax = tasmax)

# disaggregate those variables by a factor of 10 using bilinear interpolation
# might be worth thinking about whether or not to do it bilinearly...
had_bv_fine <- disaggregate(had_bv, fact = 10, fun = 'bilinear')

# name them their real names
names(had_bv_fine) <- c("AnnualTemp",
                        "MeanDiRange",
                        "Isotherm",
                        "TempSeasonality",
                        "MaxTempWarmestMonth",
                        "MinTempColdestMonth",
                        "TempAnnualRange",
                        "MeanTempWetQuarter",
                        "MeanTempDriestQuarter",
                        "MeanTempWarmQuarter",
                        "MeanTempColdQuarter",
                        "AnnualPrecip",
                        "PrecipWetMonth",
                        "PrecipDriestMonth",
                        "PrecipSeasonality",
                        "PrecipWettestQuarter",
                        "PrecipDriestQuarter",
                        "PrecipWarmQuarter",
                        "PrecipColdQuarter")



# # save the output as a RASTER
# writeRaster(had_bv_fine, filename = "Data/HadUK_dat/had_bv_fine.grd", format = "raster")
had_bv_fine <- raster::stack("Data/HadUK_dat/had_bv_fine.grd")

plot(had_bv_fine[[1]])

# ## This was from before, disaggregated the raw weather variables to 100m resolution
# ## but then realised the biovars function wouldn't work
# rain_fine <- disaggregate(rain, fact = 10, fun = 'bilinear')
# tasmax_fine <- disaggregate(tasmax, fact = 10, fun = 'bilinear')
# tasmin_fine <- disaggregate(tasmin, fact = 10, fun = 'bilinear')


#######################
### National grid CRS

# convert to the 19 bioclimatic variables
had_bv <- dismo::biovars(prec = out_var$rainfall, tmin = out_var$tasmin, tmax = out_var$tasmax)

# name real names
names(had_bv) <- c("AnnualTemp",
                   "MeanDiRange",
                   "Isotherm",
                   "TempSeasonality",
                   "MaxTempWarmestMonth",
                   "MinTempColdestMonth",
                   "TempAnnualRange",
                   "MeanTempWetQuarter",
                   "MeanTempDriestQuarter",
                   "MeanTempWarmQuarter",
                   "MeanTempColdQuarter",
                   "AnnualPrecip",
                   "PrecipWetMonth",
                   "PrecipDriestMonth",
                   "PrecipSeasonality",
                   "PrecipWettestQuarter",
                   "PrecipDriestQuarter",
                   "PrecipWarmQuarter",
                   "PrecipColdQuarter")

# # store the coarse UK
# writeRaster(had_bv, filename = "Data/HadUK_dat/had_bv_1km_national_grid.grd", format = "raster", overwrite = T)

plot(had_bv)

# disaggregate those variables by a factor of 10 using bilinear interpolation
# might be worth thinking about whether or not to do it bilinearly...
had_bv_fine <- disaggregate(had_bv, fact = 10, fun = 'bilinear')

# name them their real names
names(had_bv_fine) <- c("AnnualTemp",
                        "MeanDiRange",
                        "Isotherm",
                        "TempSeasonality",
                        "MaxTempWarmestMonth",
                        "MinTempColdestMonth",
                        "TempAnnualRange",
                        "MeanTempWetQuarter",
                        "MeanTempDriestQuarter",
                        "MeanTempWarmQuarter",
                        "MeanTempColdQuarter",
                        "AnnualPrecip",
                        "PrecipWetMonth",
                        "PrecipDriestMonth",
                        "PrecipSeasonality",
                        "PrecipWettestQuarter",
                        "PrecipDriestQuarter",
                        "PrecipWarmQuarter",
                        "PrecipColdQuarter")

# writeRaster(had_bv_fine, filename = "Data/HadUK_dat/had_bv_fine_national_grid.grd", format = "raster")
hbv <- raster::stack("Data/HadUK_dat/had_bv_fine_national_grid.grd")
hbv

######################################
#####     ELEVATION


el1 <- raster("Data/Copernicus_Elevation/eu_dem_v11_E30N30/eu_dem_v11_E30N30.TIF")
el2 <- raster("Data/Copernicus_Elevation/eu_dem_v11_E30N40/eu_dem_v11_E30N40.TIF")

ext_no_ice <- extent(matrix(c(3000000,3000000, 3900000, 4500000), ncol = 2))
el2_n_ice <- crop(el2, ext_no_ice)
plot(el2_n_ice)


# aggregate to get 100m resolution - mean across all cells
el1_coarse <- raster::aggregate(el1, fact = 4, fun = mean)
el2_coarse <- raster::aggregate(el2_n_ice, fact = 4, fun = mean)

# merge the rasters together
elev_100 <- merge(el1_coarse, el2_coarse)

# remove iceland
ext_proj <- extent(matrix(c(3000000,3000000, 3900000, 4500000), ncol = 2))
elev_100_crop <- crop(elev_100, ext_proj) 
plot(elev_100_crop)

# writeRaster(elev_100_crop, filename = "Data/Copernicus_Elevation/elevation_UK.tif", format = "GTiff")
elev_100_crop <- raster("Data/Copernicus_Elevation/elevation_UK.tif")

##############################
##### Land cover map

# percentage cover of each land cover map
lcm_p <- raster::stack("Data/CEH_land_cov/lcm_2015/GB/data/lcm2015gb100perc.tif")
names(lcm_p) <-  c("sea",
                   "broad_wood",
                   "conif_wood",
                   "arable",
                   "impr_grass",
                   "neutr_grass",
                   "calc_grass",
                   "acid_grass",
                   "fen_marsh_swamp",
                   "heather",
                   "heath_grass",
                   "bog",
                   "inland_rock",
                   "saltwater",
                   "freshwater",
                   "sup_lit_rock",
                   "sup_lit_sed",
                   "lit_rock",
                   "lit_sed",
                   "saltmarsh",
                   "urban",
                   "suburban")

# aggregate(lcm_m, fact = 10, fun = mean)
plot(lcm_p[[1]])

# get landcover and had vars on same extent
cr_uk <- extent(lcm_p)
had_bv_fine_cropped <- crop(hbv, cr_uk)
had_bv_fine_cropped
lcm_p

# still not the same, lcm goes further north
cr_uk2 <- extent(had_bv_fine_cropped)
lcm_p_crop <- crop(lcm_p, cr_uk2)

# stack them together
lcm_had_national_grid <- raster::stack(lcm_p_crop, had_bv_fine_cropped)


# writeRaster(lcm_had_national_grid, filename = "Data/lcm_had_national_grid.grd", format = "raster")
# load raster
lcm_had <- raster::stack("Data/lcm_had_national_grid.grd")
lcm_had

plot(lcm_had[[24]], main = names(lcm_had[[24]]))

# project elevation onto national grid and same extent as LCM and hadUK variables
beginCluster()
elev <- raster::projectRaster(elev_100_crop, lcm_had, method = 'bilinear')
endCluster()

lcm_had_elev_national_grid <- raster::stack(lcm_had, elev)
names(lcm_had_elev_national_grid)

slope_asp <- terrain(lcm_had_elev_national_grid[[41]], opt = c('slope', 'aspect'), unit = 'degrees', neighbors = 8)
lcm_had_elev_national_grid <- raster::stack(lcm_had_elev_national_grid, slope_asp)


# writeRaster(lcm_had_elev_national_grid, filename = "Data/lcm_had_elev_national_grid.grd", format = "raster")
test <- raster::stack("Data/lcm_had_elev_national_grid.grd")
names(test)
test

# high elev and low precip
par(mfrow = c(1,2));plot(test[[34]], main = names(test[[34]])); plot(test[[42]], main = names(test[[42]])); par(mfrow = c(1,1))

p <- as.data.frame(test[[34]], xy = T)
e <- as.data.frame(test[[42]], xy = T)
head(p)
head(e)

df <- cbind(p, elev =e[,3]) %>% na.omit()
head(df)

hist(df$elev)
hist(df$AnnualPrecip)


library(tidyverse)

ggplot(subset(df, elev>200), aes(y = AnnualPrecip, x = elev)) +
  geom_point()





#############################################
#### combine them all OLD
beginCluster()
elev <- projectRaster(elev_100_crop, lcm_had_national_grid)
endCluster()

beginCluster()
lcm <- projectRaster(lcm_p, had_bv_fine)
endCluster()

# stack all the files together
e_dat_fine <- raster::stack(had_bv_fine, lcm, elev)
e_dat_fine

writeRaster(e_dat_fine, filename = "Data/finescale_hadbiovars_LCMperc_elev.grd", format = "raster")
?writeRaster



e_100_proj <- projectRaster(elev_100_crop, 
                            crs = "+proj=longlat +datum=WGS84 +no_defs", res = 0.0008333333)
e_100_proj
plot(e_100_proj)

ext <- extent(matrix(c(-9,4.98, 1.98,59.8), ncol = 2))
e_100_proj_crop <- crop(e_100_proj, ext, snap = 'in')
had_bv_fine_cr2 <- crop(had_bv_fine, ext, snap = 'in')


hbv_elev <- raster::stack(e_100_proj_crop, had_bv_fine_cr2)
