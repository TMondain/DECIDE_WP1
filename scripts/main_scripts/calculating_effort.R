

####    Making effort layer for butterflies and moths

library(tidyverse)
library(raster)
library(BRCmap)


## read in raster and change to data frame
r <- raster::stack("data/environmental_data/edat_nocorrs_nosea_cropped.grd")
r <- r[[31]] # choose elevation cos likely to be most complete
r_df <- as.data.frame(r, xy=T, centroids=TRUE)[,1:2]


## read in moths
dfm <- read_csv("data/edited_insect_data/moth/DayFlyingMoths_EastNorths_no_duplicates.csv")
dfm <- cbind(dfm, gr_let2num(dfm$TO_GRIDREF, centre=T))


## get the effort layer at the current resolution (100m)
m1_eff <- dfm %>%
  group_by(EASTING,NORTHING) %>%
  summarise(effort = length(unique(date))) %>%
  arrange(-effort)
m1_eff
hist(m1_eff$effort)


## match the effort values to corresponding grid cells from the raster
r_df$effort <- m1_eff$effort[match(paste(r_df$x, r_df$y), paste(m1_eff$EASTING, m1_eff$NORTHING))]
hist(r_df$effort) # looks the same as for the data frame - good!


## convert to raster
rast_eff <- rasterFromXYZ(r_df)
# writeRaster(rast_eff, filename = 'data/edited_insect_data/moth/effort_layers/100m_effort_moths.grd',
#             format = 'raster', overwrite = T)


## aggregate error raster to 1km (factor of)
km1_eff <- raster::aggregate(rast_eff, 10, fun=sum)
# writeRaster(km1_eff, filename = 'data/edited_insect_data/moth/effort_layers/1km_effort_moths.grd',
#             format = 'raster', overwrite = T)


## aggregate raster to 10km (factor of 10)
km10_eff <- raster::aggregate(km1_eff, 10, fun=sum)
# writeRaster(km10_eff, filename = 'data/edited_insect_data/moth/effort_layers/10km_effort_moths.grd',
#             format = 'raster', overwrite = T)


#### Spatial Kriging to get an effort layer to downweight DECIDE score 
eff_100_krig <- focal(x, 
                      w = matrix(c(0.3,0.3,0.3,0.3,1,
                                   0.3,0.3,0.3,0.3), nrow = 3, ncol = 3))



#####     Now for butterflies

## read in butterflies
bdf <- read_csv("data/edited_insect_data/butterfly/BNM_NRMS_No_Duplicates_east_north.csv")
bdf <- cbind(bdf, gr_let2num(bdf$TO_GRIDREF))


## get the effort layer at the current resolution (100m)
m1_eff_b <- bdf %>%
  group_by(EASTING,NORTHING) %>%
  summarise(effort = length(unique(date))) %>%
  arrange(-effort)
hist(m1_eff_b$effort)


## match the effort values to corresponding grid cells from the raster
## replace the moth effort column - overwrite it
r_df$effort <- m1_eff_b$effort[match(paste0(r_df$x, r_df$y), paste0(m1_eff_b$EASTING, m1_eff_b$NORTHING))]
hist(r_df$effort) # looks the same as for the data frame - good!

## convert to raster
rast_eff_butt <- rasterFromXYZ(r_df)
# writeRaster(rast_eff_butt, filename = 'data/edited_insect_data/butterfly/effort_layers/100m_effort_butterfly.grd',
#             format = 'raster', overwrite = T)


## aggregate error raster to 1km (factor of)
km1_eff_butt <- raster::aggregate(rast_eff_butt, 10, fun=sum)
# writeRaster(km1_eff_butt, filename = 'data/edited_insect_data/butterfly/effort_layers/1km_effort_butterfly.grd',
#             format = 'raster', overwrite = T)


## aggregate raster to 10km (factor of 10)
km10_eff_butt <- raster::aggregate(km1_eff_butt, 10, fun=sum)
# writeRaster(km10_eff_butt, filename = 'data/edited_insect_data/butterfly/effort_layers/10km_effort_butterfly.grd',
#             format = 'raster', overwrite = T)

