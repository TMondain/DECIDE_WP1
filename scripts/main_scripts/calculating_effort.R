

####    Making effort layer for butterflies and moths

library(tidyverse)
library(raster)
library(BRCmap)


## read in raster and change to data frame
r <- raster::stack("data/environmental_data/edat_nocorrs_nosea.grd")
r <- r[[31]] # choose elevation cos likely to be most complete
r_df <- as.data.frame(r, xy=T)[,1:2]


## read in moths
dfm <- read_csv("data/edited_insect_data/moth/DayFlyingMoths.csv")
dfm$lon


## get the effort layer at the current resolution (100m)
m1_eff <- dfm %>%
  group_by(lon,lat) %>%
  summarise(effort = length(unique(date))) %>%
  arrange(-effort)
m1_eff
hist(m1_eff$effort)


## match the effort values to corresponding grid cells from the 
r_df$effort <- m1_eff$effort[match(paste0(r_df$x, r_df$y), paste0(km1_eff$lon, km1_eff$lat))]
hist(r_df$effort) # looks the same as for the data frame - good!

## convert to raster
rast_eff <- rasterFromXYZ(r_df)
writeRaster(rast_eff, filename = 'data/edited_insect_data/moth/100m_effort_moths.grd',
            format = 'raster', overwrite = T)

## aggregate error raster to 1km (factor of)
km1_eff <- raster::aggregate(rast_eff, 10, fun=sum)
writeRaster(rast_eff, filename = 'data/edited_insect_data/moth/100m_effort_moths.grd',
            format = 'raster', overwrite = T)


## aggregate raster to 10km (factor of 10)
km10_eff <- raster::aggregate(km1_eff, 10, fun=sum)
writeRaster(rast_eff, filename = 'data/edited_insect_data/moth/100m_effort_moths.grd',
            format = 'raster', overwrite = T)

