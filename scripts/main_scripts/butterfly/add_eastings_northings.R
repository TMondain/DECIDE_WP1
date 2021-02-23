
## transforming butterfly data to east-norths

library(readr)
library(BRCmap)
library(tidyverse)
library(sp)
library(rgdal)

source('scripts/functions/Edited_Rob_Functions.R')


## load in butterfly data with no duplicates
bdf <- read_csv('data/edited_insect_data/butterfly/BNM_NRMS_No_Duplicates.csv')
head(bdf)

bdf <- c_en(bdf)

# as a hack to not have to create a new function for east and north
# create a data frame that has East and north renamed as lon lat
ndf <- bdf %>% mutate(lon = EASTING,
                      lat = NORTHING) %>% 
  dplyr::select(-EASTING, -NORTHING)

head(ndf)

# write file
write.csv(ndf, file = 'data/edited_insect_data/butterfly/BNM_NRMS_No_Duplicates_east_north.csv')
