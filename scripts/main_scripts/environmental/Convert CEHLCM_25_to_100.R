
### Convert CEH LCM from 25m mode raster to 100m mode and 100m percentage cover


library(tidyverse)
library(raster)
library(gtools)
library(rgdal)
library(foreach)
library(doParallel)

list.files()

lcm_25 <- raster("Data/CEH_land_cov/lcm_2015/GB/data/lcm2015gb25m.tif")
lcm_25


xmin = 200000
xmax = 300000
ymin = 600000
ymax = 700000

xmin = 200000
xmax = 205000
ymin = 598000
ymax = 600000

ext <- raster::extent(xmin, xmax, ymin, ymax)

t <- raster::crop(lcm_25, ext)
plot(t)
sort(unique(values(t)))

# 0 is sea but need it to have an index
values(t) <- values(t) + 1

# create mode raster at 100m res
t_100_m <- raster::aggregate(t, fact = 4, FUN = mode, na.rm = T)
t_100_m

# percentage cover function using apply - takes a LONG time on test raster
# so probably need to parallelise it to make it better
t_100_p <- lapply(sort(unique(t)), function(land_class) {
  aggregate(t, fact=4, fun=function(vals, na.rm) {
    sum(vals==land_class, na.rm=na.rm)/length(vals)*100
  })
})

# write lapply as a for loop to get for the foreach
# function to work
o <- list()
# try to parellise
# first do normal for loop

# for each of the unique land cover variables (numbers 1-22 (number 1 is sea which is normally 0))
# get the sum of all of its value in the are of interest (4 times bigger than original) and
# divide by the number of cells in that block. Gives percentage of each land class in the UK
for(i in unique(t)) {
  
  ag <- aggregate(t, fact=4, expand = T, fun=function(vals, na.rm) {
    sum(vals==i, na.rm=na.rm)/length(vals)*100
  })
  
  o[[i]] <- ag
}

lc_p <- do.call("stack",o)
unique(lc_p)

sum(values(t) == 8)


# I think this works
# and outputs a list of aggregated rasterlayers
registerDoParallel(cores = detectCores() - 1)

system.time(
  ag <- foreach(i = unique(t)) %dopar% {
    
    
    a <- aggregate(t, fact=4, fun=function(vals, na.rm) {
      sum(vals==i, na.rm=na.rm)/length(vals)*100
      
      
    })
    
  })

ag
registerDoSEQ()


lcm_25

# create mode raster at 100m res
lcm_100_mode <- raster::aggregate(lcm_25, fact = 4, FUN = mode, na.rm = T)

# writeRaster(lcm_100_mode, filename = "Data/CEH_land_cov/lcm_2015/GB/data/lcm2015gb100mode.tif",
#             format = "GTiff")


registerDoParallel(cores = detectCores() - 1)

system.time(
  ag <- foreach(i = unique(lcm_25)) %dopar% {
    
    
    a <- aggregate(lcm_25, fact=4, fun=function(vals, na.rm) {
      sum(vals==i, na.rm=na.rm)/length(vals)*100
      
      
    })
    
  })

registerDoSEQ()


lc_p <- do.call("stack",ag)

names(lc_p) <- c("sea",
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

# writeRaster(lc_p, filename = "Data/CEH_land_cov/lcm_2015/GB/data/lcm2015gb100perc.tif",
#             format = "GTiff", bylayer = T)

t <- raster::stack("Data/CEH_land_cov/lcm_2015/GB/data/lcm2015gb100perc.tif")
names(t) <-  c("sea",
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

# whooop it works

plot(t[[1]])

