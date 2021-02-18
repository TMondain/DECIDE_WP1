
######      Test script for the SDM functions      ######

library(readr)
library(tidyverse)
library(raster)
library(foreach)
library(doParallel)
library(dismo)
library(randomForest)
library(rfinterval)
library(BRCmap)
library(ranger)
library(mgcv)
source('scripts/functions/Edited_Rob_Functions.R')


ed <- raster::stack("data/environmental_data/lcm_had_elev_national_grid.gri")
names(ed)

## big problem with slope and aspect!!!
ed <- dropLayer(ed, i = match(c("slope", "aspect"), names(ed)))

slope <- terrain(ed[[42]], opt = "slope", unit = "degrees")
aspect <- terrain(ed[[42]], opt = "aspect", unit = "degrees")

ed <- raster::stack(ed, slope, aspect)
names(ed)

plot(ed[[44]], main = names(ed[[44]]), ylim=c(780000, 860000), xlim = c(60000, 200000))
plot(ed[[44]])

getwd()
# save(ed, file = "AllEnvironmentalData.rdata")
# load("AllEnvironmentalData.rdata")

# slp_asp <- raster::stack(slope, aspect)
# save(slp_asp, file = "slope_aspect.rdata")
# load("slope_aspect.rdata")
# plot(slp_asp)


## cropping to a small extent
ext_h <- extent(matrix(c(-1,53, 0.2,54.5), ncol = 2))
e <- as(ext_h, "SpatialPolygons")
sp::proj4string(e) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

e.geo <- sp::spTransform(e, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))

hbv_y <- raster::crop(ed, e.geo)

plot(hbv_y[[43]], main = names(hbv_y[[43]]))

names(hbv_y)


####    Remove correlated variables    ####
# exclude variables with >0.7 correlation
whichVars <- usdm::vifcor(hbv_y[[26:41]], th = 0.7)
whichVars

# # excluding those with >10 vif stepwise
# whichVars2 <- usdm::vifstep(hbv_y, th = 10)
# whichVars2

whichVars@excluded

hbv_y <- dropLayer(x=hbv_y, i = match(whichVars@excluded, names(hbv_y)))
# plot(hbv_y)
names(hbv_y)
ht <- hbv_y

plot(ht[[1]])
names(ht)

#####    Load in day flying moth data     #####
dfm_df <- read_csv("data/edited_insect_data/moth/DayFlyingMoths.csv")
head(dfm_df)

# get eastings and northings
dfm_df <- c_en(dfm_df)


# numbers all uk
dfm_df %>% group_by(sp_n) %>% tally %>% 
  ggplot(aes(x = reorder(sp_n, n), y = n)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 50, hjust=1)) +
  xlab("species") + ggtitle("All UK")

## subset to AOI
sp_y <- subset(dfm_df, lat > 53.1 & lat <= 54.4 &
                 lon > -0.8 & lon < 0.2) %>% 
  mutate(species = sp_n,
         year = yr)

# New tally
sp_y %>% group_by(sp_n) %>% tally %>% 
  ggplot(aes(x = reorder(sp_n, n), y = n)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 50, hjust=1)) +
  xlab("species") + ggtitle("Subset UK")

# View(sp_y %>% group_by(sp_n) %>% tally %>%
#        arrange(n))

#### Create presence/absence
spp <- unique(sp_y$species)

# as a hack to not have to create a new function for east and north
# create a data frame that has East and north renamed as lon lat
ndf <- sp_y %>% mutate(lon = EASTING,
                       lat = NORTHING) %>% 
  dplyr::select(-EASTING, -NORTHING)


registerDoParallel(cores = detectCores() - 1)

system.time( 
  ab1 <- foreach(i = 1 : length(spp)) %dopar% {
    
    cpa(spdat = ndf, species = spp[i], matchPres = TRUE,
        minYear = 2000, maxYear = 2017, recThresh = 5,
        screenRaster = ht)
    
  }
)

registerDoSEQ()


## test functions with for loop

names(ab1) <- spp

spp_lr_out <- list()

system.time(
  for(s in 2:3){
    print(paste(s, spp[s], sep = " "))
    
    if(is.null(ab1[[s]])){
      print(paste("No data for species", spp[s]))
      next
    }
    
    sdm_lr <- fsdm(species = spp[s], model = "lr",
                   climDat = ht, spData = ab1, knots = -1,
                   k = 10,
                   write =  F, outPath = "C:/Users/thoval/Documents/Analyses/lr_outs/")
    
    # se_out <- predict(ht, sdm_lr$Model, type = "response", se.fit = TRUE, index = 2)
    # save(se_out, file = paste0("C:/Users/thoval/Documents/Analyses/lr_outs/", spp[s], "_lr_se_error_prediction.rdata"))
    
    spp_lr_out[[s]] <- sdm_lr
    
  }
)

spp_lr_out[[3]]$AUC


## predict from each of the bootstrapped models
boots_2 <- lapply(spp_lr_out[[2]]$Bootstrapped_models, 
                  FUN = function(x) predict(ht, x, type='response', index=NULL))
boot_out <- do.call("stack", boots_2)

## quantiles
print(paste0('#####   getting quantiles   #####'))
mean_preds <- calc(boot_out, fun = mean, na.rm = T)
quant_preds <- calc(boot_out, fun = function(x) {quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)})
rnge <- quant_preds[[2]]-quant_preds[[1]]

par(mfrow=c(1,2))
plot(mean_preds)
plot(rnge)
par(mfrow=c(1,1))


#######################################################
## testing the bootstrapped model errors for logistic regression

pred_out_boot <- list()

for(b in 1:length(spp_lr_out[[1]]$Bootstrapped_models)) {
  print(b)
  
  b_t <- spp_lr_out[[1]]$Bootstrapped_models[[b]]
  p_t <- predict(ht, b_t, index = NULL, type = 'response', ext = NULL)
  pred_out_boot[[b]] <- p_t
  
}

boots <- do.call("stack", pred_out_boot)
b_r <- range(boots)
plot(b_r[[2]]-b_r[[1]])



# do all species in parallel to save time
registerDoParallel(cores = detectCores() - 1)

names(ab1) <- spp
system.time(
  spp_lr_out <- foreach(s = 1:length(spp), #.combine='comb', .multicombine=TRUE,
                        # .init=list(list()), 
                        .packages = c('tidyverse', 'raster', 'readr', 'dismo', 'ranger'),
                        .errorhandling = 'pass') %dopar% {
                          
                          # print(paste(s, spp[s], sep = " "))
                          
                          if(is.null(ab1[s])){
                            # print(paste("No data for species", spp[s]))
                            next
                          }
                          
                          sdm_lr <- fsdm(species = spp[s], model = "lr", 
                                         climDat = ht, spData = ab1, knots = -1,
                                         k = 20,
                                         prediction = T,
                                         inters = F, 
                                         write =  F, outPath = "C:/Users/thoval/Documents/Analyses/lr_outs/")
                          
                          se_out <- predict(ht, sdm_lr$Model, type = "response", se.fit = TRUE, index = 2)
                          # save(se_out, file = paste0("C:/Users/thoval/Documents/Analyses/lr_outs/", spp[s], "_lr_se_error_prediction.rdata"))
                          
                          list(sdm_lr, se_out)
                        }
)

registerDoSEQ()
par(mfrow = c(1,1))
plot(spp_lr_out[[1]][[1]]$Predictions,main = spp_lr_out[[1]][[1]]$Species)

spp_out_boot <- list()

for(s in 1:length(spp_lr_out)){
  
  # some entries in the list are NULL because of too few data points in the kfold cross validation part of the model above
  skip_to_next <- FALSE
  
  # find the errror
  tryCatch({ 
    bt_mods_loc <- spp_lr_out[[s]][[1]]$Bootstrapped_models
  },
  error = function(e) {skip_to_next <<- TRUE})
  
  # skip to next if error
  if(skip_to_next){print(paste("SKIPPING LIST ENTRY", s))
    next}
  pred_out_boot <- list()
  for(b in 1:length(spp_lr_out[[s]][[1]]$Bootstrapped_models)) {
    print(paste(s, spp_lr_out[[s]][[1]]$Species))
    print(b)
    
    # predict from each bootsrapped model
    b_t <- spp_lr_out[[s]][[1]]$Bootstrapped_models[[b]]
    p_t <- predict(ht, b_t, index = NULL, type = 'response', ext = NULL)
    pred_out_boot[[b]] <- p_t
    
  }
  
  # stack all the predictions and store as outputs
  boots <- do.call("stack", pred_out_boot)
  spp_out_boot[[spp_lr_out[[s]][[1]]$Species]] <- boots
  
}

spp_out_boot$`Zygaena lonicerae` 

zl_q <- calc(spp_out_boot$`Zygaena lonicerae` , fun=quantile, na.rm = T)
zl_q

seq(0, 1, 0.25)

zl_q_25_75 <- zl_q[[4]]-zl_q[[1]]

par(mfrow = c(1,2))
plot(spp_lr_out[[s]][[1]]$Predictions)
# points(x = spp_lr_out[[s]][[1]]$Data$lon[spp_lr_out[[s]][[1]]$Data$val == 1], 
#        spp_lr_out[[s]][[1]]$Data$lat[spp_lr_out[[s]][[1]]$Data$val == 1],
#        cex = 0.6, col = "red", pch = 20)
plot(zl_q_25_75)
# points(x = spp_lr_out[[s]][[1]]$Data$lon[spp_lr_out[[s]][[1]]$Data$val == 1], 
#        spp_lr_out[[s]][[1]]$Data$lat[spp_lr_out[[s]][[1]]$Data$val == 1],
#        cex = 0.6, col = "red", pch = 20)
par(mfrow = c(1,1))










x <- 
  foreach(s = 1:length(spp_lr_out), .packages = 'raster') %:%
  foreach(b = 1:20, .combine = 'stack') %dopar% {
    
    b_t <- spp_lr_out[[s]][[1]]$Bootstrapped_models[[b]]
    p_t <- predict(ht, b_t, index = NULL, type = 'response', ext = NULL)
    p_t
  }

# plot all the saved linear regression outputs
# assumes they are loaded in alphabetically!
fls_lr <- list.files("C:/Users/thoval/Documents/Analyses/lr_outs/", full.names = T)
par(mfrow = c(2,2))
for(j in seq(1,40, by =2)){
  load(fls_lr[j])
  load(fls_lr[j+1])
  
  plot(out$Predictions, main = paste(out$Species, "AUC", round(out$AUC,4)))
  points(x = out$Data$lon[out$Data$val == 1], out$Data$lat[out$Data$val == 1],
         cex = 0.6, col = "red", pch = 20)
  
  plot(se_out, main = paste(out$Species, "SE"))
  # points(x = out$Data$lon[out$Data$val == 1], out$Data$lat[out$Data$val == 1],
  #        cex = 0.6, col = "red", pch = 20)
  
}
par(mfrow = c(1,1))
plot(ht[[33]])

length(spp)

par(mfrow = c(4,4))

for(p in 1:length(spp)) {
  if(is.null(spp_lr_out[[p]])){
    next
  }
  
  print(plot(spp_lr_out[[p]]$Predictions, main = spp_lr_out[[p]]$Species))
  
}




########################
### Random Forest

names(ab1) <- spp
ab1$`Tyria jacobaeae`

## remove null items from list
## ones that don't have enough data
ab1[sapply(ab1, is.null)] <- NULL
names(ab1)


cores<-7
cl <- makeCluster(cores)
registerDoParallel(cl)
# registerDoParallel(cores = detectCores() - 1)

system.time(
  spp_rf_out <- 
    foreach(s = 1:length(spp), #.combine='comb', .multicombine=TRUE,
            # .init=list(list()), 
            .packages = c('tidyverse', 'raster', 'readr', 'dismo', 'randomForest', 'ranger'),
            .errorhandling = 'pass') %dopar% {
              
              # print(paste(s, spp[s], sep = " "))
              
              if(is.null(ab1[[s]])){
                # print(paste("No data for species", spp[s]))
                next
              }
              
              sdm_rf <- fsdm(species = spp[s], model = "rf", 
                             climDat = ht, spData = ab1, knots = -1,
                             k = 5, 
                             inters = F, 
                             write =  T, outPath = "C:/Users/thoval/Documents/Analyses/rf_outs/")
              
              # se_out <- predict(ht, sdm_lr$Model, type = "response", se.fit = TRUE, index = 2)
              # save(se_out, file = paste0("C:/Users/thoval/Documents/Analyses/lr_outs/", spp[s], "_lr_se_error_prediction.rdata"))
              
              list(sdm_rf)#, se_out)
            }
)

registerDoSEQ()
length(spp_rf_out)

par(mfrow = c(1,2))

for(j in 30){
  
  plot(spp_rf_out[[j]][[1]]$Predictions, main = spp_rf_out[[j]][[1]]$Species) 
  plot(spp_rf_out[[j]][[2]], main = paste(spp_rf_out[[j]][[1]]$Species, "SE")) 
  
}

plot(spp_rf_out[[2]][[1]]$Predictions, main = spp_rf_out[[2]][[1]]$Species) 
plot(spp_rf_out[[2]][[2]], main = spp_rf_out[[2]][[1]]$Species) 


plot(spp_rf_out[[3]][[1]]$Predictions, main = spp_rf_out[[3]][[1]]$Species) 
plot(spp_rf_out[[3]][[2]], main = spp_rf_out[[3]][[1]]$Species) 



fls <- list.files("C:/Users/thoval/Documents/Analyses/rf_outs/", full.names = T)[c(1:11, 13:26)]

par(mfrow = c(2,2))

for(i in fls){
  load(i)
  plot(out$Predictions, main = paste(out$Species, "AUC", round(out$AUC,4)))
  points(x = out$Data$lon[out$Data$val == 1], out$Data$lat[out$Data$val == 1],
         cex = 0.6, col = "red", pch = 20)
}


data.frame(imp = importance(out$Model)) %>% 
  rownames_to_column(var = 'var') %>% 
  ggplot(aes(x = reorder(var, imp), y = imp)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  ylab('Variable importance') +
  xlab('Variable')



#############################################
#### GAMs

cores<-7
cl <- makeCluster(cores)
registerDoParallel(cl)
# registerDoParallel(cores = detectCores() - 1)

system.time(
  spp_gam_out <- 
    foreach(s = 1:length(spp), #.combine='comb', .multicombine=TRUE,
            # .init=list(list()), 
            .packages = c('tidyverse', 'raster', 'readr', 'dismo', 'randomForest', 'mgcv'),
            .errorhandling = 'pass') %dopar% {
              
              # print(paste(s, spp[s], sep = " "))
              
              if(is.null(ab1[[s]])){
                # print(paste("No data for species", spp[s]))
                next
              }
              
              sdm_gam <- fsdm(species = spp[s], model = "gam", 
                              climDat = ht, spData = ab1, knots = -1,
                              k = 5, 
                              inters = F, 
                              write =  T, outPath = "C:/Users/thoval/Documents/Analyses/gam_outs/")
              
              se_out_gam <- predict(ht, sdm_gam$Model, type = "response", se.fit = TRUE, index = 2)
              save(se_out_gam, file = paste0("C:/Users/thoval/Documents/Analyses/gam_outs/", spp[s], "_gam_se_error_prediction.rdata"))
              
              list(sdm_gam, se_out_gam)
            }
)

registerDoSEQ()

(spp_gam_out)

par(mfrow = c(1,2))

for(j in 1:30){
  
  plot(spp_gam_out[[j]][[1]]$Predictions, main = spp_gam_out[[j]][[1]]$Species) 
  plot(spp_gam_out[[j]][[2]], main = paste(spp_gam_out[[j]][[1]]$Species, "SE")) 
  
}


fls_gam <- list.files("C:/Users/thoval/Documents/Analyses/gam_outs/", full.names = T)

par(mfrow = c(2,2))

for(i in fls_gam){
  load(i)
  plot(out$Predictions, main = paste(out$Species, "AUC", round(out$AUC,4), "n =", out$`Number of records`))
  points(x = out$Data$lon[out$Data$val == 1], out$Data$lat[out$Data$val == 1],
         cex = 0.6, col = "red", pch = 20)
}


# plot all the saved GAM outputs
# assumes they are loaded in alphabetically!
fls_gam <- list.files("C:/Users/thoval/Documents/Analyses/gam_outs/", full.names = T)

par(mfrow = c(2,2))
for(j in seq(1,40, by =2)){
  load(fls_gam[j])
  load(fls_gam[j+1])
  
  plot(out$Predictions, main = paste(out$Species, "AUC", round(out$AUC,4)))
  points(x = out$Data$lon[out$Data$val == 1], out$Data$lat[out$Data$val == 1],
         cex = 0.6, col = "red", pch = 20)
  
  plot(se_out_gam, main = paste(out$Species, "SE"))
  # points(x = out$Data$lon[out$Data$val == 1], out$Data$lat[out$Data$val == 1],
  #        cex = 0.6, col = "red", pch = 20)
  
}
par(mfrow = c(1,1))
