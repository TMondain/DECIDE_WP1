

#### Outline for SDMs
library(readr)
library(raster)
library(tidyverse)
library(BRCmap)
library(parallel)
library(foreach)
library(dismo)
library(randomForest)
library(mgcv)
library(glmnet)

source('scripts/functions/Edited_Rob_Functions.R')


#' ## load in environmental data
#' do whatever you want to do to it

ed <- raster::stack("data/environmental_data/edat_nocorrs_nosea.gri")

#' e,g. crop to a smaller extent
#' 

## cropping to a small extent
## based on lat lon values (which make more sense to me)
ext_h <- extent(matrix(c(-1,53, 0.2,54.5), ncol = 2))
e <- as(ext_h, "SpatialPolygons")
sp::proj4string(e) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

e.geo <- sp::spTransform(e, CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs"))

ht <- raster::crop(ed, e.geo)


#' ##    Load in day flying moth data
#' 
dfm_df <- read_csv("data/edited_insect_data/moth/DayFlyingMoths.csv")
head(dfm_df)

#' Get eastings and northings from the grid references using this function
#'  
dfm_df <- c_en(dfm_df)

#' # crop to the same small extent
#' 

sp_y <- subset(dfm_df, lat > 53.1 & lat <= 54.4 &
                 lon > -0.8 & lon < 0.2) %>%
  mutate(species = sp_n,
         year = yr)

#' This is a bit of a hack. I started using these functions when I was coding 
#' everything in lat lon and changed the functions accordingly.
#' Eventually, I need to go and convert all the functions to use eastings and
#' northings (like they were coded originally).
#' But for now, create a data frame that has East and north renamed as lon lat
ndf <- sp_y %>% mutate(lon = EASTING,
                       lat = NORTHING) %>% 
  dplyr::select(-EASTING, -NORTHING)


#' #### Create presence/absence
spp <- unique(sp_y$species)

#' # # if you want to remove duplicates from the data like we discussed
#' ### removing duplicates

ndf <- ndf %>%
  mutate(thinned_id = paste(species, TO_GRIDREF, lon, lat))

ndf_thinned <- ndf[!duplicated(ndf$thinned_id),]

#' Use the cpa() function to create presence and absences
#' I do this using parallelisation because it can take a while depending on 
#' the dataset. I don't think it does for this dataset though.
#' 
#' I have set the function to choose 10000 pseudoabsences for each species
#' In reality, because we have subsetted the data to a smaller area, the number
#' of pseudoabsences is set to the maximum number of sightings of other species
#' possible in the dataset. There are warning messaged stating this to be the
#' case.
#' 
#' Also, if the number of presence points > than the number of pseudoabsences 
#' asked for, then the number of pseudoabsences is increased to match the number
#' on presences.
#' 
#' There are year dates to remove data outside of the desired range.
#' screenRaster checks the presence and absence points to make sure none fall
#' outside of the raster layers we are using. This is the thing that takes 
#' longest in the function so you can disable it by setting 
#' screenRaster = NULL


registerDoParallel(cores = detectCores() - 1)

system.time( 
  ab1 <- foreach(i = 1 : length(spp),
                 .packages = c('raster')) %dopar% {
                   
                   
                   cpa(spdat = ndf, species = spp[i], matchPres = FALSE,
                       nAbs = 10000,
                       minYear = 2000, maxYear = 2017, recThresh = 5,
                       screenRaster = ht)
                   
                 }
)

registerDoSEQ()

#' name the list with the original species names
#' this is needed for the next function to work...
names(ab1) <- spp


#' Just going to do this for one species but it would be easy enough
#' to go through each species in turn or for each model in turn.
#' On lotus, all the models are run separately (i.e. in different jobs)
#' and then combined afterwards.
#' 

#' ## choose the model and species

# choose model 
model = "lr" # one of ("lr", "lrReg", "rf", "gam", "me")

# choose species
species = spp[1]
species

# choose number of validations/bootstraps
k = 2

# run the model
sdm <- fsdm(species = species, 
            model = model,
            climDat = ht, 
            spData = ab1, 
            knots_gam = 4, # number of knots to use for GAM, = -1 for default
            k = k, # number of k-fold validations to do
            write = FALSE)


#' ## Predicting code
#' 
#' This is a big chunk of code that behaves differently for each model because
#' of the behaviours of the different predict() functions.
#' It takes a long time to run on my computer. Partly because of the predict()
#' function taking a long time, but also th calc() function trying to get
#' summary values across all of the layers in the raster stack
#' 
#' Having a look at this with fresh eyes has made me realise that I could make 
#' this into a function - so I might do that soon.


# choose the type and index for predict function
if(model == 'lr'|model == 'gam'){
  
  type <- "response"
  index <- NULL
  
} else if(model == 'rf'){
  
  type <- "prob"
  index <- 2
  
}


## bootstrapped models
print(paste0('#####   predicting from bootstrapped models   #####')) 

## predict from each of the bootstrapped models
## different workflow for lrReg and other methods
if(model != 'lrReg') {
  
  # predict from each of the bootstrapped models and stack them together
  boots_out <- raster::stack(lapply(sdm$Bootstrapped_models, FUN = function(x) predict(ht, x, type=type, index=index)))
  
  ## quantiles
  print(paste0('#####   getting quantiles   #####'))
  mean_preds <- calc(boots_out, fun = mean, na.rm = T) # the mean
  quant_preds <- calc(boots_out, fun = function(x) {quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)}) # get the quantile variation
  rnge <- quant_preds[[2]]-quant_preds[[1]] # get the range of max - min
  
  
} else if(model == 'lrReg') { 
  
  print(paste0('#####   predicting for lrReg bootstrapped models   #####'))
  
  ## convert variables to matrix
  covsMat <- as.matrix(rasterToPoints(env_dat)) # convert variables to matrix
  
  ## predict from lrReg model
  boots <- stack(lapply(sdm$Bootstrapped_models, FUN = function(x) {
    
    pred <- predict(x, covsMat[, 3:ncol(covsMat)], type = "response") # predict from matrix
    
    pred <- as.matrix(cbind(covsMat[, 1:2], pred)) # combine predictions with east - norths
    
    if (any(is.na(pred[, 3]))) pred <- pred[-which(is.na(pred[,3])), ] # get rid of NAs
    
    pred_rast <- rasterize(pred[, 1:2], env_dat[[1]], field = pred[, 3]) # turn into a raster of probabilities
    
    ## return predictions
    return(pred_rast)
    
  }))
  
  ## check for intercept-only models which are 0.5 probability in all cells
  uniqueVals <-lapply(1:k,
                      function(x) { length(cellStats(boots[[x]], unique)) })
  
  drop <- which(uniqueVals <= 2) ## i.e. the mean and NA, because with intercept-only models the only values are 0.5 and NA
  
  # assign intercept only models an AUC of NULL - important for weighted average later
  if(any(drop)){
    
    sdm$AUC[drop] <- NA
    
    print(paste("Dropping", length(drop), "intercept-only model(s). Intercept-only models are given an AUC value of NA so they can be identified.
                        Where 1:(k-1) models are intercept only, only the non-intercept models are included in the final average. Where all models are intercept-only,
                        their predictions are returned but should not be used."))
    
  }
  
  
  if (length(drop) == k | length(drop) == 0) {
    
    ## quantiles
    print(paste0('#####   getting quantiles lrReg   #####'))
    mean_preds <- mean(boots) # where all models are intercept-only, takes the mean to avoid errors later but AUC scores are NA which means they are dropped for final ensembles
    quant_preds <- calc(boots, fun = function(x) {quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)})
    rnge <- quant_preds[[2]]-quant_preds[[1]]
    
    
  } else {
    
    print(paste0('#####   getting quantiles lrReg   #####'))
    mean_preds <- mean(boots[[-drop]])
    quant_preds <- calc(boots[[-drop]], fun = function(x) {quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)})
    rnge <- quant_preds[[2]]-quant_preds[[1]]
    
  }
  
  ## if some models were intercept-only then recalculate k (number of models used)
  k <- k - length(drop)
  
}



#' ## Get the decide score
#' For the moment, the decide score for a single species is just the probability
#' of presence * the quantile range. So for a single species and single
#' model, the decide score is:
#' 
decide <- mean_preds*rnge
plot(decide)

