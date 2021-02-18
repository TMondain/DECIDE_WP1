
###' Fit SMDs - main function for fitting SDMs to data from create presence absence
###'

fsdm <- function(species, model, climDat, spData, k, write, outPath, #inters = F, prediction = TRUE,
                 knots = -1){ 
  
  # select species
  ind <- which(names(spData) == species)
  spDat <- spData[[ind]]
  
  if (is.null(spDat)) {
    out <- NULL
  }
  else {
    
    #### extract data ####
    pres <- data.frame(val = 1, raster::extract(x = climDat, y = spDat$Presence, na.rm = T), spDat$Presence)
    
    if (any(is.na(pres))) { # this might not be needed with the new screenraster argument for the create pseudoabs function
      pres <- na.omit(pres)
    }
    
    nRec <- nrow(pres)
    
    print(paste("Occurrence records:", nRec))
    
    ab <- data.frame(val = 0, raster::extract(x = climDat, y = spDat$pseudoAbsence, na.rm = T), spDat$pseudoAbsence)
    
    if (any(is.na(ab))) {
      ab <- na.omit(ab)
    }
    
    
    allDat <- rbind(pres[!names(pres) %in% c("lon", "lat")], ab[!names(ab) %in% c("lon", "lat")])
    allDat_loc <- rbind(pres, ab)
    
    # if (model == "lr") {
    #   
    #   fullMod <- glm(val ~ ., data = allDat, family = binomial(link = "logit"))
    #   
    #   type <- "response"
    #   index <- NULL
    # }
    # 
    # else if (model == "rf") {
    #   fullMod <- randomForest(x = allDat[, 2:ncol(allDat)], 
    #                           y = as.factor(allDat[, 1]), importance = T, 
    #                           norm.votes = TRUE)
    #   type <- "prob"
    #   index <- 2
    # }  else if(model == "me"){
    #   fullMod <- maxent(x = climDat, p = data.frame(spDat$Presence)[,1:2], a = data.frame(spDat$pseudoAbsence)[,1:2])
    # } 
    # else if(model == "gam"){
    #   
    #   # gams need variables with ~ >10 knots to fit them automatically
    #   # so need to remove the variables that appear infrequently in the dataset
    #   # get a dataframe with the number of unique values for each variable
    #   l <- sapply(allDat, unique)
    #   ks <- rownames_to_column(data.frame(k = round(sapply(l, length))[-1]), 
    #                            var = "variable")
    #   
    #   "%!in%" <- Negate("%in%")
    #   # drop variables accordiong to number of knots asked for
    #   # -1 is basically 9 knots
    #   if(knots == -1) {
    #     v_keep <- ks[ks$k > 11,]
    #     
    #     print(paste("variable dropped =", ks$variable[ks$variable %!in% v_keep$variable]))
    #   }
    #   
    #   # any others just keep the variables with over the number of knots
    #   if(knots > 0) {
    #     v_keep <- ks[ks$k > (knots+3),]
    #     print(paste("variable dropped =", ks$variable[ks$variable %!in% v_keep$variable]))
    #   }
    #   v_keep
    #   
    #   form <- as.formula(paste0("val ~ s(", paste(v_keep$variable, 
    #                                               ", k = ", knots) %>% #, 
    #                               # v_keep$knts) %>% 
    #                               paste0(collapse = ") + s("), ")"))
    #   form
    #   
    #   
    #   system.time(
    #     fullMod <- gam(formula = form, data = allDat, 
    #                    family = binomial(link = 'logit'), 
    #                    select = T, method = 'REML', gamma = 1.4))
    #   
    #   type <- "response"
    #   index <- NULL
    # }
    # 
    # if(prediction == TRUE){
    #   pred <- predict(climDat, fullMod, type = type, index = index)
    # } else if(prediction == FALSE){ 
    #   pred <- NULL 
    # }
    
    # raster::plot(pred, col = matlab.like(30))
    # points(spDat$Presence, pch = "+", cex = 0.4)
    
    ###### Move straight to 'bootstrapping' ######
    
    folds <- c(kfold(pres, k), kfold(ab, k))
    folds_me_pres <- kfold(spDat$Presence, k)
    folds_me_ab <- kfold(spDat$pseudoAbsence, k)
    
    e <- vector('list', length = k)
    mods_out <- vector('list', length = k)
    
    for (i in 1:k) {
      
      
      if (model == 'me') {
        
        train_me_pres <- spDat$Presence[folds_me_pres != i, ]
        train_me_abs <- spDat$pseudoAbsence[folds_me_ab != i, ]
        
        test_me_pres <- spDat$Presence[folds_me_ab == i, ]
        test_me_abs <- spDat$pseudoAbsence[folds_me_ab == i, ]
        
        mod <- maxent(x = climDat, p = data.frame(train_me_pres)[,1:2], 
                      a = data.frame(train_me_abs)[,1:2])
        
        e[[i]] <- dismo::evaluate(p = test_me_pres, 
                                  a = test_me_abs, 
                                  x = climDat,
                                  mod, tr = seq(0, 1, length.out = 200))
        
      }
      
      if(model == 'lr' | model == 'rf' | model == 'gam'){
        
        train <- allDat[folds != i, ]
        test <- allDat[folds == i, ]
        
        if (model == "lr") {
          mod <- glm(val ~ ., data = train, family = binomial(link = "logit"))
        }
        else if (model == "rf") {
          mod <- randomForest(x = train[, 2:ncol(train)], 
                              y = as.factor(train[, 1]), 
                              importance = T, 
                              norm.votes = TRUE)
          
        }
        else if (model == "gam"){
          
          ## create formula for gam
          l <- sapply(allDat, unique)
          ks <- rownames_to_column(data.frame(k = round(sapply(l, length))[-1]),
                                   var = "variable")
          
          "%!in%" <- Negate("%in%")
          # drop variables according to number of knots asked for
          # -1 is basically 9 knots
          if(knots == -1) {
            v_keep <- ks[ks$k > 11,]
            
            print(paste("variable dropped =", ks$variable[ks$variable %!in% v_keep$variable]))
          }
          
          # any others just keep the variables with over the number of knots
          if(knots > 0) {
            v_keep <- ks[ks$k > (knots+3),]
            print(paste("variable dropped =", ks$variable[ks$variable %!in% v_keep$variable]))
          }
          # v_keep
          
          form <- as.formula(paste0("val ~ s(", paste(v_keep$variable,
                                                      ", k = ", knots) %>% #,
                                      # v_keep$knts) %>%
                                      paste0(collapse = ") + s("), ")"))
          # form
          
          
          mod <- gam(formula = form, data = train, 
                     family = binomial(link = 'logit'), 
                     select = TRUE, method = 'REML', gamma = 1.4)
        }
        
        e[[i]] <- evaluate(p = test[test$val == 1, ], 
                           a = test[test$val == 0, ], 
                           mod, tr = seq(0, 1, length.out = 200))
        
        
      }
      
      # store all the bootstrapped models
      mods_out[[i]] <- mod
      
    }
    
    auc <- mean(sapply(e, function(x) {
      slot(x, "auc")
    }))
    
    ## old storage when predicting from full model
    # out <- NULL
    # out <- list(species, nRec, fullMod, auc, k, pred, allDat_loc, e, mods_out)
    # names(out) <- c("Species", "Number of records", "Model", 
    #                 "AUC", "Number of folds for validation", "Predictions", 
    #                 "Data", "Model_evaluation", "Bootstrapped_models")
    
    out <- NULL
    out <- list(species, nRec, 
                auc, k, allDat_loc, e, mods_out)
    names(out) <- c("Species", "Number of records", 
                    "AUC", "Number of folds for validation",
                    "Data", "Model_evaluation", "Bootstrapped_models")
    
    if (write == TRUE) {
      print(species)
      save(out, file = paste0(outPath, species, "_", model, 
                              ".rdata"))
    }
  }
  return(out)
}



#### Presence absence function
cpa <- function (spdat, species, minYear, maxYear, nAbs, matchPres = FALSE, 
                 recThresh, replace = F, screenRaster = NULL) 
{
  dat <- spdat[spdat$year >= minYear & spdat$year <= maxYear,]
  pres <- dat[dat$species == species, c("lon", "lat")]
  
  if (nrow(pres) < recThresh) {
    warning("Number of records does not exceed recThresh")
    out <- NULL
  }
  
  else {
    
    ab <- dat[dat$species != species, c("lon", "lat")]
    
    "%!in%" <- Negate("%in%")
    ab <- ab[ab %!in% pres]
    
    if (nrow(ab) < nrow(pres)) {
      warning("More presences than possible locations for absences. Consider lowering the number of pseudo absences.")
    }
    
    sampInd <- sample(1:nrow(ab), nAbs, replace = replace)
    if (matchPres == TRUE) {
      sampInd <- sampInd[1:nrow(pres)]
    } else {
      if (nAbs <= nrow(ab)) {
        sampInd <- sample(1:nrow(ab), nAbs)
      } else {
        warning(paste0("Fewer than 10,000 locations available for pseudo absences when using the target group approach. Setting nAbs to the maximum number possible (", nrow(ab), ")."))
        sampInd <- 1:nrow(ab)
      }
    }
    ab <- ab[sampInd, ]
    out <- list(pres, ab)
    names(out) <- c("Presence", "pseudoAbsence")
    
    ## if screenRaster is specified, check if any presence or absence points fall outside of the raster extent (i.e. they are NA).
    ## If some data fall outside of the extent of the covariates, drop them, and drop the equivalent number of absences or presences
    ## to ensure they are equal in number.
    
    if (!is.null(screenRaster)) {
      
      presDrop <- raster::extract(screenRaster, out$Presence)
      
      abDrop <- raster::extract(screenRaster, out$pseudoAbsence)
      
      if (any(is.na(presDrop))) out$Presence <- out$Presence[-which(is.na(presDrop)), ]
      
      if (any(is.na(abDrop))) out$pseudoAbsence <- out$pseudoAbsence[-which(is.na(abDrop)), ]
      
      if (nrow(out$Presence) > nrow(out$pseudoAbsence)) {
        
        out$Presence <- out$Presence[1:nrow(out$pseudoAbsence), ]
        
      } else if (nrow(out$Presence) < nrow(out$pseudoAbsence)) {
        
        out$pseudoAbsence <- out$pseudoAbsence[1:nrow(out$Presence), ]
        
      }
      
    }
    
  }
  
  return(out)
  
}



#### Convert grid reference from BRC data to eastings and northings
c_en <- function(occ){
  
  coords <- gr_let2num(gridref = occ$TO_GRIDREF, centre = T)
  
  occ$cn <- gr_det_country(occ$TO_GRIDREF)
  
  occ <- cbind(occ, coords)
  
  if (length(unique(occ$cn) > 1)) {
    
    GBCRS <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs")
    
    NICRS <- CRS("+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=1 +x_0=200000 +y_0=250000 +ellps=airy +towgs84=482.5,-130.6,564.6,-1.042,-0.214,-0.631,8.15 +units=m +no_defs")
    
    NIocc <- occ[occ$cn == "OSNI",]
    
    GBocc <- occ[occ$cn == "OSGB",]
    
    NIcoords <- NIocc[, c("EASTING", "NORTHING")]
    
    coordinates(NIcoords) <- c("EASTING", "NORTHING")
    
    proj4string(NIcoords) <- NICRS
    
    NIcoords <- spTransform(NIcoords, GBCRS)
    
    NIocc[,c("EASTING", "NORTHING")] <- data.frame(NIcoords)
    
    occ <- rbind(GBocc, NIocc)
    
  }
  occ
}

### Thin records from Richard Hassall
thin_records = function(occs, templateGrid){
  require(raster)
  require(sp)
  
  # set where to store temporary raster files...
  oldTempDir = rasterOptions(setfileext=T)$tmpdir # include a default option to stop it printing information
  myTempDir = paste(oldTempDir, "thin_records", sep="/")
  assign("rasterOptions", rasterOptions(tmpdir=myTempDir), envir = .GlobalEnv)
  
  # convert records to presence on the template grid...
  rownames(occs) = 1:nrow(occs)
  R = rasterize(occs[,1:2], templateGrid, field=1, max)
  xy = coordinates(R)[!is.na(getValues(R)),] # presence grid cell coords
  
  # reset the temporary raster files and delete any created in the function
  assign("rasterOptions", rasterOptions(tmpdir=oldTempDir), envir = .GlobalEnv)
  unlink(myTempDir, recursive=TRUE)
  
  # return SpatialPoints of the presence grid cell coords...
  SpatialPoints(xy, proj4string=CRS(projection(templateGrid)))
  
}


### Get AUCs for different models
getAUC <- function(file) { # list of different SDMs
  
  f <- list()
  for(i in 1:length(file)){
    df <- file[[i]]
    f[[i]] <- data.frame(species = df$Species, 
                         auc = df$AUC, model = class(df$Model)[1])
  }
  
  f <- do.call("rbind", f)
  return(f)
}


### Average models with AUC weighting
sdm_av <- function (sdm_out, skillDat, species){
  
  rfSkill <- skillDat$auc[skillDat$species == species & skillDat$model == 
                            "randomForest"]
  lrSkill <- skillDat$auc[skillDat$species == species & skillDat$model == 
                            "glm"]
  meSkill <- skillDat$auc[skillDat$species == species & skillDat$model == 
                            "MaxEnt"]
  
  if (rfSkill >= 0.7 | lrSkill >= 0.7| meSkill >= 0.7) {
    
    lr_m <- sdm_out[[1]]
    rf_m <- sdm_out[[2]]
    me_m <- sdm_out[[3]]
    
    if(class(lr_m$Model)[1] != "glm" | 
       class(rf_m$Model)[1] != "randomForest" | 
       class(me_m$Model)[1] != "MaxEnt") {
      
      stop("list not in right order. Order must be -> glm, ranfor, maxent")
      
    }
    
    lrRast <- lr_m$Predictions
    rfRast <- rf_m$Predictions
    meRast <- me_m$Predictions
    
    stack <- stack(rfRast, lrRast, meRast)
    print(c(rfSkill, lrSkill, meSkill))
    if (rfSkill < 0.7) {
      rfSkill <- 0
    }
    if (lrSkill < 0.7) {
      lrSkill <- 0
    }
    if (meSkill < 0.7) {
      meSkill <- 0
    }
    print(c(rfSkill, lrSkill, meSkill))
    ensemble <- weighted.mean(x = stack, w = c(rfSkill, lrSkill, meSkill))
    ensemble
  }
}


# function to calculate partial dependence for a single predictor from eBird tutorial
calculate_pd <- function(predictor, model, data, 
                         x_res = 25, n = 1000) {
  # create prediction grid
  rng <- range(data[[predictor]], na.rm = TRUE)
  x_grid <- seq(rng[1], rng[2], length.out = x_res)
  grid <- data.frame(covariate = predictor, x = x_grid, 
                     stringsAsFactors = FALSE)
  names(grid) <- c("covariate", predictor)
  
  # subsample training data
  n <- min(n, nrow(data))
  s <- sample(seq.int(nrow(data)), size = n, replace = FALSE)
  data <- data[s, ]
  
  # drop focal predictor from data
  data <- data[names(data) != predictor]
  grid <- merge(grid, data, all = TRUE)
  
  # predict
  p <- predict(model, data = grid)
  
  # summarize
  pd <- grid[, c("covariate", predictor)]
  names(pd) <- c("covariate", "x")
  pd$pred <- p$predictions[, 2]
  pd <- dplyr::group_by(pd, covariate, x) %>% 
    dplyr::summarise(pred = mean(pred, na.rm = TRUE)) %>% 
    dplyr::ungroup()
  
  return(pd)
}