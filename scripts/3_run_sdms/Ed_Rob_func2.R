fsdm <- function(species, model, climDat, spData, k, write, outPath, inters = F, knots = -1){#, filename){
  
  # species = spp[s]
  # model = "rf"
  # climDat = ht
  # spData = ab1
  # knots = -1
  # k = 5
  # boot = 1000
  
  ind <- which(names(spData) == species)
  spDat <- spData[[ind]]
  if (is.null(spDat)) {
    out <- NULL
  }
  else {
    pres <- data.frame(val = 1, raster::extract(x = climDat, y = spDat$Presence, na.rm = T), spDat$Presence)
    # pres$land_cov <- as.factor(pres$land_cov)
    if (any(is.na(pres))) {
      pres <- na.omit(pres)
    }
    nRec <- nrow(pres)
    ab <- data.frame(val = 0, raster::extract(x = climDat, y = spDat$pseudoAbsence, na.rm = T), spDat$pseudoAbsence)
    # ab$land_cov <- as.factor(ab$land_cov)
    if (any(is.na(ab))) {
      ab <- na.omit(ab)
    }
    allDat <- rbind(pres[!names(pres) %in% c("lon", "lat")], ab[!names(ab) %in% c("lon", "lat")])
    allDat_loc <- rbind(pres, ab) # have a data frame with lat lon at the end to output
    if (model == "lr") {
      
      # no interactions
      if(inters == F) {
        fullMod <- glm(val ~ ., data = allDat, family = binomial(link = "logit"))
      }
      
      # .*. is how to include all two-way interactions
      if(inters == T) {
        fullMod <- glm(val ~ .*., data = allDat, family = binomial(link = "logit"))
      }
      
      type <- "response"
      index <- NULL
    }
    
    else if (model == "rf") {
      
      fullMod <- ranger(x = allDat[!names(allDat) %in% c("val")], 
                        y = as.factor(allDat$val),
                        importance = "impurity",
                        probability = TRUE,
                        replace = TRUE,
                        num.threads = 1) 
      type <- "response"
      index <- 2
    }  else if(model == "me"){
      fullMod <- maxent(x = climDat, p = data.frame(spDat$Presence)[,1:2], a = data.frame(spDat$pseudoAbsence)[,1:2])
    } 
    else if(model == "gam"){
      
      # gams need variables with ~ >10 knots to fit them automatically
      # so need to remove the variables that appear infrequently in the dataset
      # get a dataframe with the number of unique values for each variable
      l <- sapply(allDat, unique)
      ks <- rownames_to_column(data.frame(k = round(sapply(l, length))[-1]), 
                               var = "variable")
      
      "%!in%" <- Negate("%in%")
      # drop variables accordiong to number of knots asked for
      # -1 is basically 9 knots
      if(knots == -1) {
        v_keep <- ks[ks$k > 11,]
        
        print(paste("variable dropped =", ks$variable[ks$variable %!in% v_keep$variable]))
      }
      
      # any others just keep the variables with over the number of knots
      if(knots > 0) {
        v_keep <- ks[ks$k > knots,]
        print(paste("variable dropped =", ks$variable[ks$variable %!in% v_keep$variable]))
      }
      v_keep
      
      form <- as.formula(paste0("val ~ s(", paste(v_keep$variable, 
                                                  ", k = ", knots) %>% #, 
                                  # v_keep$knts) %>% 
                                  paste0(collapse = ") + s("), ")"))
      form
      
      
      system.time(
        fullMod <- gam(formula = form, data = allDat, 
                       family = binomial(link = 'logit'), 
                       select = T, method = 'REML', gamma = 1.4))
      
      type <- "response"
      index <- NULL
    }
    
    ###### PREDICTION
    # separate predict function for random forests from ranger
    if(model == 'rf'){
      pred <- predict(climDat, fullMod, type = type, index = index,
                      predict.all = F,
                      fun = function(model, ...) predict(model, ...)$predictions,
                      num.threads = 1)
    }
    if(model == 'lr'  | model == 'gam' | model == 'me'){
      pred <- predict(climDat, fullMod, type = type, index = index)
    }
    
    # raster::plot(pred, col = matlab.like(30))
    # points(spDat$Presence, pch = "+", cex = 0.4)
    
    folds <- c(kfold(pres, k), kfold(ab, k))
    folds_me_pres <- kfold(spDat$Presence, k)
    folds_me_ab <- kfold(spDat$pseudoAbsence, k)
    e <- list()
    rf_int <- list()
    for (i in 1:k) {
      print(i)
      
      if (model == 'me') {
        
        train_me_pres <- spDat$Presence[folds_me_pres != i, ]
        train_me_abs <- spDat$pseudoAbsence[folds_me_ab != i, ]
        
        test_me_pres <- spDat$Presence[folds_me_ab == i, ]
        test_me_abs <- spDat$pseudoAbsence[folds_me_ab == i, ]
        
        mod <- maxent(x = climDat, p = data.frame(train_me_pres)[,1:2], a = data.frame(train_me_abs)[,1:2])
        
        e[[i]] <- dismo::evaluate(p = test_me_pres, a = test_me_abs, x = climDat,
                                  mod, tr = seq(0, 1, length.out = 200))
        
      }
      
      if(model == 'lr' | model == 'rf' | model == 'gam'){
        
        train <- allDat[folds != i, ]
        test <- allDat[folds == i, ]
        
        if (model == "lr") {
          mod <- glm(val ~ ., data = train, family = binomial(link = "logit"))
        }
        else if (model == "rf") {
          
          mod <- ranger(x = train[!names(train) %in% c("val")], # want to keep latlon in the data frame for later so choosing all columns except for these three
                        y = as.factor(train[, 1]),
                        importance = "impurity",
                        probability = TRUE,
                        replace = TRUE,
                        num.threads = 1)
          
          # rf_int[[i]] <- rfinterval(formula = val~., train_data = train, test_data = test, alpha = 0.05)
          
        }
        else if (model == "gam"){
          mod <- gam(formula = form, data = allDat, 
                     family = binomial(link = 'logit'), 
                     select = T, method = 'REML', gamma = 1.4)
        }
        
        if(model != 'rf'){
          e[[i]] <- dismo::evaluate(p = test[test$val == 1, ], 
                                    a = test[test$val == 0, ], model = mod, 
                                    tr = seq(0, 1, length.out = 200))
        }
        else if(model == 'rf'){
          
          # needs to be different because of weird predict function behaviour with
          # the predict.ranger function
          train <- allDat_loc[folds != i, ]
          test <- allDat_loc[folds == i, ]
          
          e[[i]] <- dismo::evaluate(predict(mod, data = test[test$val == 1, ],
                                            num.threads = 1)$predictions[ , 2, drop = TRUE],
                                    predict(mod, data = test[test$val == 0, ],
                                            num.threads = 1)$predictions[ , 2, drop = TRUE], 
                                    tr = seq(0, 1, length.out = 200))
          
        }
        
        
      }
      
    }
    
    
    # ###########
    # ####    Bootstrap random forest to get Confidence Intervals
    # ####    So am going to randomly sample from both presence and absence points
    # ####    This is different problem to trying to make the RF model as accurate as possible
    # ####    by subsampling absence points
    # if(model == 'rf'){
    #   boot_out <- list()
    #   
    #   for(b in 1:boot){
    #     
    #     # use presence and absence DFs from above
    #     # sample from presences and absences, match the number of presences with absences
    #     t_boot_p <- pres[sample(nrow(pres), nrow(pres), replace = T), ]
    #     t_boot_a <- ab[sample(nrow(ab), nrow(pres), replace = T), ]
    #     
    #     boot_df <- rbind(t_boot_p, t_boot_a)
    #     
    #     boot_pred <- ranger(x = boot_df[!names(boot_df) %in% c("val", "lon", "lat")], # want to keep latlon in the data frame for later so choosing all columns except for these three
    #                   y = as.factor(boot_df$val),
    #                   importance = "impurity",
    #                   probability = TRUE,
    #                   replace = TRUE)$predictions[,1]
    #     
    #     # type <- "response"
    #     # index <- 2
    #     # pred <- predict(climDat, boot_pred, type = type, index = index,
    #     #                 predict.all = F,
    #     #                 fun = function(model, ...) predict(model, ...)$predictions)
    #     
    #     boot_out[[b]] <- boot_pred
    #     
    #   }
    #   
    #   t <- data.frame(do.call("cbind", boot_out))
    #   head(t)
    #   
    #   lwr <- apply(t, 1, function(x) quantile(x, probs=0.05))
    #   upr <- apply(t, 1, function(x) quantile(x, probs=0.95))
    #   pred_ci <- predict(fullMod, data = rbind(pres,ab), type = type, index = index,
    #                   predict.all = F)$predictions[,1]
    #   # still wrong, predicted values and bootstrapped values have different lengths...
    #   b_o <- data.frame(cbind(pred_ci, lwr, upr))
    #   
    #   b_o$ind <- seq(1, length(b_o$lwr))
    #   
    #   ggplot(data = b_o) +
    #     geom_errorbar(aes(x = reorder(ind, pred_ci), ymin = lwr, ymax = upr), colour = "red", size = 0.1) +
    #     geom_point(aes(x = reorder(ind, pred_ci), y = pred_ci), size = 0.7)
    #     
    #   
    #   bs <- data.frame(apply(t, 1, quantile, probs = c(0.05, 0.95),  na.rm = TRUE))
    #   head(bs)
    #   bootstrap <- do.call('rbind', lapply(simplify2array(boot_out), 2 , FUN = quantile, probs = c(0.05, 0.95)))
    #   colnames(bootstrap) <- c("lwr", "upr")
    #   head(bootstrap)
    #   bs <- rowname
    #   ggplot(bootstrap, aes(x = ))
    # }
    
    
    auc <- mean(sapply(e, function(x) {
      slot(x, "auc")
    }))
    out <- NULL
    out <- list(species, nRec, fullMod, auc, k, pred, allDat_loc, e)
    names(out) <- c("Species", "Number of records", "Model", 
                    "AUC", "Number of folds for validation", "Predictions", 
                    "Data", "Model_evaluation")
    
    if(model == "rf") {
      out <- list(species, nRec, fullMod, auc, k, pred, allDat_loc, e,  rf_int)
      names(out) <- c("Species", "Number of records", "Model", 
                      "AUC", "Number of folds for validation", "Predictions", 
                      "Data", "Model_evaluation", "RF_interval")
    }
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
                 recThresh, replace = F) 
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
    }
    ab <- ab[sampInd, ]
    out <- list(pres, ab)
    names(out) <- c("Presence", "pseudoAbsence")
  }
  return(out)
}


# presence absence on National grid
cpa_grid <- function (spdat, species, minYear, maxYear, nAbs, matchPres = FALSE, 
                      recThresh, replace = F) 
{
  dat <- spdat[spdat$year >= minYear & spdat$year <= maxYear,]
  pres <- dat[dat$species == species, c("EASTING", "NORTHING")]
  if (nrow(pres) < recThresh) {
    warning("Number of records does not exceed recThresh")
    out <- NULL
  }
  else {
    ab <- dat[dat$species != species, c("EASTING", "NORTHING")]
    "%!in%" <- Negate("%in%")
    ab <- ab[ab %!in% pres]
    if (nrow(ab) < nrow(pres)) {
      warning("More presences than possible locations for absences. Consider lowering the number of pseudo absences.")
    }
    sampInd <- sample(1:nrow(ab), nAbs, replace = replace)
    if (matchPres == TRUE) {
      sampInd <- sampInd[1:nrow(pres)]
    }
    ab <- ab[sampInd, ]
    out <- list(pres, ab)
    names(out) <- c("Presence", "pseudoAbsence")
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


# Move points
move.points <- function(r, pts, spatial=FALSE) {
  require(raster)
  require(sp)
  
  if (is(pts, 'SpatialPoints')) pts <- coordinates(pts)
  if (is(!r, 'Raster')) r <- raster(r)
  
  loc <- colSums(sapply(pts[, 1], '>', bbox(r)[1, ])) * 3 + 
    colSums(sapply(pts[, 2], '>', bbox(r)[2, ]))
  
  L <- split(as.data.frame(pts), loc)
  
  new.pts <- lapply(names(L), function(x) {
    switch(x, 
           '0' = xyFromCell(r, ncell(r) - ncol(r) + 1)[rep(1, nrow(L[[x]])), ],
           '1' = xyFromCell(r, cellFromXY(r, cbind(xmin(r), L[[x]][, 2]))),
           '2' = xyFromCell(r, 1)[rep(1, nrow(L[[x]])), ],
           '3' = xyFromCell(r, cellFromXY(r, cbind(L[[x]][, 1], ymin(r)))),
           '4' = {
             xy <- as.matrix(L[[x]])
             dimnames(xy) <- list(NULL, c('x', 'y'))
             xy
           },
           '5' = xyFromCell(r, cellFromXY(r, cbind(L[[x]][, 1], ymax(r)))),
           '6' = xyFromCell(r, ncell(r))[rep(1, nrow(L[[x]])), ],
           '7' = xyFromCell(r, cellFromXY(r, cbind(xmax(r), L[[x]][, 2]))),
           '8' = xyFromCell(r, ncol(r))[rep(1, nrow(L[[x]])), ]
    )
  })
  
  new.pts <- unsplit(mapply(function(x, y) {
    row.names(x) <- row.names(y)
    as.data.frame(x)
  }, new.pts, L, SIMPLIFY=FALSE), loc)
  
  colnames(new.pts) <- colnames(pts)
  if(isTRUE(spatial)) new.pts <- SpatialPoints(new.pts)  
  return(new.pts)
}





