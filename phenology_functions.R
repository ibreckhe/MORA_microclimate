####Functions for the phenology analysis####
####Author: Ian Breckheimer
####April 9, 2014.

# Function to classify the observations into "successes" and 
# "failures" based on the value of a factor or a character vector.
# The function returns a data frame identical to the input, but with
# with an extra column of logical data indicating whether the
# observation met the criteria.

class_obs <- function(data,spec_col,species,stage_col,stage){
  
  # Checks data to make sure both columns exist
  allcols <- colnames(data)
  col_test <- all(spec_col %in% allcols & stage_col %in% allcols)
  if (col_test == FALSE){stop("One of those columns doesn't exist!")}
  # Creates the output vector.
  criterion <- data[,spec_col] %in% species & data[,stage_col] == stage
  
  # Checks to make sure that at least one row is true,
  # otherwise print a warning
  if (any(criterion)==FALSE){warning("No rows met the criteria!")}
  
  # Binds to input.
  output <- cbind(data,criterion)
  return(output)
}

## Function to bin observations in a numeric column and append to a data frame.
bin_obs <- function(data,column,breaks,new_colname){
  
  # Checks data to make sure all the columns exist
  allcols <- colnames(data)
  col_test <- all(column %in% allcols)
  if (col_test == F){stop("That column doesn't exist!")}
  
  # Makes sure the data column is numeric.
  num_test <- is.numeric(data[,c(column)])
  if (num_test == F){stop("All of the data columns must be numeric.")}
  
  # Calculates bin centers for labels.
  centers <- c()
  for (i in 1:(length(breaks)-1)){
    centers[i] <- (breaks[i] + breaks[i+1])/2
  }
  
  # Uses the cut function to categorize the numeric column.
  cutdat <- data.frame(cut(data[,column],breaks=breaks,ordered_result=TRUE,
                           labels=centers))
  names(cutdat) <- new_colname
  
  # Appends that column to the data frame.
  output <- cbind(data,cutdat)
}

## Function to summarize the number of successes and trials data frame 
## by combinations of binned values.
sum_obs <- function(data,bin_cols,success_colname,trial_colname){
  
  # Uses the ddply function in plyr.
  require(plyr)
  
  # Checks data to make sure all the columns exist
  allcols <- colnames(data)
  col_test <- all(bin_cols %in% allcols)
  if (col_test == F){stop("One of those columns doesn't exist!")}
  # Does the summation.
  out <- ddply(data,bin_cols,.fun=summarise,
                  successes = sum(criterion),
                  trials = length(criterion))

  # Renames the columns
  names(out)[names(out)=="successes"] <- success_colname
  names(out)[names(out)=="trials"] <- trial_colname
  return(out)
}

## Function to extract AIC from a list of model objects.
aic_list <- function(model_list){
  aic_all <- c()
  for (i in 1:length(model_list)){
    aic_i <- AIC(model_list[[i]])
    names(aic_i) <- names(model_list)[i]
    aic_all <- c(aic_all,aic_i)
  }
  return(aic_all)
}

# Function to find the values that define 68.4% of the area under the curve.
obs_intervals <- function(preds,bin_width=0.1,threshold=0.158){
  
  # Temporary variables.
  dss <- preds$dss
  pred <- preds$pred
  
  # Assume all NA predictions are zero
  pred[is.na(pred)] <- 0
  
  # Total area under the curve.
  area <- sum(pred*bin_width,na.rm=TRUE)
  
  # Computes cumulative proportions in both directions
  cumprop_up <- cumsum(pred)/area
  cumprop_down <- rev(cumsum(rev(pred)))/area
  
  # Finds the indices of the first and last values greater than 0.158
  lwr_index <- min(which(cumprop_up >= threshold))
  upr_index <- max(which(cumprop_down >= threshold))
  
  # Finds the corresponding values of dss.
  lwr_bound <- dss[lwr_index]
  upr_bound <- dss[upr_index]
  bounds <- c(lwr_bound,upr_bound)
  names(bounds) <- c("lwr_bound","upr_bound")
  
  # Output
  return(bounds)
}

# Function to return 2.5 and 97.5% prediction quantiles from each boostrap sample.
boot_preds <- function(model,data,dss_breaks,pred.data,n_replicates,clustervar,intervals=FALSE,optimums=FALSE){
  # Distributes jobs to nodes.
  require(foreach)
  require(doParallel)
  require(data.table)
  
  cl <- makeCluster(4)
  registerDoParallel(cl)
  spp <- model$data$Species[1]
  
  # Reduces complexity of the data, keeping only columns we need
  data <- data[,c("id","SPECIES","PHEN_PHASE","datetaken_DOY","dss","ownerdate")]     
  
  # Creates a list with clusters.
  # get a vector with all clusters
  clust <- sort(unique(data[,clustervar]))
  
  # group the data points per cluster
  clust.group <- function(clust) {
    data[data[,clustervar]==clust,]
  }
  
  clust.list <- lapply(clust,clust.group)
  
  boot.pred.df <- foreach(i=1:n_replicates,.combine=rbind) %dopar% {
    # parallel processes can't see these objects in the global env.
    quad_form <- formula(prop_flower~dss_num+I(dss_num^2))
    log_form <- formula(prop_flower~log(dss_num)+I(log(dss_num)^2))
    source("~/code/MORA_microclimate/phenology_functions.R")
    require(data.table)
    
    # Resamples the cluster list and reassembles it as a data frame.
    boot.ind <- sample(clust.list,length(clust.list),replace=TRUE)
    boot.resamp <- as.data.frame(rbindlist(boot.ind)) #speedy list to data frame
    
    # Re-bins the data.
    boot.classed <- class_obs(boot.resamp,"SPECIES",spp,"PHEN_PHASE","FLOWERING")
    boot.binned <- bin_obs(boot.classed,column="dss",breaks=dss_breaks,
                      new_colname="dss_bin")
    boot.data <- sum_obs(boot.binned,"dss_bin",success_colname=spp,
                      trial_colname="All_Photos")
    boot.data <- boot.data[-c(dim(boot.data)[1]),] #Drops observations outside of the bin range.
    boot.data$Species <- spp
    boot.data$dss_num <- as.numeric(as.character(boot.data$dss_bin))
    colnames(boot.data) <- c("dss_bin","Flower_Phot","All_Phot","Species","dss_num")
    boot.data$prop_flower <- boot.data$Flower_Phot / boot.data$All_Phot
    
    #Updates the model
    boot.model <- update(model,data=subset(boot.data,Species==spp))
    boot.pred <- predict(boot.model,newdata=data.frame(dss_num=pred.data),type="response")
    boot.pred.df <- data.frame(Species=model$data$Species[1],dss=pred.data,boot.pred)
  }
  
  # Summarize output using dplyr function
  bnds <- group_by(boot.pred.df,"dss") %.%
          summarise(Species = Species[1],
                    lwr_bnd = quantile(boot.pred,0.025,na.rm=TRUE),
                    upr_bnd = quantile(boot.pred,0.975,na.rm=TRUE))
  
  return(bnds)
}

# Function to return 2.5 and 97.5% prediction quantiles from each boostrap sample.
boot_preds2 <- function(model,data,dss_breaks,pred.data,n_replicates,
                        clustervar,colnames,species_col,stage_col,stage){
  # Distributes jobs to nodes.
  require(foreach)
  require(plyr)
  require(doParallel)
  require(data.table)
  
  cl <- makeCluster(4)
  registerDoParallel(cl)
  spp <- model$data$Species[1]
  
  # Reduces complexity of the data, keeping only columns we need
  data <- data[,colnames]     
  
  # Classifies observations
  boot.classed <- class_obs(data,species_col,spp,stage_col,stage)
  
  ## Removes other records of the particular photos with target species
  photos <- unique(boot.classed$photo[boot.classed$criterion==TRUE])
  boot.classed2 <- boot.classed[-c(which((boot.classed$photo %in% photos & boot.classed$criterion==FALSE))),]
  
  # Creates a list with clusters.
  # get a vector with all clusters
  clust <- sort(unique(boot.classed2[,clustervar]))
  
  # group the data points per cluster
  clust.group <- function(clust) {
    boot.classed2[boot.classed2[,clustervar]==clust,]
  }
  
  clust.list <- lapply(clust,clust.group)
  
  boot.pred.df <- foreach(i=1:n_replicates,.combine=rbind) %dopar% {
    # parallel processes can't see these objects in the global env.
    quad_form <- formula(prop_flower~dss_num+I(dss_num^2))
    log_form <- formula(prop_flower~log(dss_num)+I(log(dss_num)^2))
    source("~/code/MORA_microclimate/phenology_functions.R")
    require(data.table)
    
    # Resamples the cluster list
    boot.clust <- sample(clust.list,length(clust.list),replace=TRUE)
    
    # Resamples all the photos in that cluster
    boot.fun <- function(x){x[sample(1:(dim(x)[1]),replace=TRUE),]}
    boot.clust2 <- lapply(boot.clust,boot.fun)
    
    # Assembles everything back into a data frame
    boot.resamp <- as.data.frame(rbindlist(boot.clust2)) #speedy list to data frame
    
    # Re-bins the data.
    boot.binned <- bin_obs(boot.resamp,column="dss",breaks=dss_breaks,
                           new_colname="dss_bin")
    boot.data <- sum_obs(boot.binned,"dss_bin",success_colname=spp,
                         trial_colname="All_Photos")
    boot.data <- boot.data[-c(dim(boot.data)[1]),] #Drops observations outside of the bin range.
    boot.data$Species <- spp
    boot.data$dss_num <- as.numeric(as.character(boot.data$dss_bin))
    colnames(boot.data) <- c("dss_bin","Flower_Phot","All_Phot","Species","dss_num")
    boot.data$prop_flower <- boot.data$Flower_Phot / boot.data$All_Phot
    
    #Updates the model
    boot.model <- update(model,data=subset(boot.data,Species==spp))
    boot.pred <- predict(boot.model,newdata=data.frame(dss_num=pred.data),type="response")
    boot.pred.df <- data.frame(rep=i,Species=model$data$Species[1],dss=pred.data,pred=boot.pred)
  }

  # Gets rid of NaN predictions.
  boot.pred.df$pred[is.nan(boot.pred.df$pred)] <- NA

  # Summarize output using the ddply function
  bnds <- ddply(boot.pred.df,~dss,.fun=summarise,
              Species = Species[1],
              lwr_bnd = quantile(pred,0.025,na.rm=TRUE),
              upr_bnd = quantile(pred,0.975,na.rm=TRUE))
  
  opt <- ddply(boot.pred.df,~rep,.fun=summarise,
               max=dss[which.max(pred)])
  opt_lwr <- quantile(opt$max,0.025,na.rm=TRUE)
  opt_upr <- quantile(opt$max,0.975,na.rm=TRUE)
  opt_int <- c(opt_lwr,opt_upr)
  
  intervals <- ddply(boot.pred.df,~rep,.fun=obs_intervals,
                     threshold=0.158,bin_width=1)
  intervals$tol <- (intervals$upr_bound - intervals$lwr_bound) / 2
  tol_lwr <- quantile(intervals$tol,0.025,na.rm=TRUE)
  tol_upr <- quantile(intervals$tol,0.975,na.rm=TRUE)
  tol_int <- c(tol_lwr,tol_upr)
  
  return(list(prediction_intervals=bnds,
              optimum_intervals=opt_int,
              tolerance_intervals=tol_int))
}