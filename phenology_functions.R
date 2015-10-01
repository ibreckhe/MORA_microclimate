####Functions for the phenology analysis####
####Author: Ian Breckheimer
####April 9, 2014.

##Logit and antilogit function.####
antilogit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}


# Function to classify the observations into "successes" and 
# "failures" based on the value of a factor or a character vector.
# The function returns a data frame identical to the input, but with
# with an extra column of logical data indicating whether the
# observation met the criteria.####

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

## Function to bin observations in a numeric column and append to a data frame.####
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

#### Function to summarize the number of successes and trials data frame 
## by combinations of binned values.####
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

# Function to find the values that define 68.4% of the area under the curve.####
obs_intervals <- function(preds,threshold=0.1586553){
  
  # Temporary variables.
  dss <- preds$dss
  pred <- preds$pred
  
  # Gets the bin width from the first pred interval.
  bin_width <- dss[2] - dss[1]
  
  # Assume all NA predictions are zero
  pred[is.na(pred)] <- 0
  
  # Total area under the curve.
  total_area <- sum(pred*bin_width,na.rm=TRUE)
  
  # Computes cumulative proportions in both directions
  cumprop_up <- cumsum(pred)*bin_width/total_area
  cumprop_down <- rev(cumsum(rev(pred))*bin_width)/total_area
  
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

# Function to return 2.5 and 97.5% prediction quantiles from each boostrap sample.####
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

# Function to return 2.5 and 97.5% prediction quantiles from each boostrap sample.####
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

####Function to return optimum, bounds, and range from a fit lme4 merMod object.####
glmer_params <- function(model){
  log <- names(fixef(model))[2]!="Spdss"
  if(log){
    newdata <- expand.grid(log_Spdss=seq(from=-2.5, to=2.5, by=0.01))
  }else{
    newdata <- expand.grid(Spdss=seq(from=-2.5, to=2.5, by=0.01))}
  pred <- predict(object=model,newdata=newdata,re.form=~0,type="response")
  if(log){
    preds <- data.frame(dss=newdata$log_Spdss,pred=pred)
    opt <- newdata$log_Spdss[which.max(pred)]
  }else{
    preds <- data.frame(dss=newdata$Spdss,pred=pred)
    opt <- newdata$Spdss[which.max(pred)]
  }
  range <- obs_intervals(preds,threshold=0.1586553)
  out <- c(opt=opt,range[1],range[2])
  return(out)
}

####Function to fit a Bayesian mixed-effects model in JAGS.####
fit.jags.mixed <- function(x,y,sdd,groups,species,nsamples=10000){
  require(rjags)
  
  ##Removes negative and NA values from x, y and group vectors.
  data <- data.frame(x,y,groups,sdd)
  data_complete <- data[complete.cases(data),]
  
  ##Prepares input.
  group <- as.numeric(as.factor(data_complete$groups))
  ngroups <- length(unique(group))
  x <- data_complete$x
  y <- as.numeric(data_complete$y > 0)
  n <- length(x)
  sdd <- as.vector(tapply(data_complete$sdd, data_complete$groups, function(x) x[1]))
  
  ##Fits model.
  cat(
    "
    model{
    # priors
    height ~ dnorm(0,0.001)
    width.mu ~ dnorm(-10,0.001)T(,0)
    width.sigma ~ dunif(5,20)
    width.tau <- pow(width.sigma,-2)
    opt.mu ~ dnorm(0,0.001)T(-3,3)
    opt.sigma ~ dunif(0.01,5)
    opt.tau <- pow(opt.sigma,-2)
    opt.b1 ~ dnorm(1,10)
    
    
    for (j in 1:ngroups){
    opt.hat[j] <- opt.mu + opt.b1 * sdd[j]
    opt.g[j] ~ dnorm(opt.hat[j], opt.tau)
    width.g[j] ~ dnorm(width.mu,width.tau)T(,0)
    }
    
    # likelihood
    for (i in 1:n){
    y[i] ~ dbinom(p[i],1)
    p[i] <- 1 / (1 + exp(-z[i]))
    z[i] <- width.g[group[i]] * (x[i] - opt.g[group[i]])^2 + height
    }
    }  
    ", file="jagsmodel_log_vertexform.txt"
  )
  #   inits <-expression(list(height=rnorm(1,0,1),
  #                           width=rnorm(1,-10,1),
  #                           opt.mu=rnorm(1,0,1),
  #                           opt.sigma=runif(1,0.01,20),
  #                           b.opt=rnorm(1,0,1)))
  #   init.list <- list(eval(inits),eval(inits),eval(inits))
  jd <- list(x=x, y=y, n=n,group=group,ngroups=ngroups,sdd=sdd)
  mod <- jags.model("jagsmodel_log_vertexform.txt", data= jd, n.chains=3, n.adapt=2000)
  update(mod,n.iter=nsamples)
  out <- coda.samples(mod, c("height","opt.mu","opt.sigma","width.mu","width.sigma","opt.b1"),
                      n.iter=nsamples,thin=10)
  #diags <- gelman.diag(out)
  return(list(mod=mod,out=out))
  }

####Function to fit a Bayesian mixed-effects model in JAGS.####
fit.jags.mixed.allspp <- function(x,y,sdd,groups,species,years,nsamples=10000){
  require(rjags)
  
  ##Removes negative and NA values from x, y and group vectors.
  data <- data.frame(x,y,groups,species,years,sdd)
  data_complete <- data[complete.cases(data),]
  
  ##Prepares input.
  group <- as.numeric(as.factor(data_complete$groups))
  ngroups <- length(unique(group))
  species <- as.numeric(as.factor(data_complete$species))
  nspp <- length(unique(species))
  years <- as.numeric(as.factor(data_complete$years))
  nyears <- length(unique(years))
  x <- data_complete$x
  y <- as.numeric(data_complete$y > 0)
  n <- length(x)
  sdd <- as.vector(tapply(data_complete$sdd, data_complete$groups, function(x) x[1]))
  
  ##Fits model.
  cat(
    "
    model{
    # priors
    height.g.sigma ~ dunif(0,3.5)
    height.g.tau <- pow(height.g.sigma,-2)
    width.g.sigma ~ dunif(0.01,15)
    width.g.tau <- pow(width.g.sigma,-2)
    opt.g.sigma ~ dunif(0.01,5)
    opt.g.tau <- pow(opt.g.sigma,-2)
    
    height.s.mu ~ dnorm(0,0.001)
    height.s.sigma ~ dunif(0,100)
    height.s.tau <- pow(height.s.sigma,-2)
    opt.s.mu ~ dnorm(-0.5,0.001)
    opt.s.sigma ~ dunif(0,100)
    opt.s.tau <- pow(opt.s.sigma,-2)
    width.s.mu ~ dnorm(-10,0.001)T(,-1)
    width.s.sigma ~ dunif(0,80)
    width.s.tau <- pow(width.s.sigma,-2)
    
    for (j in 1:ngroups){
    height.g[j] ~ dnorm(0,height.g.tau)
    opt.g[j] ~ dnorm(0, opt.g.tau)
    width.g[j] ~ dnorm(0,width.g.tau)
    }
    
    for (l in 1:nspp){
    height.s[l] ~ dnorm(height.s.mu,height.s.tau)
    opt.s[l] ~ dnorm(opt.s.mu, opt.s.tau)
    width.s[l] ~ dnorm(width.s.mu,width.s.tau)T(,-1)
    }
    
    
    # likelihood
    for (i in 1:n){
    y[i] ~ dbinom(p[i],1)
    logit(p[i]) <- (width.s[spp[i]] + width.g[group[i]]) * 
                    (x[i] - ( opt.s[spp[i]] + opt.g[group[i]]))^2 + 
                    height.s[spp[i]] + height.g[group[i]]
    }
    }  
    ", file="jagsmodel_log_vertexform_allspp.txt"
  )
  #   inits <-expression(list(height=rnorm(1,0,1),
  #                           width=rnorm(1,-10,1),
  #                           opt.mu=rnorm(1,0,1),
  #                           opt.sigma=runif(1,0.01,20),
  #                           b.opt=rnorm(1,0,1)))
  #   init.list <- list(eval(inits),eval(inits),eval(inits))
  jd <- list(x=x, y=y, n=n,group=group,spp=species,ngroups=ngroups,nspp=nspp)
  mod <- jags.model("jagsmodel_log_vertexform_allspp.txt", data= jd, n.chains=3, n.adapt=1000)
  update(mod,n.iter=nsamples)
  out <- coda.samples(mod, c("height.s","opt.s","width.s","width.s.mu","width.s.sigma","height.s.mu",
                             "height.s.sigma","opt.s.mu","opt.s.sigma"),
                      n.iter=nsamples,thin=10)
  #diags <- gelman.diag(out)
  return(list(mod=mod,out=out))
  }

####Function to fit a Bayesian mixed-effects model in JAGS.####
fit.jags.mixed.allspp.mw <- function(x,y,sdd,groups,observers,species,years,nsamples=10000){
  require(rjags)
  
  ##Removes negative and NA values from x, y and group vectors.
  data <- data.frame(x,y,groups,observers,species,years,sdd)
  data_complete <- data[complete.cases(data),]
  
  ##Prepares input.
  group <- as.numeric(as.factor(data_complete$groups))
  ngroups <- length(unique(group))
  observers <- as.numeric(as.factor(data_complete$observers))
  nobs <- length(unique(observers))
  species <- as.numeric(as.factor(data_complete$species))
  nspp <- length(unique(species))
  years <- as.numeric(as.factor(data_complete$years))
  nyears <- length(unique(years))
  x <- data_complete$x
  y <- as.numeric(data_complete$y > 0)
  n <- length(x)
  sdd <- as.vector(tapply(data_complete$sdd, data_complete$groups, function(x) x[1]))
  
  ##Fits model.
  cat(
    "
    model{
    # priors
    height.g.sigma ~ dunif(0,3)
    height.g.tau <- pow(height.g.sigma,-2)
    opt.g.sigma ~ dunif(0.01,20)
    opt.g.tau <- pow(opt.g.sigma,-2)
    width.g.sigma ~ dunif(0.01,10)
    width.g.tau <- pow(width.g.sigma,-2)
    
    height.s.mu ~ dnorm(0,0.01)
    height.s.sigma ~ dunif(0,100)
    height.s.tau <- pow(height.s.sigma,-2)
    opt.s.mu ~ dnorm(-0.5,0.01)
    opt.s.sigma ~ dunif(0,1)
    opt.s.tau <- pow(opt.s.sigma,-2)
    width.s.mu ~ dnorm(-10,0.1)T(-100,-5)
    width.s.sigma ~ dunif(0,20)
    width.s.tau <- pow(width.s.sigma,-2)

    height.o.sigma ~ dunif(0,20)
    height.o.tau <- pow(height.o.sigma,-2)
    
    for (j in 1:ngroups){
    height.g[j] ~ dnorm(0,height.g.tau)
    opt.g[j] ~ dnorm(0, opt.g.tau)
    width.g[j] ~ dnorm(0,width.g.tau)
    }

    for (k in 1:nobs){
    height.o[k] ~ dnorm(0,height.o.tau)
    }
    
    for (l in 1:nspp){
    height.s[l] ~ dnorm(height.s.mu,height.s.tau)
    opt.s[l] ~ dnorm(opt.s.mu, opt.s.tau)
    width.s[l] ~ dnorm(width.s.mu,width.s.tau)T(-100,-5)
    }
    
    
    # likelihood
    for (i in 1:n){
    y[i] ~ dbinom(p[i],1)
    logit(p[i]) <- (width.s[spp[i]] + width.g[group[i]]) * 
    (x[i] - ( opt.s[spp[i]] + opt.g[group[i]]))^2 + 
    height.s[spp[i]] + height.g[group[i]] + height.o[observer[i]]
    }
    }  
    ", file="jagsmodel_log_vertexform_allspp.txt"
  )
  #   inits <-expression(list(height=rnorm(1,0,1),
  #                           width=rnorm(1,-10,1),
  #                           opt.mu=rnorm(1,0,1),
  #                           opt.sigma=runif(1,0.01,20),
  #                           b.opt=rnorm(1,0,1)))
  #   init.list <- list(eval(inits),eval(inits),eval(inits))
  jd <- list(x=x, y=y, n=n,group=group,spp=species,observer=observers,ngroups=ngroups,nspp=nspp,nobs=nobs)
  mod <- jags.model("jagsmodel_log_vertexform_allspp.txt", data= jd, n.chains=3, n.adapt=1000)
  update(mod,n.iter=nsamples)
  out <- coda.samples(mod, c("height.s","opt.s","width.s","width.s.mu","width.s.sigma","height.s.mu",
                             "height.s.sigma","opt.s.mu","opt.s.sigma"),
                      n.iter=nsamples,thin=10)
  #diags <- gelman.diag(out)
  return(list(mod=mod,out=out))
  }

####Function to update a JAGS model.
update.jags.mixed <- function(jagsmodel,n.update,n.iter,thin,
                              params=c("height.s","opt.s","width.s","width.s.mu","width.s.sigma","height.s.mu",
                                                   "height.s.sigma","opt.s.mu","opt.s.sigma")){
  mod <- jagsmodel$mod
  update(mod,n.iter=n.update)
  out <- coda.samples(mod, params,
                      n.iter=n.iter,thin=thin)
  return(list(mod=mod,out=out))
}
                              
####Function to fit a Bayesian mixed-effects model in JAGS.####
fit.jags.mixed.mw <- function(x,y,sdd,groups,species,nsamples=10000){
  require(rjags)
  
  ##Removes negative and NA values from x, y and group vectors.
  data <- data.frame(x,y,groups,sdd)
  data_complete <- data[complete.cases(data),]
  
  ##Prepares input.
  group <- as.numeric(as.factor(data_complete$groups))
  ngroups <- length(unique(group))
  x <- data_complete$x
  y <- as.numeric(data_complete$y > 0)
  n <- length(x)
  sdd <- as.vector(tapply(data_complete$sdd, data_complete$groups, function(x) x[1]))
  
  ##Fits model.
  cat(
    "
    model{
    # priors
    height ~ dnorm(0,0.001)
    width.mu ~ dnorm(-10,0.01)T(,-1)
    width.sigma ~ dunif(0.01,100)
    width.tau <- pow(width.sigma,-2)
    opt.mu ~ dnorm(0,0.001)T(-3,3)
    opt.sigma ~ dunif(0.01,20)
    opt.tau <- pow(opt.sigma,-2)
    
    for (j in 1:ngroups){
    opt.g[j] ~ dnorm(opt.mu, opt.tau)
    width.g[j] ~ dnorm(width.mu,width.tau)T(,-1)
    }
    
    # likelihood
    for (i in 1:n){
    y[i] ~ dbinom(p[i],1)
    p[i] <- 1 / (1 + exp(-z[i]))
    z[i] <- width.g[group[i]] * (x[i] - opt.g[group[i]])^2 + height
    }
    }  
    ", file="jagsmodel_log_vertexform.txt"
  )
  #   inits <-expression(list(height=rnorm(1,0,1),
  #                           width=rnorm(1,-10,1),
  #                           opt.mu=rnorm(1,0,1),
  #                           opt.sigma=runif(1,0.01,20),
  #                           b.opt=rnorm(1,0,1)))
  #   init.list <- list(eval(inits),eval(inits),eval(inits))
  jd <- list(x=x, y=y, n=n,group=group,ngroups=ngroups,sdd=sdd)
  mod <- jags.model("jagsmodel_log_vertexform.txt", data= jd, n.chains=3, n.adapt=2000)
  update(mod,n.iter=nsamples)
  out <- coda.samples(mod, c("height","opt.mu","opt.sigma","width.mu","width.sigma"),
                      n.iter=nsamples,thin=10)
  #diags <- gelman.diag(out)
  return(list(mod=mod,out=out))
  }

####Function to fit a Bayesian mixed-effects model in JAGS.####
fit.jags.mixed.fl <- function(x,y,groups,species,nsamples=10000,thin=10){
  require(rjags)
  
  ##Removes negative and NA values from x, y and group vectors.
  data <- data.frame(x,y,groups)
  data_complete <- data[complete.cases(data),]
  
  ##Prepares input.
  group <- as.numeric(as.factor(data_complete$groups))
  ngroups <- length(unique(group))
  x <- data_complete$x
  y <- as.numeric(data_complete$y > 0)
  n <- length(x)
  
  ##Fits model.
  cat(
    "
    model{
    # priors
    height ~ dnorm(0,0.001)
    width.mu ~ dnorm(-10,0.01)T(,-0.1)
    opt.mu ~ dnorm(0,0.001)T(-3,3)
    opt.sigma ~ dunif(0.01,5)
    opt.tau <- pow(opt.sigma,-2)
    
    for (j in 1:ngroups){
    opt.g[j] ~ dnorm(opt.mu, opt.tau)
    }
    
    # likelihood
    for (i in 1:n){
    y[i] ~ dbinom(p[i],1)
    p[i] <- 1 / (1 + exp(-z[i]))
    z[i] <- width.mu * (x[i] - opt.g[group[i]])^2 + height
    }
    }  
    ", file="jagsmodel_vertexform.txt"
  )
  #   inits <-expression(list(height=rnorm(1,0,1),
  #                           width=rnorm(1,-10,1),
  #                           opt.mu=rnorm(1,0,1),
  #                           opt.sigma=runif(1,0.01,20),
  #                           b.opt=rnorm(1,0,1)))
  #   init.list <- list(eval(inits),eval(inits),eval(inits))
  jd <- list(x=x, y=y, n=n,group=group,ngroups=ngroups)
  mod <- jags.model("jagsmodel_vertexform.txt", data= jd, n.chains=3, n.adapt=2000)
  update(mod,n.iter=nsamples)
  out <- coda.samples(mod, c("height","opt.mu","opt.sigma","width.mu"),
                      n.iter=nsamples,thin=thin)
  #diags <- gelman.diag(out)
  return(list(mod=mod,out=out))
  }

####Function to fit a Bayesian mixed-effects model in JAGS.####
fit.jags.mixed.allspp.fl <- function(x,y,groups,observers,species,years,nsamples=10000){
  require(rjags)
  
  ##Removes negative and NA values from x, y and group vectors.
  data <- data.frame(x,y,groups,observers,species,years)
  data_complete <- data[complete.cases(data),]
  
  ##Prepares input.
  group <- as.numeric(as.factor(data_complete$groups))
  ngroups <- length(unique(group))
  observers <- as.numeric(as.factor(data_complete$observers))
  nobs <- length(unique(observers))
  species <- as.numeric(as.factor(data_complete$species))
  nspp <- length(unique(species))
  years <- as.numeric(as.factor(data_complete$years))
  nyears <- length(unique(years))
  x <- data_complete$x
  y <- as.numeric(data_complete$y > 0)
  n <- length(x)

  ##Fits model.
  cat(
    "
    model{
    # priors
    opt.g.sigma ~ dgamma(5,10)
    opt.g.tau <- pow(opt.g.sigma,-2)
    
    height.s.mu ~ dnorm(-1,0.01)
    height.s.sigma ~ dunif(0.001,100)
    height.s.tau <- pow(height.s.sigma,-2)
    opt.s.mu ~ dnorm(-0.6,10)
    opt.s.sigma ~ dunif(0.001,100)
    opt.s.tau <- pow(opt.s.sigma,-2)
    width.s.mu ~ dnorm(-20,10)T(-120,-0.5)
    width.s.sigma ~ dgamma(5,1)
    width.s.tau <- pow(width.s.sigma,-2)
    
    height.o.sigma ~ dgamma(1,0.1)
    height.o.tau <- pow(height.o.sigma,-2)
    
    for (j in 1:ngroups){
    opt.g[j] ~ dnorm(0, opt.g.tau)
    }
    
    for (k in 1:nobs){
    height.o[k] ~ dnorm(0,height.o.tau)
    }
    
    for (l in 1:nspp){
    height.s[l] ~ dnorm(height.s.mu,height.s.tau)
    opt.s[l] ~ dnorm(opt.s.mu, opt.s.tau)
    width.s[l] ~ dnorm(width.s.mu,width.s.tau)T(-120,-0.5)
    }
    
    
    # likelihood
    for (i in 1:n){
    y[i] ~ dbinom(p[i],1)
    logit(p[i]) <- (width.s[spp[i]]) * 
    (x[i] - ( opt.s[spp[i]] + opt.g[group[i]]))^2 + 
    height.s[spp[i]] + height.o[observer[i]]
    }
    }  
    ", file="jagsmodel_log_vertexform_allspp.txt"
  )
  #   inits <-expression(list(height=rnorm(1,0,1),
  #                           width=rnorm(1,-10,1),
  #                           opt.mu=rnorm(1,0,1),
  #                           opt.sigma=runif(1,0.01,20),
  #                           b.opt=rnorm(1,0,1)))
  #   init.list <- list(eval(inits),eval(inits),eval(inits))
  jd <- list(x=x, y=y, n=n,group=group,spp=species,observer=observers,ngroups=ngroups,nspp=nspp,nobs=nobs)
  mod <- jags.model("jagsmodel_log_vertexform_allspp.txt", data= jd, n.chains=3, n.adapt=1000)
  update(mod,n.iter=nsamples)
  out <- coda.samples(mod, c("height.s","opt.s","width.s","width.s.mu",
                                    "width.s.sigma","height.s.mu",
                                    "height.s.sigma","opt.s.mu","opt.s.sigma",
                                    "opt.g.sigma","height.o.sigma"),
                      n.iter=nsamples,thin=5)
  #diags <- gelman.diag(out)
  return(list(mod=mod,out=out))
  }


####Gets credible intervals for the parameters of interest.####
get.param.creds <- function(out,params=c("width.mu","opt.mu","height","opt.b1"),intervals=c(0.025,0.5,0.975)){
  gg_out <- ggmcmc::ggs(out)
  samples <- matrix(NA,ncol=length(params),nrow=(length(out[[1]])*length(out)))
  for (i in 1:length(params)){
    samples[,i] <- subset(gg_out,Parameter==params[i])$value
  }
  cred.fun <- function(x){quantile(x,probs=intervals)}
  creds <- apply(samples,FUN=cred.fun,MARGIN=2)
  return(creds)
}

####Gets 95% credible intervals on the fit mean curve.####
get.curve.creds <- function(out,xnew,
                            params=c("width.mu","opt.mu","height","opt.b1"),
                            probs=c(0.025,0.5,0.975)){
  gg_out <- ggmcmc::ggs(out)
  samples <- matrix(NA,ncol=length(params),nrow=(length(out[[1]])*length(out)))
  for (i in 1:length(params)){
    samples[,i] <- subset(gg_out,Parameter==params[i])$value
  }
  pred.fun <- function(x){antilogit(x[1] * (xnew - x[2])^2 + x[3])}
  preds <- apply(samples,FUN=pred.fun,MARGIN=1)
  preds_t <- t(preds)
  cred.fun <- function(x){quantile(x,probs=probs,na.rm=TRUE)}
  creds <- apply(preds_t,FUN=cred.fun,MARGIN=2)
  creds_t <- t(creds)
  return(creds_t)
}

####Gets 95% credible intervals on the fit mean curve.####
get.curve.intervals <- function(out,xnew,x_unscaled,
                                params=c("width.mu","opt.mu","height"),
                                probs=c(0.025,0.5,0.975)){
  gg_out <- ggmcmc::ggs(out)
  samples <- matrix(NA,ncol=length(params),nrow=(length(out[[1]])*length(out)))
  for (i in 1:length(params)){
    samples[,i] <- subset(gg_out,Parameter==params[i])$value
  }
  pred.fun <- function(x){antilogit(x[1] * (xnew - x[2])^2 + x[3])}
  preds <- apply(samples,FUN=pred.fun,MARGIN=1)
  preds_t <- t(preds)
  interval.fun <- function(x){
    preds_df <- data.frame(dss=x_unscaled,pred=x)
    intervals <- obs_intervals(preds_df)
    width <- intervals[2] - intervals[1]
    names(width) <- "width"
    max <- preds_df$dss[which.max(preds_df$pred)]
    return(c(intervals,opt=max,width))
  }
  creds <- apply(preds_t,FUN=interval.fun,MARGIN=1)
  creds_t <- t(creds)
  cred.fun <- function(x){quantile(x,probs=probs,na.rm=TRUE)}
  interval_creds <- apply(creds_t,FUN=cred.fun,MARGIN=2)
  return(interval_creds)
}

##Function to parse JAGS output into a matrix of functions representing posterior samples of response curves
##for each species.
make.fun.matrix <- function(jags.out,par.names=c("width.s","opt.s","height.s"),
                            n.samples=100){
  gg.out <- ggmcmc::ggs(jags.out)
  ngroups <- sum(grepl(paste(par.names[1],"[",sep=""),unique(gg.out$Parameter),fixed=TRUE))
  samples <- sample(1:dim(jags.out[[1]])[1],size=n.samples,replace=FALSE)
  chains <- sample(1:3,size=n.samples,replace=TRUE)
  gg.sample <- gg.out[gg.out$Iteration %in% samples & gg.out$Chain == chains,]
  fun.list <- replicate(ngroups*n.samples,function(x){x})
  fun.matrix <- matrix(fun.list,ncol=ngroups)
  print("Making function matrix for species ")
  for(j in 1:ngroups){
    print(paste(j,"of ",ngroups))
    pname1 <- paste(par.names[1],"[",j,"]",sep="")
    pname2 <- paste(par.names[2],"[",j,"]",sep="")
    pname3 <- paste(par.names[3],"[",j,"]",sep="")
    p1.gr <- gg.sample[gg.sample$Parameter == pname1,4]
    p2.gr <- gg.sample[gg.sample$Parameter == pname2,4]
    p3.gr <- gg.sample[gg.sample$Parameter == pname3,4]
    np <- length(p3.gr)
    fun.body <- paste(rep("antilogit(",np),p1.gr,rep(" * (x - ",np),p2.gr,rep(")^2 + ",np),
                      p3.gr,rep(")",np),sep="")
    for (i in 1:n.samples){
      body(fun.matrix[i,j][[1]]) <- parse(text=fun.body[i])
    }
  }
  return(fun.matrix)
}

##Function to parse JAGS output into a matrix of functions representing posterior samples of response curves
##for each species.
make.med.fun.matrix <- function(jags.out,par.names=c("width.s","opt.s","height.s"),
                            n.samples=100,spp=TRUE){
  gg.out <- ggmcmc::ggs(jags.out)
  if(spp==FALSE){
    ngroups <- 1
  }else{
    ngroups <- sum(grepl(paste(par.names[1],"[",sep=""),unique(gg.out$Parameter),fixed=TRUE))
  }
  fun.list <- replicate(ngroups*n.samples,function(x){x})
  fun.matrix <- matrix(fun.list,ncol=ngroups)
  print("Making function matrix for species ")
  for(j in 1:ngroups){
    print(paste(j,"of ",ngroups))
    if(spp==FALSE){
      pname1 <- par.names[1]
      pname2 <- par.names[2]
      pname3 <- par.names[3]
    }else{
      pname1 <- paste(par.names[1],"[",j,"]",sep="")
      pname2 <- paste(par.names[2],"[",j,"]",sep="")
      pname3 <- paste(par.names[3],"[",j,"]",sep="")
    }
    p1.gr <- quantile(gg.out[gg.out$Parameter == pname1,4],probs=c(0.5))
    p2.gr <- quantile(gg.out[gg.out$Parameter == pname2,4],probs=c(0.5))
    p3.gr <- quantile(gg.out[gg.out$Parameter == pname3,4],probs=c(0.5))
    np <- length(p3.gr)
    fun.body <- paste(rep("antilogit(",np),p1.gr,rep(" * (x - ",np),p2.gr,rep(")^2 + ",np),
                      p3.gr,rep(")",np),sep="")
    body(fun.matrix[1,j][[1]]) <- parse(text=fun.body)
    }
  return(fun.matrix)
}

make.community.fun.matrix <- function(jags.out,par.names=c("width.s.mu","opt.s.mu","height.s.mu"),
                                      n.samples=1000){
    gg.out <- ggmcmc::ggs(jags.out)
    samples <- sample(1:dim(jags.out[[1]])[1],size=n.samples,replace=FALSE)
    chains <- sample(1:3,size=n.samples,replace=TRUE)
    gg.sample <- gg.out[gg.out$Iteration %in% samples & gg.out$Chain == chains,]
    fun.list <- replicate(n.samples,function(x){x})
    fun.matrix <- matrix(fun.list,ncol=1)
    print("Making function matrix for iteration.")
    p1 <- gg.sample[gg.sample$Parameter == par.names[1],4]
    p2 <- gg.sample[gg.sample$Parameter == par.names[2],4]
    p3 <- gg.sample[gg.sample$Parameter == par.names[3],4]
    np <- length(p3)
    fun.body <- paste(rep("antilogit(",np),p1,rep(" * (x - ",np),p2,rep(")^2 + ",np),
                  p3,rep(")",np),sep="")
    for (i in 1:n.samples){
      print(paste(i,"of ",n.samples))
      body(fun.matrix[i,1][[1]]) <- parse(text=fun.body[i])
    }
    return(fun.matrix)
}





