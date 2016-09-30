##Functions for Mt. Rainier microclimate analysis.

##Logit and inverse logit function
logit <- function(x) {log(x/(1-x))}
inv.logit <- function(x) {exp(x)/(exp(x)+1)}

## Function for fitting an OLS model for each day.
model.temps.lm <- function(dataset,daterange=750:2800){
  
  date_range <- daterange
  
  ## Data frame to hold linear model coefficients and residual spatial parameters.
  results <- data.frame(DATE=dataset$DATE[date_range],
                        intercepts = rep(NA,length(date_range)),
                        lapses = rep(NA,length(date_range)),
                        canopies = rep(NA,length(date_range)),
                        relevs = rep(NA,length(date_range)),
                        rsqs = rep(NA,length(date_range)),
                        aics = rep(NA,length(date_range)),
                        maes = rep(NA,length(date_range)),
                        nuggets = rep(NA,length(date_range)),
                        ranges = rep(NA,length(date_range)))
  
  for (i in date_range){
    index <- (i - min(date_range))+1
    day <- t(dataset[i,2:dim(dataset)[2]])
    colnames(day) <- "temp"
    day_meta <- merge(data.frame(alt_code=rownames(day),tmax=day),meta)
    #plot(temp~elev,data=day_meta,main=dataset$DATE[date_range[index]])
    #text(x=day_meta$elev,y=day_meta$temp,labels=day_meta$alt_code)
    test_lm <- lm(temp~elev+I(asin(ccov_81m))+relev_729m+utmx+utmy,data=day_meta)
    results$intercepts[index] <- test_lm$coefficients[1]
    results$lapses[index] <- test_lm$coefficients[2]
    results$canopies[index] <- test_lm$coefficients[3]
    results$relevs[index] <- test_lm$coefficients[4]
    results$rsqs[index] <- summary(test_lm)$r.squared
    results$aics[index] <- AIC(test_lm)
    summary(test_lm)
    temp_resid <- data.frame(day_meta[complete.cases(day_meta),],residuals(test_lm))
    results$maes[index] <- mean(abs(temp_resid$residuals.test_lm.))
    coordinates(temp_resid) <- ~utmx+utmy
    
    ##Residual variogram.
    resid_variog <- variog(coords=coordinates(temp_resid),data=temp_resid$residuals.test_lm.,
                           max.dist=20000,option="bin",breaks=seq(0,20000,by=1000),messages=FALSE)
    resid_fit <- variofit(resid_variog,cov.model="spherical",fix.nugget=FALSE,max.dist=20000,messages=FALSE)
    results$nuggets[index] <- resid_fit$cov.pars[1]
    results$ranges[index] <- resid_fit$cov.pars[2]
    #plot(resid_variog)
    #lines(resid_fit)
  }
  return(results)
}

##Function for diagnostic plotting.
plot.diags <- function(dataset,maintext=""){
  
  tickpos <- seq(as.POSIXct("2008-01-01", tz="GMT"),
                 as.POSIXct("2014-01-01", tz="GMT"),
                 by="year")
  par(mfrow=c(8,1),mar=c(0,4,1,1),oma=c(4,0,2,0))
  plot(dataset$DATE,type="l",xaxt="n",dataset$lapses*1000,ylab="Lapse Rate",ylim=c(-9,9))
  lapse_sm <- loess((dataset$lapses*1000)~as.numeric(dataset$DATE),span=0.2)
  lines(dataset$DATE,predict(lapse_sm),col="red",lwd=2)
  abline(h=0,lty=2)
  abline(v=tickpos,col="grey60")
  axis.POSIXct(side = 1,at=tickpos,labels=FALSE)
  text(x=as.numeric(max(dataset$DATE)),y=7,labels=round(mean(dataset$lapses*1000),digits=3))
  
  plot(dataset$DATE,type="l",xaxt="n",dataset$canopies,ylab="Can. Cov.",ylim=c(-8,8))
  can_sm <- loess(dataset$canopies~as.numeric(dataset$DATE),span=0.2)
  lines(dataset$DATE,predict(can_sm),col="red",lwd=2)
  abline(h=0,lty=2)
  abline(v=tickpos,col="grey60")
  axis.POSIXct(side = 1,at=tickpos,labels=FALSE)
  text(x=as.numeric(max(dataset$DATE)),y=3,labels=round(mean(dataset$canopies),digits=3))
  
  plot(dataset$DATE,type="l",xaxt="n",dataset$relevs,ylab="R. elev.",ylim=c(-0.1,0.1))
  relev_sm <- loess(dataset$relevs~as.numeric(dataset$DATE),span=0.2)
  lines(dataset$DATE,predict(relev_sm),col="red",lwd=2)
  abline(h=0,lty=2)
  abline(v=tickpos,col="grey60")
  axis.POSIXct(side = 1,at=tickpos,labels=FALSE)
  text(x=as.numeric(max(dataset$DATE)),y=0.08,labels=round(mean(dataset$relevs),digits=3))
  
  plot(dataset$DATE,type="l",xaxt="n",dataset$rsqs,ylim=c(0,1),ylab=expression(r^2))
  text(x=as.numeric(max(dataset$DATE)),y=0.1,labels=round(mean(dataset$rsqs),digits=3))
  abline(v=tickpos,col="grey60")
  axis.POSIXct(side = 1,at=tickpos,labels=FALSE)
  
  plot(dataset$DATE,type="l",xaxt="n",dataset$aics,ylim=c(0,400),ylab="AIC")
  text(x=as.numeric(max(dataset$DATE)),y=50,labels=round(mean(dataset$aics),digits=3))
  abline(v=tickpos,col="grey60")
  axis.POSIXct(side = 1,at=tickpos,labels=FALSE)
  
  plot(dataset$DATE,type="l",xaxt="n",dataset$maes,ylab="Mean abs. err. (C)",ylim=c(0,2))
  text(x=as.numeric(max(dataset$DATE)),y=0.25,labels=round(mean(dataset$maes),digits=3))
  abline(v=tickpos,col="grey60")
  axis.POSIXct(side = 1,at=tickpos,labels=FALSE)
  
  plot(dataset$DATE,type="l",xaxt="n",dataset$nuggets,ylab="Resid. Nugget")
  text(x=as.numeric(max(dataset$DATE)),y=2,labels=round(mean(dataset$nuggets),digits=3))
  abline(v=tickpos,col="grey60")
  axis.POSIXct(side = 1,at=tickpos,labels=FALSE)
  
  plot(dataset$DATE,type="l",y=(dataset$ranges/1000),ylab="Resid. Range (km)",ylim=c(0,100))
  abline(v=tickpos,col="grey60")
  mtext(text=maintext,side=3,outer=TRUE)
}

##Custom color ramp function from (http://stackoverflow.com/questions/13327326/r-image-function-in-r)
color.palette <- function(steps, n.steps.between=NULL, ...){
  
  if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
  if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")
  
  fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
  RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
  RGB[,fill.steps] <- col2rgb(steps)
  
  for(i in which(n.steps.between>0)){
    col.start=RGB[,fill.steps[i]]
    col.end=RGB[,fill.steps[i+1]]
    for(j in seq(3)){
      vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]  
      RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
    }
  }
  
  new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
  pal <- colorRampPalette(new.steps, ...)
  return(pal)
}


##Bootstrap functions.

ci_boot <- function(data,rows){
  boot_lm <- lm(data[rows,4]~data[rows,3])
  disp <- boot_lm$coefficients[1]
  sens <- boot_lm$coefficients[2]
  decoup <- summary(boot_lm)$sigma
  return(c(disp,sens,decoup))
}

boot_indices <- function(regional_series,local_frame,site_names,nboot=100){
  require(boot)
  par(mfrow=c(3,3))
  nsites <- length(site_names)
  indices<- data.frame(site=site_names,
                       disp=rep(NA,nsites),
                       disp_lwr95=rep(NA,nsites),
                       disp_upr95=rep(NA,nsites),
                       sens=rep(NA,nsites),
                       sens_lwr95=rep(NA,nsites),
                       sens_upr95=rep(NA,nsites),
                       decouple=rep(NA,nsites),
                       decouple_lwr95=rep(NA,nsites),
                       decouple_upr95=rep(NA,nsites))
  
  ##Examines only columns with enough data.
  compare_cols <- colSums(local_frame[1:length(regional_series),],na.rm=TRUE) != 0.00 &
    colSums(!is.na(local_frame[1:length(regional_series),]),na.rm=TRUE) > 60
  
  for (i in which(compare_cols)){
    title <- site_names[i]
    local_meas <- local_frame[1:length(regional_series),i]
    if(compare_cols[i]){
      
      ##Puts complete data into a data frame.
      comp_dat <- data.frame(regional_series,local_meas)
      colnames(comp_dat) <- c("regional_series","local_meas")
      comp_dat <- comp_dat[complete.cases(comp_dat),]
      
      ##Centers regional and local temp.
      comp_dat$regional_series_center <- scale(comp_dat$regional_series,scale=FALSE,
                                               center=TRUE)
      comp_dat$local_meas_center <- comp_dat$local_meas - attr(comp_dat$regional_series_center,"scaled:center")
      
      ##Plots the relationship
      plot(x=comp_dat$regional_series_center,y=comp_dat$local_meas_center,
           xlab="centered regional temp.",ylab="centered local temp.",main=title,pch=".",
           ylim=c(-20,30),xlim=c(-20,30),col="blue")
      abline(0,1,lty=1)
      
      ##Number of bootstrap samples.
      boot_samples <- nboot
      
      ##Computes the intercept,slope and variance.
      lm_est_pr <- lm(local_meas_center~regional_series_center,data=comp_dat)
      disp_obs_pr <-lm_est_pr$coefficients[1]
      sens_obs_pr <- lm_est_pr$coefficients[2]
      decoup_obs_pr <- summary(lm_est_pr)$sigma
      
      ##Computes boostrap confidence intervals for all parameters.
      allparam_boot <- boot(data=comp_dat,statistic=ci_boot,R=boot_samples,parallel = "multicore")
      disp_boot_ci <- boot.ci(allparam_boot,conf=0.95,type="basic",index=1)
      sens_boot_ci <- boot.ci(allparam_boot,conf=0.95,type="basic",index=2)
      decoup_boot_ci <- boot.ci(allparam_boot,conf=0.95,type="basic",index=3)
      
      ##Adds to output
      indices[i,2] <- disp_obs_pr
      indices[i,3:4] <- disp_boot_ci$basic[4:5]
      indices[i,5] <- sens_obs_pr
      indices[i,6:7] <- sens_boot_ci$basic[4:5]
      indices[i,8] <- decoup_obs_pr
      indices[i,9:10] <- decoup_boot_ci$basic[4:5]
      
      abline(lm_est_pr$coefficients[1],lm_est_pr$coefficients[2],lty=2,col="blue")
      text(-5,25,paste("Intercept:",round(disp_obs_pr,digits=3)),col="blue")
      text(-5,15,paste("Slope:",round(sens_obs_pr,digits=3)),col="blue")
      text(12,-10,paste("Sigma:",round(decoup_obs_pr,digits=3)),col="blue")
    }
  }
  return(indices)
}

## Function for plotting time-series output from the state-space model.
plot.beta.ts <-  function(y.dates,beta.0,beta.1,beta.2,beta.3,beta.4,beta.5){
  par(mfrow=c(6,1),mar=c(0,4,0,0),oma=c(4,1,1,1))
  plot(y.dates, beta.0[1,], pch=19, cex=0.5, xlab="", ylab="Intercept",
       xaxt="n",ylim=range(beta.0))
  arrows(y.dates, beta.0[1,], y.dates, beta.0[3,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  arrows(y.dates, beta.0[1,], y.dates, beta.0[2,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  abline(h=0,lty=2,col=2)
  xlabels <- seq(as.POSIXct("2014-9-1 00:00:00"),as.POSIXct("2015-10-1 00:00:00"),by='month')
  xlabels2 <- seq(as.POSIXct("2014-1-1 00:00:00"),as.POSIXct("2016-1-1 00:00:00"),by='year')
  axis.POSIXct(side=1, at=xlabels,format="%b",labels = FALSE)
  axis.POSIXct(side=1,at=xlabels2,format="%Y",tck=1,labels = FALSE,hadj=-3,padj=1.8)
  
  plot(y.dates, beta.1[1,], pch=19, cex=0.5, xlab="", ylab="Elevation",
       xaxt="n",ylim=c(-4,4))
  arrows(y.dates, beta.1[1,], y.dates, beta.1[3,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  arrows(y.dates, beta.1[1,], y.dates, beta.1[2,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  abline(h=0,lty=2,col=2)
  xlabels <- seq(as.POSIXct("2014-9-1 00:00:00"),as.POSIXct("2015-10-1 00:00:00"),by='month')
  xlabels2 <- seq(as.POSIXct("2014-1-1 00:00:00"),as.POSIXct("2016-1-1 00:00:00"),by='year')
  axis.POSIXct(side=1, at=xlabels,format="%b",labels = FALSE)
  axis.POSIXct(side=1,at=xlabels2,format="%Y",tck=1,labels = FALSE,hadj=-3,padj=1.8)
  
  plot(y.dates, beta.2[1,], pch=19, cex=0.5, xlab="", ylab="Cold-air",
       xaxt="n",ylim=c(-4,4))
  arrows(y.dates, beta.2[1,], y.dates, beta.2[3,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  arrows(y.dates, beta.2[1,], y.dates, beta.2[2,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  abline(h=0,lty=2,col=2)
  xlabels <- seq(as.POSIXct("2014-9-1 00:00:00"),as.POSIXct("2015-10-1 00:00:00"),by='month')
  xlabels2 <- seq(as.POSIXct("2014-1-1 00:00:00"),as.POSIXct("2016-1-1 00:00:00"),by='year')
  axis.POSIXct(side=1, at=xlabels,format="%b",labels = FALSE)
  axis.POSIXct(side=1,at=xlabels2,format="%Y",tck=1,labels = FALSE,hadj=-3,padj=1.8)
  
  plot(y.dates, beta.3[1,], pch=19, cex=0.5, xlab="", ylab="Can. Vol.",
       xaxt="n",ylim=c(-4,4))
  arrows(y.dates, beta.3[1,], y.dates, beta.3[3,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  arrows(y.dates, beta.3[1,], y.dates, beta.3[2,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  abline(h=0,lty=2,col=2)
  xlabels <- seq(as.POSIXct("2014-9-1 00:00:00"),as.POSIXct("2015-10-1 00:00:00"),by='month')
  xlabels2 <- seq(as.POSIXct("2014-1-1 00:00:00"),as.POSIXct("2016-1-1 00:00:00"),by='year')
  axis.POSIXct(side=1, at=xlabels,format="%b",labels = TRUE)
  axis.POSIXct(side=1,at=xlabels2,format="%Y",tck=1,labels = TRUE,hadj=-3,padj=1.8)
  
  plot(y.dates, beta.4[1,], pch=19, cex=0.5, xlab="", ylab="UTM X",
       xaxt="n",ylim=c(-4,4))
  arrows(y.dates, beta.4[1,], y.dates, beta.4[3,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  arrows(y.dates, beta.4[1,], y.dates, beta.4[2,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  abline(h=0,lty=2,col=2)
  xlabels <- seq(as.POSIXct("2014-9-1 00:00:00"),as.POSIXct("2015-10-1 00:00:00"),by='month')
  xlabels2 <- seq(as.POSIXct("2014-1-1 00:00:00"),as.POSIXct("2016-1-1 00:00:00"),by='year')
  axis.POSIXct(side=1, at=xlabels,format="%b",labels = TRUE)
  axis.POSIXct(side=1,at=xlabels2,format="%Y",tck=1,labels = TRUE,hadj=-3,padj=1.8)
  
  plot(y.dates, beta.5[1,], pch=19, cex=0.5, xlab="", ylab="UTM Y",
       xaxt="n",ylim=c(-4,4))
  arrows(y.dates, beta.4[1,], y.dates, beta.4[3,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  arrows(y.dates, beta.4[1,], y.dates, beta.4[2,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  abline(h=0,lty=2,col=2)
  xlabels <- seq(as.POSIXct("2014-9-1 00:00:00"),as.POSIXct("2015-10-1 00:00:00"),by='month')
  xlabels2 <- seq(as.POSIXct("2014-1-1 00:00:00"),as.POSIXct("2016-1-1 00:00:00"),by='year')
  axis.POSIXct(side=1, at=xlabels,format="%b",labels = TRUE)
  axis.POSIXct(side=1,at=xlabels2,format="%Y",tck=1,labels = TRUE,hadj=-3,padj=1.8)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#### Temporal Interpolation ####################################################
# Perform cell-wise linear interpolation between multiple raster layers, and 
# extrapolation beyond the upper limit of input data. Output is saved in .tif 
# format. From John Baums: https://gist.github.com/johnbaums/10465462
#
# Arguments
# s: a rasterStack containing the time slices to be interpolated 
#
# xin: a numeric vector that indicates the times associated with layers in s (in
# the same order as the layers of s - see names(s))
#
# xout: a numeric vector that indicates the times to be interpolated to (NB: if
# xout extends beyond the latest time slice in xin, it will be extrapolated
# using the rate of change from the period between the last and second to last
# time in xin.) 
#
# outdir: the directory to which files will be written (recursively created if 
# not already in existence) (character string)
#
# prefix: the output files will have pattern prefix_x.tif, where x is the
# timestep (potentially multiple digits), and prefix is a string that you 
# specify here (character string) 
#
# progress: show progress bar (TRUE/FALSE) 
#
# writechange: write the change grids that define the change in cell value per
# timestep between each pair of time slices, (TRUE/FALSE). If TRUE, these are
# written to outdir.
#
# returnstack: if TRUE, returns the interpolated grids (at timesteps xout) as a 
# rasterStack (TRUE/FALSE)
#
# ...: additional arguments passed to writeRaster

interpolateTemporal <- function(s, xin, xout, outdir, prefix, progress=TRUE, 
                                writechange=TRUE, returnstack=FALSE, ...) {
  require(raster)
  require(rgdal)
  if(nlayers(s) != length(xin)) stop('Length of xin must equal the number of layers in s.')
  if(nlayers(s) < 2) stop('stack s must have at least 2 layers.')
  if(!all(findInterval(xout, range(xin), rightmost.closed=TRUE) == 1)) {
    if(any(xout < min(xin))) {
      stop('This function does not extrapolate backwards (i.e. below the earliest element in xin). All elements of xout must be greater that min(xin).')
    } else {      
      warning('Some values of xout require extrapolation beyond the range of xin.\nThe rate of change for extrapolation is assumed to be equal to that for the period between the last and second-to-last elements of xin (after temporal ordering).')
    }
  }
  outdir <- normalizePath(sub('/$|\\\\$', '', outdir), winslash='/', 
                          mustWork=FALSE)
  if(!file.exists(outdir)) dir.create(outdir, recursive=TRUE)
  xout <- unique(xout)
  if(is.unsorted(xin)) {
    s <- s[[order(xin)]]
    xin <- sort(xin)
  }
  len <- diff(xin)
  base <- findInterval(xout, xin)
  lower <- unique(base[base < nlayers(s)])
  s.change <- stack(sapply(if(length(lower) > 0) lower else nlayers(s) - 1, 
                           function(x) {
                             message(sprintf('Calculating change grid for %s to %s.', xin[x], xin[x+1]))
                             overlay(s[[x]], s[[x+1]], fun=function(x1, x2) (x2-x1)/len[x],
                                     filename=ifelse(writechange, 
                                                     file.path(outdir, sprintf('changegrid_%s_%s', xin[x], xin[x+1])), 
                                                     ''), recycle=FALSE, format='GTiff', ...)
                           }))
  
  multi <- xout - xin[base]
  chg.ind <- ifelse(base > nlayers(s.change), nlayers(s.change), base)
  message('Calculating grids for years specified in xout...')
  if(progress) pb <- txtProgressBar(0, length(xout), style=3)
  invisible(sapply(seq_along(xout), function(x) {
    out.rast <- if(xout[x] %in% xin) {
      s[[base[x]]]
    } else {
      overlay(s[[base[x]]], s.change[[chg.ind[x]]],
              fun=function(x1, x2) x1 + (x2*multi[x]))
    }
    writeRaster(out.rast, 
                filename=file.path(outdir, sprintf('%s_%s', prefix, sprintf('%04d',xout[x]))), 
                format='GTiff', ...)
    if(progress) setTxtProgressBar(pb, x)
  }))
  if(isTRUE(returnstack)) {
    f <- file.path(outdir, paste0(prefix, '_', xout, '.tif'))
    return(stack(f[order(as.numeric(sub('.*_(\\d+)\\.tif$', '\\1', f)))]))
  }
}

##Function to fit a Bayesian spatial regression to daily data.
bayes_daily_anoms <- function(exp_brick,exp_brick_path,meas_df,meta,preds,out_prefix,nsamples,pp=TRUE,
                              overwrite=FALSE){
  require(raster)
  require(spBayes)
  
  setwd(exp_brick_path)
  temp_exp_year <- brick(exp_brick)
  
  stopifnot(nlayers(temp_exp_year) >= dim(meas_df)[1])
  stopifnot(all(colnames(meas_df[,-1]) %in% meta$alt_code))
  stopifnot("utmx_s" %in% names(preds))
  stopifnot("utmy_s" %in% names(preds))
  stopifnot("elev" %in% names(preds))
  stopifnot(all(!is.nan(rowMeans(meas_df[,-1],na.rm=TRUE))))
  
  temp_exp_year <- brick(exp_brick)
  start_num <- as.numeric(format(as.Date(meas_df$DATE[1]),format="%j"))
  end_num <- as.numeric(format(as.Date(meas_df$DATE[length(meas_df$DATE)]),format="%j"))
  
  temp_year <- meas_df[start_num:end_num,-1]
  temp_meta_year <- merge(data.frame(alt_code=colnames(temp_year)),
                          meta,by="alt_code",sort=FALSE)
  print(paste("Extracting expected values from raster",exp_brick))
  temp_year_exp <- t(extract(temp_exp_year,cbind(temp_meta_year$X_UTM,temp_meta_year$Y_UTM),
                             method="bilinear",sp=FALSE) / 10000)
  temp_year_exp <- temp_year_exp[start_num:end_num,]
  year_anom <- temp_year - temp_year_exp
  
  ##Creates an empty PDF
  pdf_name <- paste("../MORAclim_plots/Bayes_",exp_brick,".pdf",sep="")
  pdf(file=pdf_name,width=10,height=4,onefile=TRUE)
  
  for(i in 1:dim(meas_df)[1]){
    
    print(paste("Bayesian estimates for day",i))
    day_data <- data.frame(site=colnames(temp_year),
                           utmx_s=temp_meta_year$utmx_s,
                           utmy_s=temp_meta_year$utmy_s,
                           utmx=temp_meta_year$X_UTM,
                           utmy=temp_meta_year$Y_UTM,
                           elev=temp_meta_year$elev,
                           exp=as.numeric(temp_year_exp[i,]),
                           meas=as.numeric(temp_year[i,]),
                           anom=as.numeric(year_anom[i,]))
    day_data_comp <- day_data[!is.na(day_data$anom),]
    y <- day_data_comp$anom
    coords <- cbind(day_data_comp$utmx_s,day_data_comp$utmy_s)
    utmx_s <- scale(day_data_comp$utmx_s)
    utmx_scale <- attributes(utmx_s)$'scaled:scale'
    utmx_center <- attributes(utmx_s)$'scaled:center'
    utmy_s <- scale(day_data_comp$utmy_s)
    utmy_scale <- attributes(utmy_s)$'scaled:scale'
    utmy_center <- attributes(utmy_s)$'scaled:center'
    elev <- scale(day_data_comp$elev)
    elev_scale <- attributes(elev)$'scaled:scale'
    elev_center <- attributes(elev)$'scaled:center'
    
    ##Priors for Bayesian model.
    p <- 6
    starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1)
    
    tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
    
    priors.1 <- list("beta.Norm"=list(rep(0,p), diag(10,p)),
                     "phi.Unif"=c(0.1, 300), "sigma.sq.IG"=c(2, 2),
                     "tau.sq.IG"=c(2, 0.1))
    
    priors.2 <- list("beta.Flat", "phi.Unif"=c(3/1, 3/0.1),
                     "sigma.sq.IG"=c(2, 2), "tau.sq.IG"=c(2, 0.1))
    
    cov.model <- "exponential"
    
    n.samples <- nsamples * 80
    n.report <- 2000
    verbose <- TRUE
    
    if(pp==TRUE){
      m.1 <- spLM(y~utmx_s + utmy_s + elev + elev:utmx_s + elev:utmy_s, coords=coords,starting=starting,modified.pp=TRUE,
                  tuning=tuning, priors=priors.1, cov.model=cov.model,knots=c(10,10),
                  n.samples=n.samples, verbose=verbose, n.report=n.report)
    }else{
      m.1 <- spLM(y~utmx_s+utmy_s+elev+elev:utmx_s + elev:utmy_s, coords=coords,starting=starting,pp=TRUE,
                  tuning=tuning, priors=priors.1, cov.model=cov.model,
                  n.samples=n.samples, verbose=verbose, n.report=n.report)
    }
    
    
    burn.in <- 0.75*n.samples
    
    ##recover beta and spatial random effects
    #m.1.rc <- spRecover(m.1, start=burn.in, verbose=FALSE)
    #round(summary(m.1.rc$p.theta.recover.samples)$quantiles[,c(3,1,5)],4)
    #round(summary(m.1.rc$p.beta.recover.samples)$quantiles[,c(3,1,5)],4)
    #m.1.w.summary <- summary(mcmc(t(m.1.rc$p.w.recover.samples)))$quantiles[,c(3,1,5)]
    
    ##Predictions for new data locations.
    coords_pred <- xyFromCell(preds$utmx_s, cell=which(!is.na(preds$utmx_s[]))) / 1000
    int_pred <- rep(1,dim(coords_pred)[1])
    elev_pred <- (preds$elev[!is.na(preds$utmy_s[])] - elev_center) / elev_scale

    utmx_s_pred <- (preds$utmx_s[!is.na(preds$utmx_s[])] - utmx_center) / utmx_scale
    utmy_s_pred <- (preds$utmy_s[!is.na(preds$utmy_s[])] - utmy_center) / utmy_scale
    covars_pred <- cbind(int_pred,utmx_s_pred,utmy_s_pred,elev_pred,utmx_s_pred*elev_pred,utmy_s_pred*elev_pred)
    
    bayes_preds <- spPredict(m.1,pred.coords=coords_pred,pred.covars=covars_pred,
                             start=burn.in,end=n.samples,thin=20,n.report=10)
    predAnom <- rowMeans(bayes_preds$p.y.predictive.samples)
    predSD <- apply(bayes_preds$p.y.predictive.samples,FUN=sd,MARGIN = 1)
    predAnom_rast <- preds$utmx_s
    predAnom_rast[!is.na(preds$utmx_s)] <- predAnom
    predSD_rast <- preds$utmx_s
    predSD_rast[!is.na(preds$utmx_s)] <- predSD
    
    ##Resamples grid to 90m
    predAnom_rast_res <- disaggregate(predAnom_rast,fact=3,method='bilinear')
    predSD_rast_res <- disaggregate(predSD_rast,fact=3,method='bilinear')
    predEst_rast <- (temp_exp_year[[i]] / 10000) + predAnom_rast_res
    
    ##Makes a plot
    par(mfrow=c(1,3),oma=c(2,2,4,4))
    plot(predEst_rast,main=paste("Estimate, Day",i),zlim=c(-15,35),
         col=rev(rainbow(255)))
    plot(predAnom_rast,main=paste("Anomaly, Day",i),zlim=c(-10,10),
         col=colorRampPalette(c("blue", "white", "red"))(255))
    plot(predSD_rast,main=paste("Kriging Variance, Day",i),zlim=c(0,3),
         col=rev(heat.colors(255)))
    
    filename <- paste(out_prefix,"_day_",sprintf("%04d", (start_num - 1 + i)),sep="")
    predAnom_rast_res <- predAnom_rast_res * 10000
    predSD_rast_res <- predSD_rast_res * 10000
    predEst_rast <- predEst_rast * 10000
    writeRaster(predAnom_rast_res,filename=paste(filename,"_anom.tif",sep=""),
                datatype="INT4S",overwrite=overwrite)
    writeRaster(predSD_rast_res,filename=paste(filename,"_se.tif",sep=""),
                datatype="INT4S",overwrite=overwrite)
    writeRaster(predEst_rast,
                filename=paste(filename,"_est.tif",sep=""),
                datatype="INT4S",overwrite=overwrite)
  }
  ##Closes plot pdf
  dev.off()
}

krige_daily_anoms <- function(exp_brick,exp_brick_path,meas_df,meta,preds,out_prefix,
                              overwrite=FALSE){
  require(raster)
  require(automap)
  require(gstat)
  setwd(exp_brick_path)
  temp_exp_year <- brick(exp_brick)
  
  ##Checks inputs.
  stopifnot(nlayers(temp_exp_year) >= dim(meas_df)[1])
  #stopifnot(all(colnames(meas_df[,-1]) %in% meta$alt_code))
  stopifnot("utmx_s" %in% names(preds))
  stopifnot("utmy_s" %in% names(preds))
  stopifnot("elev" %in% names(preds))
  stopifnot(all(!is.nan(rowMeans(meas_df[,-1],na.rm=TRUE))))
 
  ##Extracts common days of years
  temp_year <- meas_df[,-1]
  start_num <- as.numeric(format(as.Date(meas_df$DATE[1]),format="%j"))
  end_num <- as.numeric(format(as.Date(meas_df$DATE[length(meas_df$DATE)]),format="%j"))
  temp_meta_year <- merge(data.frame(alt_code=colnames(temp_year)),
                          meta,by="alt_code",sort=FALSE)
  
  ##Extract values from raster
  print(paste("Extracting expected values from raster",exp_brick))
  temp_year_exp <- t(extract(temp_exp_year,cbind(temp_meta_year$X_UTM,temp_meta_year$Y_UTM),
                             method="bilinear",sp=FALSE) / 10000)
  temp_year_exp <- temp_year_exp[start_num:end_num,]
  year_anom <- temp_year - temp_year_exp
  
  ##Creates an empty PDF
  pdf_name <- paste("../MORAclim_plots/",exp_brick,".pdf",sep="")
  pdf(file=pdf_name,width=10,height=4,onefile=TRUE)
  
  for(i in 1:dim(meas_df)[1]){
    print(paste("Kriging estimates for day",i))
    day_data <- data.frame(site=colnames(temp_year),
                           utmx_s=temp_meta_year$utmx_s,
                           utmy_s=temp_meta_year$utmy_s,
                           utmx=temp_meta_year$X_UTM,
                           utmy=temp_meta_year$Y_UTM,
                           elev=temp_meta_year$elev,
                           coldair=temp_meta_year$cair_reg_9m,
                           meas=as.numeric(temp_year[i,]),
                           exp=as.numeric(temp_year_exp[i,]),
                           anom=as.numeric(year_anom[i,]))
    coordinates(day_data) <- ~utmx+utmy
    crs(day_data) <- pred_crs
    day_data_comp <- day_data[!is.na(day_data$anom),]
    day_data_cv <- autoKrige.cv(anom~utmx_s+utmy_s+elev,
                                model="Exp",input_data=day_data_comp)
    day_data_cv_res <- day_data_cv$krige.cv_output@data$zscore
    day_data_outl <- which(day_data_cv_res > 1.96 | day_data_cv_res < -1.96)
    if(length(day_data_outl)==0){
      day_data_red <- day_data_comp
    }else{
      day_data_red <- day_data_comp[-c(day_data_outl),]
    }
    day_data_var <- variogram(anom~utmx_s+utmy_s+I(elev/1000)+I(elev/1000):utmx_s+I(elev/1000):utmy_s,
                              width=50,cutoff=5000,data=day_data_red)
    day_data_mod <- autofitVariogram(anom~utmx_s+utmy_s+I(elev/1000)+I(elev/1000):utmx_s+I(elev/1000):utmy_s,
                                     cutoff=5000,width=50,cressie=TRUE,
                                     model="Exp",input_data=day_data_red)
    #par(mfrow=c(1,1))
    #plot(day_data_var,day_data_mod$var_model)
    day_gst <- gstat(NULL,"anom",formula=anom~utmx_s+utmy_s+I(elev/1000)+I(elev/1000):utmx_s+I(elev/1000):utmy_s,
                     data=day_data_comp,model=day_data_mod$var_model)
    day_anom <- interpolate(preds,day_gst,xyOnly=FALSE,na.rm=TRUE,index=1)
    day_se <- interpolate(preds,day_gst,xyOnly=FALSE,na.rm=TRUE,index=2)
    day_est <- (temp_exp_year[[i]] / 10000) + day_anom
    
    ##Makes a plot
    par(mfrow=c(1,3),oma=c(2,2,4,4))
    plot(day_est,main=paste("Estimate, Day",i),zlim=c(-15,35),
         col=rev(rainbow(255)))
    plot(day_anom,main=paste("Anomaly, Day",i),zlim=c(-10,10),
         col=colorRampPalette(c("blue", "white", "red"))(255))
    plot(day_se,main=paste("Kriging Variance, Day",i),zlim=c(0,3),
         col=rev(heat.colors(255)))
    
    filename <- paste(out_prefix,"_day_",sprintf("%04d", (start_num - 1 + i)),sep="")
    day_anom <- day_anom * 10000
    day_se <- day_se * 10000
    day_est <- day_est * 10000
    writeRaster(day_anom,filename=paste(filename,"_anom.tif",sep=""),
                datatype="INT4S",overwrite=overwrite)
    writeRaster(day_se,filename=paste(filename,"_se.tif",sep=""),
                datatype="INT4S",overwrite=overwrite)
    writeRaster(day_est,
                filename=paste(filename,"_est.tif",sep=""),
                datatype="INT4S",overwrite=overwrite)
  }
  ##Closes pdf
  dev.off()
}

krige_problem_anoms <- function(exp_brick,exp_brick_path,meas_df,probs,meta,preds,out_prefix,
                                overwrite=FALSE){
  require(raster)
  require(automap)
  require(gstat)
  setwd(exp_brick_path)
  temp_exp_year <- brick(exp_brick)
  
  ##Checks inputs.
  stopifnot(nlayers(temp_exp_year) >= dim(meas_df)[1])
  stopifnot(length(probs) == dim(meas_df)[1])
  stopifnot(all(colnames(meas_df[,-1]) %in% meta$alt_code))
  stopifnot("utmx_s" %in% names(preds))
  stopifnot("utmy_s" %in% names(preds))
  stopifnot("elev" %in% names(preds))
  stopifnot(all(!is.nan(rowMeans(meas_df[,-1],na.rm=TRUE))))
  
  ##Extracts common days of years
  temp_year <- meas_df[,-1]
  start_num <- as.numeric(format(as.Date(meas_df$DATE[1]),format="%j"))
  end_num <- as.numeric(format(as.Date(meas_df$DATE[length(meas_df$DATE)]),format="%j"))
  temp_meta_year <- merge(data.frame(alt_code=colnames(temp_year)),
                          meta,by="alt_code",sort=FALSE)
  
  ##Extract values from raster
  print(paste("Extracting expected values from raster",exp_brick))
  temp_year_exp <- t(extract(temp_exp_year,cbind(temp_meta_year$X_UTM,temp_meta_year$Y_UTM),
                             method="bilinear",sp=FALSE) / 10000)
  temp_year_exp <- temp_year_exp[start_num:end_num,]
  year_anom <- temp_year - temp_year_exp
  
  ##Creates an empty PDF
  #pdf_name <- paste("../MORAclim_plots/",exp_brick,".pdf",sep="")
  #pdf(file=pdf_name,width=10,height=4,onefile=TRUE)
  probs[is.na(probs)] <- FALSE
  
  for(i in 1:dim(meas_df)[1]){
    if(probs[i]==TRUE){
      print(paste("Kriging estimates for problem day",i))
      day_data <- data.frame(site=colnames(temp_year),
                             utmx_s=temp_meta_year$utmx_s,
                             utmy_s=temp_meta_year$utmy_s,
                             utmx=temp_meta_year$X_UTM,
                             utmy=temp_meta_year$Y_UTM,
                             elev=temp_meta_year$elev,
                             coldair=temp_meta_year$cair_reg_9m,
                             meas=as.numeric(temp_year[i,]),
                             exp=as.numeric(temp_year_exp[i,]),
                             anom=as.numeric(year_anom[i,]))
      coordinates(day_data) <- ~utmx+utmy
      crs(day_data) <- pred_crs
      day_data_comp <- day_data[!is.na(day_data$anom),]
      day_data_cv <- autoKrige.cv(anom~utmx_s+utmy_s+I(elev/1000)+I(elev/1000):utmx_s+I(elev/1000):utmy_s+I(elev/1000)^2+coldair,
                                  model="Exp",input_data=day_data_comp)
      day_data_cv_res <- day_data_cv$krige.cv_output@data$zscore
      day_data_outl <- which(day_data_cv_res > 1.96 | day_data_cv_res < -1.96)
      if(length(day_data_outl)==0){
        day_data_red <- day_data_comp
      }else{
        day_data_red <- day_data_comp[-c(day_data_outl),]
      }
      day_data_var <- variogram(anom~utmx_s+utmy_s+I(elev/1000)+I(elev/1000)^2+coldair,
                                width=50,cutoff=5000,data=day_data_red)
      day_data_mod <- autofitVariogram(anom~utmx_s+utmy_s+I(elev/1000)+I(elev/1000):utmx_s+I(elev/1000):utmy_s+I(elev/1000)^2+coldair,
                                       cutoff=5000,width=50,cressie=TRUE,
                                       model="Exp",input_data=day_data_red)
      #par(mfrow=c(1,1))
      #plot(day_data_var,day_data_mod$var_model)
      day_gst <- gstat(NULL,"anom",formula=anom~utmx_s+utmy_s+I(elev/1000)+I(elev/1000):utmx_s+I(elev/1000):utmy_s+I(elev/1000)^2+coldair,
                       data=day_data_red,model=day_data_mod$var_model)
      day_anom <- interpolate(preds,day_gst,xyOnly=FALSE,na.rm=TRUE,index=1)
      day_se <- interpolate(preds,day_gst,xyOnly=FALSE,na.rm=TRUE,index=2)
      day_est <- (temp_exp_year[[i]] / 10000) + day_anom
      
      ##Makes a plot
      par(mfrow=c(1,3),oma=c(2,2,4,4))
      plot(day_est,main=paste("Estimate, Day",i),zlim=c(-15,35),
           col=rev(rainbow(255)))
      plot(day_anom,main=paste("Anomaly, Day",i),zlim=c(-10,10),
           col=colorRampPalette(c("blue", "white", "red"))(255))
      plot(day_se,main=paste("Kriging Variance, Day",i),zlim=c(0,3),
           col=rev(heat.colors(255)))
      
      filename <- paste(out_prefix,"_day_",sprintf("%04d", (start_num - 1 + i)),sep="")
      day_anom <- day_anom * 10000
      day_se <- day_se * 10000
      day_est <- day_est * 10000
      writeRaster(day_anom,filename=paste(filename,"_anom_b.tif",sep=""),
                  datatype="INT4S",overwrite=overwrite)
      writeRaster(day_se,filename=paste(filename,"_se_b.tif",sep=""),
                  datatype="INT4S",overwrite=overwrite)
      writeRaster(day_est,
                  filename=paste(filename,"_est_b.tif",sep=""),
                  datatype="INT4S",overwrite=overwrite)
    }else{
      print("Not a problem day, skipping...")
    }
  }
    
  ##Closes pdf
  #dev.off()
}


##Function to remove anomalous spikes in temperature data and fill with interpolated values.
remove_spikes <- function(ts,ts_lag=2,thresh=c(-5,5)){
  ts_l1 <- lag(ts,k=ts_lag)
  ts_l2 <- lag(ts,k=ts_lag*-1)
  ts_delta1 <- ts - ts_l1
  ts_delta2 <- ts - ts_l2
  ts_spikes <- as.logical(ts_delta1 > thresh[2] | ts_delta1 < thresh[1] |
                             ts_delta2 > thresh[2] | ts_delta2 < thresh[1])
  ts_range <- as.logical(ts > 35 |ts < -25)
  ts[ts_spikes] <- NA
  ts[ts_range] <- NA
  ts_filled <- na.spline(ts,maxgap=ts_lag*4,na.rm=FALSE)
  return(ts_filled)
}

##Function to interpolate seasonal estimates of microclimate parameters.
interpolateMicro <- function(ts_frame){
  interp_frame <- ts_frame
  interp_frame[,3:dim(ts_frame)[2]] <- NA
  for(i in 3:dim(ts_frame)[2]){
    series <- xts(ts_frame[,i],order.by=ts_frame$date)
    if(all(is.na(series))){
      interp_frame[,i] <- rep(NA,length(series))
    }else{
      interp_frame[,i] <- as.numeric(na.spline(series,na.rm=FALSE))
    }
  }
  return(interp_frame)
}
