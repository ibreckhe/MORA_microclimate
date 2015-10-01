##Functions for Mt. Rainier microclimate analysis.

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
  compare_cols <- colSums(local_frame[1:length(regional_series),-1],na.rm=TRUE) != 0.00 &
    colSums(!is.na(local_frame[1:length(regional_series),-1]),na.rm=TRUE) > 200
  
  for (i in which(compare_cols)){
    title <- site_names[i]
    local_meas <- local_frame[1:length(regional_series),(i+1)]
    if(compare_cols[i]){
      
      ##Puts complete data into a data frame.
      comp_dat <- data.frame(regional_series,local_meas)
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
plot.beta.ts <-  function(y.dates,beta.0,beta.1,beta.2,beta.3,beta.4){
  par(mfrow=c(5,1),mar=c(0,4,0,0),oma=c(4,1,1,1))
  plot(y.dates, beta.0[1,], pch=19, cex=0.5, xlab="", ylab="Intercept",
       xaxt="n",ylim=range(beta.0))
  arrows(y.dates, beta.0[1,], y.dates, beta.0[3,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  arrows(y.dates, beta.0[1,], y.dates, beta.0[2,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  abline(h=0,lty=2,col=2)
  xlabels <- seq(as.POSIXct("2008-10-1 00:00:00"),as.POSIXct("2014-10-1 00:00:00"),by='month')
  xlabels2 <- seq(as.POSIXct("2009-1-1 00:00:00"),as.POSIXct("2014-1-1 00:00:00"),by='year')
  axis.POSIXct(side=1, at=xlabels,format="%b",labels = FALSE)
  axis.POSIXct(side=1,at=xlabels2,format="%Y",tck=1,labels = FALSE,hadj=-3,padj=1.8)
  
  plot(y.dates, beta.1[1,], pch=19, cex=0.5, xlab="", ylab="Elevation",
       xaxt="n",ylim=c(-4,4))
  arrows(y.dates, beta.1[1,], y.dates, beta.1[3,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  arrows(y.dates, beta.1[1,], y.dates, beta.1[2,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  abline(h=0,lty=2,col=2)
  xlabels <- seq(as.POSIXct("2008-10-1 00:00:00"),as.POSIXct("2014-10-1 00:00:00"),by='month')
  xlabels2 <- seq(as.POSIXct("2009-1-1 00:00:00"),as.POSIXct("2014-1-1 00:00:00"),by='year')
  axis.POSIXct(side=1, at=xlabels,format="%b",labels = FALSE)
  axis.POSIXct(side=1,at=xlabels2,format="%Y",tck=1,labels = FALSE,hadj=-3,padj=1.8)
  
  plot(y.dates, beta.2[1,], pch=19, cex=0.5, xlab="", ylab="Cold-air",
       xaxt="n",ylim=c(-4,4))
  arrows(y.dates, beta.2[1,], y.dates, beta.2[3,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  arrows(y.dates, beta.2[1,], y.dates, beta.2[2,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  abline(h=0,lty=2,col=2)
  xlabels <- seq(as.POSIXct("2008-10-1 00:00:00"),as.POSIXct("2014-10-1 00:00:00"),by='month')
  xlabels2 <- seq(as.POSIXct("2009-1-1 00:00:00"),as.POSIXct("2014-1-1 00:00:00"),by='year')
  axis.POSIXct(side=1, at=xlabels,format="%b",labels = FALSE)
  axis.POSIXct(side=1,at=xlabels2,format="%Y",tck=1,labels = FALSE,hadj=-3,padj=1.8)
  
  plot(y.dates, beta.3[1,], pch=19, cex=0.5, xlab="", ylab="Elev:Cold-air",
       xaxt="n",ylim=c(-4,4))
  arrows(y.dates, beta.3[1,], y.dates, beta.3[3,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  arrows(y.dates, beta.3[1,], y.dates, beta.3[2,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  abline(h=0,lty=2,col=2)
  xlabels <- seq(as.POSIXct("2008-10-1 00:00:00"),as.POSIXct("2014-10-1 00:00:00"),by='month')
  xlabels2 <- seq(as.POSIXct("2009-1-1 00:00:00"),as.POSIXct("2014-1-1 00:00:00"),by='year')
  axis.POSIXct(side=1, at=xlabels,format="%b",labels = TRUE)
  axis.POSIXct(side=1,at=xlabels2,format="%Y",tck=1,labels = TRUE,hadj=-3,padj=1.8)
  
  plot(y.dates, beta.4[1,], pch=19, cex=0.5, xlab="", ylab="Canopy Volume",
       xaxt="n",ylim=c(-4,4))
  arrows(y.dates, beta.4[1,], y.dates, beta.4[3,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  arrows(y.dates, beta.4[1,], y.dates, beta.4[2,], length=0.02, angle=90,col=rgb(0.25,0.25,0.25,0.25))
  abline(h=0,lty=2,col=2)
  xlabels <- seq(as.POSIXct("2008-10-1 00:00:00"),as.POSIXct("2014-10-1 00:00:00"),by='month')
  xlabels2 <- seq(as.POSIXct("2009-1-1 00:00:00"),as.POSIXct("2014-1-1 00:00:00"),by='year')
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