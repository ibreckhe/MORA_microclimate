##Script to estimate optimum MAT for MORA's tree species using a large network of forest plots.
##Author: Ian Breckheimer
##Date: November 6th, 2014.

##Sets up workspace and brings in data.
library(rjags)
library(ggmcmc)
library(reshape2)

##Sets working directory.
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/")
source("~/code/MORA_microclimate/community_optimum_functions.R")

##Loads data.
data <- read.csv("USFS_PMR_Franklin_YNdata_ClimateYS.csv")

##Filters out the franklin data.
data <- data[data$DATASET != "Franklin",]

##Fits models for each species using glm.
abam_glm <- glm(ABAM_YN~MAT+I(MAT^2),data=data,family="binomial")
abla_glm <- glm(ABLA_YN~MAT+I(MAT^2),data=data,family="binomial")
tsme_glm <- glm(TSME_YN~MAT+I(MAT^2),data=data,family="binomial")
tshe_glm <- glm(TSHE_YN~MAT+I(MAT^2),data=data,family="binomial")
cano_glm <- glm(CANO_YN~MAT+I(MAT^2),data=data,family="binomial")
psme_glm <- glm(PSME_YN~MAT+I(MAT^2),data=data,family="binomial")

##Predicts based on new data.
newdata <- data.frame(MAT=seq(-10,20,by=0.01))
abam_pred <- predict(abam_glm,newdata=newdata,type="response")
abla_pred <- predict(abla_glm,newdata=newdata,type="response")
tsme_pred <- predict(tsme_glm,newdata=newdata,type="response")
tshe_pred <- predict(tshe_glm,newdata=newdata,type="response")
cano_pred <- predict(cano_glm,newdata=newdata,type="response")
psme_pred <- predict(psme_glm,newdata=newdata,type="response")

##Plot predictions
par(mfrow=c(2,3),mar=c(2,2.5,2.5,0),oma=c(4,4,1,1))
plot(jitter(ABLA_YN,amount=0.05)~MAT,data=data,pch=".",xlim=c(-10,25),
     col=rgb(0,0,0,0.05),main="Abies lasiocarpa",xlab="",ylab="")
points(newdata$MAT,abla_pred,type="l",col="green")
plot(jitter(TSME_YN,amount=0.05)~MAT,data=data,pch=".",xlim=c(-10,25),
     col=rgb(0,0,0,0.05),main="Tsuga mertensiana",xlab="",ylab="")
points(newdata$MAT,tsme_pred,type="l",col="green")
plot(jitter(CANO_YN,amount=0.05)~MAT,data=data,pch=".",xlim=c(-10,25),
     col=rgb(0,0,0,0.05),main="Callitropsis nootkatensis",xlab="",ylab="")
points(newdata$MAT,cano_pred,type="l",col="green")
plot(jitter(ABAM_YN,amount=0.05)~MAT,data=data,pch=".",xlim=c(-10,25),
     col=rgb(0,0,0,0.05),main="Abies amabilis",xlab="",ylab="")
points(newdata$MAT,abam_pred,type="l",col="green")
plot(jitter(TSHE_YN,amount=0.05)~MAT,data=data,pch=".",xlim=c(-10,25),
     col=rgb(0,0,0,0.05),main="Tsuga heterophylla",xlab="",ylab="")
points(newdata$MAT,tshe_pred,type="l",col="green")
plot(jitter(PSME_YN,amount=0.05)~MAT,data=data,pch=".",xlim=c(-10,25),
     col=rgb(0,0,0,0.05),main="Pseudotsuga menziesii",xlab="",ylab="")
points(newdata$MAT,psme_pred,type="l",col="green")
mtext(text="Presence Probability",side=2,padj=-1,adj=0.5,outer=TRUE)
mtext(text="Mean Annual Temperature (C)",side=1,padj=1,adj=0.5,outer=TRUE)

##Subsets data and gets a subset in a long format.
data_s <- data[sample(1:dim(data)[1],size=1000,
                      replace=FALSE),c("REGION","AGENCY","PLOTNBR",
                                      "PLOTTYPE","DATASET","UTMEAST",
                                      "Lat","Long","ABAM_YN","ABLA_YN",
                                      "CANO_YN","PSME_YN","TSHE_YN","TSME_YN",
                                      "MAT","CMD")]
data_s$PLOT_UNIQUE <- 1:dim(data_s)[1]
data_l <- melt(data_s,id.vars=c("PLOT_UNIQUE","REGION","AGENCY","PLOTNBR","PLOTTYPE","DATASET",
                                "UTMEAST","Lat","Long","MAT","CMD"))

##Preps data for JAGS.
x <- scale(data_l$MAT)
y <- data_l$value
species <- as.numeric(data_l$variable)
plot <- data_l$PLOT_UNIQUE

jagsmodel <- fit.jags.mixed.tree(x=x,y=y,species=species,plot=plot,nsamples=10000)

##Saves posterior samples.
save(jagsmodel,file="../results/6spp_temp_curves_posterior_samples.Rdata")
load(file="../results/6spp_temp_curves_posterior_samples.Rdata")

##Computes credible intervals for parameters.
jagsparams <- get.param.creds(jagsmodel$out,params=c("width.g[1]",
                                                     "width.g[2]",
                                                     "width.g[3]",
                                                     "width.g[4]",
                                                     "width.g[5]",
                                                     "width.g[6]",
                                                     "opt.g[1]",
                                                     "opt.g[2]",
                                                     "opt.g[3]",
                                                     "opt.g[4]",
                                                     "opt.g[5]",
                                                     "opt.g[6]",
                                                     "height.g[1]",
                                                     "height.g[2]",
                                                     "height.g[3]",
                                                     "height.g[4]",
                                                     "height.g[5]",
                                                     "height.g[6]"))


##Plot predictions
new_x <- seq(-7,5,by=0.01)
new_x_scale <- attr(x,"scaled:scale")
new_x_center <- attr(x,"scaled:center")
new_x_unscaled <- new_x * new_x_scale + new_x_center
par(mfrow=c(2,3),mar=c(2,2.5,2.5,0),oma=c(4,4,1,1))
plot(jitter(ABLA_YN,amount=0.05)~MAT,data=data,pch=".",xlim=c(-10,25),
     col=rgb(0,0,0,0.1),main="Abies lasiocarpa",xlab="",ylab="")
points(new_x_unscaled,abla_fun(new_x),type="l",col="green")
plot(jitter(TSME_YN,amount=0.05)~MAT,data=data,pch=".",xlim=c(-10,25),
     col=rgb(0,0,0,0.1),main="Tsuga mertensiana",xlab="",ylab="")
points(new_x_unscaled,tsme_fun(new_x),type="l",col="green")
plot(jitter(CANO_YN,amount=0.05)~MAT,data=data,pch=".",xlim=c(-10,25),
     col=rgb(0,0,0,0.1),main="Callitropsis nootkatensis",xlab="",ylab="")
points(new_x_unscaled,cano_fun(new_x),type="l",col="green")
plot(jitter(ABAM_YN,amount=0.05)~MAT,data=data,pch=".",xlim=c(-10,25),
     col=rgb(0,0,0,0.1),main="Abies amabilis",xlab="",ylab="")
points(new_x_unscaled,abam_fun(new_x),type="l",col="green")
plot(jitter(TSHE_YN,amount=0.05)~MAT,data=data,pch=".",xlim=c(-10,25),
     col=rgb(0,0,0,0.1),main="Tsuga heterophylla",xlab="",ylab="")
points(new_x_unscaled,tshe_fun(new_x),type="l",col="green")
plot(jitter(PSME_YN,amount=0.05)~MAT,data=data,pch=".",xlim=c(-10,25),
     col=rgb(0,0,0,0.1),main="Pseudotsuga menziesii",xlab="",ylab="")
points(new_x_unscaled,psme_fun(new_x),type="l",col="green")
mtext(text="Presence Probability",side=2,padj=-1,adj=0.5,outer=TRUE)
mtext(text="Mean Annual Temperature (C)",side=1,padj=1,adj=0.5,outer=TRUE)

##Combo functions.
combo_fun <- function(x, pres_vec=rep(FALSE,6)) {
  nx <- length(x)
  (ifelse(rep(pres_vec[1],nx),abla_fun(x),1-abla_fun(x)) *
  ifelse(rep(pres_vec[2],nx),tsme_fun(x),1-tsme_fun(x)) *
  ifelse(rep(pres_vec[3],nx),cano_fun(x),1-cano_fun(x)) *
  ifelse(rep(pres_vec[4],nx),abam_fun(x),1-abam_fun(x)) *
  ifelse(rep(pres_vec[5],nx),tshe_fun(x),1-tshe_fun(x)) *
  ifelse(rep(pres_vec[6],nx),psme_fun(x),1-psme_fun(x)))
}

combo_pred0 <- combo_fun(new_x,pres_vec=c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE))
combo_pred1 <- combo_fun(new_x,pres_vec=c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE))
combo_pred2 <- combo_fun(new_x,pres_vec=c(FALSE,TRUE,TRUE,FALSE,FALSE,FALSE))
combo_pred3 <- combo_fun(new_x,pres_vec=c(FALSE,FALSE,TRUE,TRUE,FALSE,FALSE))
combo_pred4 <- combo_fun(new_x,pres_vec=c(FALSE,FALSE,FALSE,TRUE,TRUE,FALSE))
combo_pred5 <- combo_fun(new_x,pres_vec=c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE))
combo_pred6 <- combo_fun(new_x,pres_vec=c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE))

pal <- rainbow(7,v=0.8)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(1,1,1,1))
plot(new_x_unscaled,combo_pred0,type="l",ylim=c(0,0.8),xlim=c(-5,20),
     col=pal[1],xlab="Mean Annual Temperature (C)",ylab="Joint Probability")
points(new_x_unscaled,combo_pred1,type="l",col=pal[2])
points(new_x_unscaled,combo_pred2,type="l",col=pal[3])
points(new_x_unscaled,combo_pred3,type="l",col=pal[4])
points(new_x_unscaled,combo_pred4,type="l",col=pal[5])
points(new_x_unscaled,combo_pred5,type="l",col=pal[6])
points(new_x_unscaled,combo_pred6,type="l",col=pal[7])
legend("topleft",legend=c("ABLA","ABLA-TSME","TSME-CANO","CANO-ABAM","ABAM-TSHE","TSHE-PSME","PSME"),
       ncol=2,lty=1,col=pal,bty="n")

##Integrates prediction functions.
combo_dens <- function(x,pres_vec) { 
  y <- combo_fun(x,pres_vec)
  yi <- integrate(combo_fun, -Inf, +Inf, pres_vec=pres_vec)
  return(y/yi[[1]])
}

combo_dens0 <- combo_dens(new_x,pres_vec=c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE))
combo_dens1 <- combo_dens(new_x,pres_vec=c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE))
combo_dens2 <- combo_dens(new_x,pres_vec=c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE))
combo_dens3 <- combo_dens(new_x,pres_vec=c(FALSE,TRUE,TRUE,TRUE,FALSE,FALSE))
combo_dens4 <- combo_dens(new_x,pres_vec=c(FALSE,FALSE,TRUE,TRUE,TRUE,FALSE))
combo_dens5 <- combo_dens(new_x,pres_vec=c(FALSE,FALSE,FALSE,TRUE,TRUE,TRUE))
combo_dens6 <- combo_dens(new_x,pres_vec=c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE))
combo_dens7 <- combo_dens(new_x,pres_vec=c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE))


pdf("../results/3spp_densities.pdf",width=8,height=4.5)
pal <- rev(rainbow(8,v=0.8))
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(1,1,1,1))
plot(new_x_unscaled,combo_dens0,type="l",ylim=c(-0.03,1.5),xlim=c(-5,18),lwd=1.5,
     col=pal[1],xlab=expression(Mean~Annual~Temperature~(degree*C),sep=""),
     ylab="Probability Density")
points(new_x_unscaled,combo_dens1,type="l",col=pal[2],lwd=1.5)
points(new_x_unscaled,combo_dens2,type="l",col=pal[3],lwd=1.5)
points(new_x_unscaled,combo_dens3,type="l",col=pal[4],lwd=1.5)
points(new_x_unscaled,combo_dens4,type="l",col=pal[5],lwd=1.5)
points(new_x_unscaled,combo_dens5,type="l",col=pal[6],lwd=1.5)
points(new_x_unscaled,combo_dens6,type="l",col=pal[7],lwd=1.5)
points(new_x_unscaled,combo_dens7,type="l",col=pal[8],lwd=1.5)
points(jitter(data_l$MAT,amount=0.1),y=jitter(rep(-0.03,length(data_l$MAT)),amount=0.03),pch=".",col=rgb(0,0,0,0.5))
legend("topleft",legend=c("ABLA","ABLA-TSME","ABLA-TSME-CANO","TSME-CANO-ABAM"),
       ncol=1,lty=1,col=pal[1:4],bty="n",cex=1,y.intersp = 1.1,lwd=1.5)
legend("topright",legend=c("CANO-ABAM-TSHE","ABAM-TSHE-PSME","TSHE-PSME","PSME"),
       ncol=1,lty=1,col=pal[5:8],bty="n",cex=1,y.intersp = 1.1,lwd=1.5)
dev.off()

