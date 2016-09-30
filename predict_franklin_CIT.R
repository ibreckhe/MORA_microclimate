##Script to predict community inferred temperature based on spatial covariates.
##Author: Ian Breckheimer
##Date: 3 December 2014

####Sets up workspace and loads data. ####
require(rjags)
require(ggmcmc)
source('~/code/MORA_microclimate/community_optimum_functions.R')
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/")

CIT_samps <- read.csv("Franklindata_CIT_samples.csv",strip.white=TRUE,
                      stringsAsFactors=FALSE)
CIT_samps$X <- gsub("\\s{2,}"," ",CIT_samps$X)
CIT_samps_all <- read.csv("Franklindata_CIT_samples_all.csv",strip.white=TRUE,
                          stringsAsFactors=FALSE)
CIT_samps_all$X <- gsub("\\s{2,}"," ",CIT_samps_all$X)

fsites <- read.csv("franklindata_utm_covars.csv",strip.white = TRUE,
                   stringsAsFactors=FALSE)
fsites$DESC[fsites$DESC==""] <- NA
fsites$DESC <- as.factor(as.character(fsites$DESC))
fsites$PLOT <- as.factor(paste(fsites$SITE_NAME,fsites$PLOT_NO,sep=" "))


####Munges data ####
quantfun <- function(x) {c(CIT_mean=mean(x),
                           CIT_sd=sd(x),
                           CIT_lwr=quantile(x,probs=0.025),
                           CIT_lwr_quart=quantile(x,probs=0.25),
                           CIT_med=quantile(x,probs=0.5),
                           CIT_upr_quart=quantile(x,probs=0.75),
                           CIT_upr=quantile(x,probs=0.975))}
CIT_quants <- data.frame(PLOT=CIT_samps$X,t(apply(CIT_samps[,-1],FUN=quantfun,MARGIN=1)))
CIT_quants_all <- data.frame(PLOT=CIT_samps_all$X,t(apply(CIT_samps_all[,-1],FUN=quantfun,MARGIN=1)))
avg_sd <- colMeans(CIT_quants[,-1])[2]
avg_sd_all <- colMeans(CIT_quants_all[,-1])[2]
mean_sd <- sd(CIT_quants$mean)
mean_sd_all <- sd(CIT_quants_all$mean)
fsites_CIT <- merge(fsites,CIT_quants)
fsites_CIT_all <- merge(fsites,CIT_quants_all)

##Finds sites that do not match.
fsites_nomatch <- unique(CIT_samps$X)[! unique(CIT_samps$X) %in% unique(fsites$PLOT)]

##Removes outlier solar radiation.
fsites_CIT$MORA_srad_yearsum_9m[fsites_CIT$MORA_srad_yearsum_9m>5e08] <- NA
fsites_CIT_all$MORA_srad_yearsum_9m[fsites_CIT_all$MORA_srad_yearsum_9m>5e08] <- NA

####Plots data ####
par(mfrow=c(1,2))
cols <- c("slateblue","orange","grey60","green")
plot(CIT_mean~MORA_elev_3m,data=fsites_CIT,col=cols[fsites_CIT$DESC],xlab="Elevation (m)",
     ylab="Community Inferred Temp. (C)",pch=20)
arrows(x0=fsites_CIT$MORA_elev_3m,x1=fsites_CIT$MORA_elev_3m,
       y0=fsites_CIT$CIT_lwr_quart.25.,y1=fsites_CIT$CIT_upr_quart.75.,
       code=3,angle=90,length=0,col=cols[fsites_CIT$DESC],lwd=0.5)
legend("bottomleft",title="Site Desc.",bty="n",legend=levels(fsites$DESC),
       pch=20,col=cols)
plot(CIT_mean~MORA_elev_3m,data=fsites_CIT_all,col=cols[fsites_CIT_all$DESC],xlab="Elevation (m)",
     ylab="Community Inferred Temp. (C)",pch=20)
arrows(x0=fsites_CIT_all$MORA_elev_3m,x1=fsites_CIT_all$MORA_elev_3m,
       y0=fsites_CIT_all$CIT_lwr_quart.25.,y1=fsites_CIT_all$CIT_upr_quart.75.,
       code=3,angle=90,length=0,col=cols[fsites_CIT$DESC],lwd=0.5)

##Plots
par(mfrow=c(1,3),mar=c(2,2,2,0),oma=c(2,2,0,1))
plot(fsites_CIT$CIT_mean~fsites_CIT$MAT,main="ClimateWNA",xlim=c(1,9),ylim=c(1,9),
     xlab="",ylab="")
abline(0,1,lty=2)
plot(fsites_CIT$CIT_mean~fsites_CIT$MAT_c,main="PRISM / WNA",xlim=c(1,9),ylim=c(1,9),
     xlab="",ylab="",yaxt="n")
axis(2,labels=FALSE)
abline(0,1,lty=2)
plot(fsites_CIT_all$CIT_mean~fsites_CIT$MAT_c,main="PRISM",xlim=c(1,9),ylim=c(1,9),
     xlab="",ylab="",yaxt="n")
axis(2,labels=FALSE)
abline(0,1,lty=2)
mtext(side=1,outer=TRUE,text="Interpolated MAT estimate (C)",padj=1)
mtext(side=2,outer=TRUE,text="Community Inferred MAT estimate (C)",padj=-1)

####Builds a predictive model incorporating error in CIT and elevation.####

## Preps data using ClimateWNA CIT
##Complete cases
fsites_CIT <- fsites_CIT[complete.cases(fsites_CIT),]

##Scales elevation and CIT samples.
cit_samp_mean <- mean(as.matrix(CIT_samps[,-1]))
cit_samp_sd <- sd(as.matrix(CIT_samps[,-1]))
cit_samp_scaled <- (CIT_samps[,-1] - cit_samp_mean) / cit_samp_sd
cit_samp_avg_sd <- mean(apply(cit_samp_scaled,FUN=sd,MARGIN=1))

desc <- as.numeric(fsites_CIT$DESC)
elev <- scale(fsites_CIT$MORA_elev_3m)
elev_cat <- as.factor(elev>0.0)
cit <- scale(fsites_CIT$CIT_mean)
cair <- scale(fsites_CIT$MORA_coldair_index)
cair_cat <- cut(cair,breaks=c(-5,-0.5,0.5,5),include.lowest=TRUE)
srad <- scale(fsites_CIT$MORA_srad_yearsum_9m)
srad_cat <- as.factor(srad>0)
relev <- scale(fsites_CIT$MORA_elev_3m - fsites_CIT$MORA_elev__focal_729m)
dry <- scale(fsites_CIT$MORA_dry_index_81m)
dry_cat <- cut(dry,breaks=c(-5,-0.5,0.5,5),include.lowest=TRUE)
topo <- scale(fsites_CIT$MORA_topoidx_81m_log)
topo_cat <- as.factor(topo >= 0)
snow <- scale(fsites_CIT$MORA_snow_pred_2013_3m_parkbound_est)
snow_cat <- as.factor(snow > 1)
dry_snow_cat <- factor(paste(snow_cat,dry_cat),
                       labels=c("Rainy-Int","Rainy-Dry","Rainy-Wet","Snowy-Int","Snowy-Dry","Snowy-Wet"))
levels(dry_snow_cat) <- c("Rainy-Int","Rainy-Dry","Rainy-Wet","Snowy","Snowy","Snowy")
dry_snow_cat <- factor(dry_snow_cat,levels(dry_snow_cat)[c(4,2,1,3)])
dry_snow_num <- as.numeric(dry_snow_cat)
cair_snow_cat <- factor(paste(snow_cat,cair_cat),
                            labels=c("Rainy-Int","Rainy-Cold","Rainy-Warm","Snowy-Int","Snowy-Cold","Snowy-Warm"))
levels(cair_snow_cat) <- c("Rainy-Int","Rainy-Cold","Rainy-Warm","Snowy","Snowy","Snowy")
cair_snow_cat <- factor(cair_snow_cat,levels(cair_snow_cat)[c(4,2,1,3)])
cair_snow_num <- as.numeric(cair_snow_cat)
utmx <- scale(fsites_CIT$utmx)
n <- length(cit)

# Parameters for measurement error model
obs_tauy <- 1 / (cit_samp_avg_sd  * cit_samp_avg_sd)
obs_taux <- 1 / (0.2 * 0.2)
pop_sdy <- sd(cit,na.rm=TRUE)
pop_sdx <- sd(elev,na.rm=TRUE)

##Preps data using PRISM CIT
fsites_CIT_all <- fsites_CIT_all[complete.cases(fsites_CIT_all),]

##Scales elevation and CIT samples.
cit_samp_all_mean <- mean(as.matrix(CIT_samps_all[,-1]))
cit_samp_all_sd <- sd(as.matrix(CIT_samps_all[,-1]))
cit_samp_all_scaled <- (CIT_samps_all[,-1] - cit_samp_all_mean) / cit_samp_all_sd
cit_samp_all_avg_sd <- mean(apply(cit_samp_all_scaled,FUN=sd,MARGIN=1))

desc_all <- as.numeric(fsites_CIT_all$DESC)
elev_all <- scale(fsites_CIT_all$MORA_elev_3m)
elev_cat_all <- as.factor(elev_all>0.0)
cit_all <- scale(fsites_CIT_all$CIT_mean)
cair_all <- scale(fsites_CIT_all$MORA_coldair_index)
cair_cat_all <- cut(cair_all,breaks=c(-5,-0.5,0.5,5),include.lowest=TRUE)
srad_all <- scale(fsites_CIT_all$MORA_srad_yearsum_9m)
srad_cat_all <- as.factor(srad_all>0)
relev_all <- scale(fsites_CIT_all$MORA_elev_3m - fsites_CIT_all$MORA_elev__focal_729m)
dry_all <- scale(fsites_CIT_all$MORA_dry_index_81m)
dry_cat_all <- cut(dry_all,breaks=c(-5,-0.5,0.5,5),include.lowest=TRUE)
topo_all <- scale(fsites_CIT_all$MORA_topoidx_81m_log)
topo_cat_all <- as.factor(topo_all >= 0)
snow_all <- scale(fsites_CIT_all$MORA_snow_pred_2013_3m_parkbound_est)
snow_cat_all <- as.factor(snow_all > 1)
dry_snow_cat_all <- factor(paste(snow_cat_all,dry_cat_all),
                       labels=c("Rainy-Int","Rainy-Dry","Rainy-Wet","Snowy-Int","Snowy-Dry","Snowy-Wet"))
levels(dry_snow_cat_all) <- c("Rainy-Int","Rainy-Dry","Rainy-Wet","Snowy","Snowy","Snowy")
dry_snow_cat_all <- factor(dry_snow_cat_all,levels(dry_snow_cat_all)[c(4,2,1,3)])
dry_snow_num_all <- as.numeric(dry_snow_cat_all)
cair_snow_cat_all <- factor(paste(snow_cat_all,cair_cat_all),
                           labels=c("Rainy-Int","Rainy-Cold","Rainy-Warm","Snowy-Int","Snowy-Cold","Snowy-Warm"))
levels(cair_snow_cat_all) <- c("Rainy-Int","Rainy-Cold","Rainy-Warm","Snowy","Snowy","Snowy")
cair_snow_cat_all <- factor(cair_snow_cat_all,levels(cair_snow_cat_all)[c(4,2,1,3)])
cair_snow_num_all <- as.numeric(cair_snow_cat_all)
utmx_all <- scale(fsites_CIT_all$utmx)
n_all <- length(cit_all)

# Parameters for measurement error model
obs_tauy_all <- 1 / (cit_samp_all_avg_sd  * cit_samp_all_avg_sd)
obs_taux_all <- 1 / (0.2 * 0.2)
pop_sdy_all <- sd(cit_all,na.rm=TRUE)
pop_sdx_all <- sd(elev_all,na.rm=TRUE)

####Plots data with alternative factor.
cols <- c("blue","green","grey60","orange")
par(mfrow=c(1,2))
plot(CIT_mean~MORA_elev_3m,data=fsites_CIT,col=cols[cair_snow_cat],xlab="Elevation (m)",
     ylab="Community Inferred Temp. (C)",pch=20,ylim=c(3,9))
arrows(x0=fsites_CIT$MORA_elev_3m,x1=fsites_CIT$MORA_elev_3m,
       y0=fsites_CIT$CIT_lwr_quart.25.,y1=fsites_CIT$CIT_upr_quart.75.,
       code=3,angle=90,length=0,col=cols[cair_snow_cat],lwd=0.5)
abline(lm(fsites_CIT$MAT~fsites_CIT$MORA_elev_3m),lty=2)
legend("topright",title="Site Class",bty="n",legend=levels(cair_snow_cat),
       pch=20,col=cols)
plot(CIT_mean~MORA_elev_3m,data=fsites_CIT_all,col=cols[cair_snow_cat_all],xlab="Elevation (m)",
     ylab="Community Inferred Temp. (C)",pch=20,ylim=c(3,9))
arrows(x0=fsites_CIT_all$MORA_elev_3m,x1=fsites_CIT_all$MORA_elev_3m,
       y0=fsites_CIT_all$CIT_lwr_quart.25.,y1=fsites_CIT_all$CIT_upr_quart.75.,
       code=3,angle=90,length=0,col=cols[cair_snow_cat_all],lwd=0.5)
abline(lm(fsites_CIT_all$MAT_c~fsites_CIT$MORA_elev_3m),lty=2)

#### linear model for comparison.
fsites_lm <- lm(cit~elev_all*cair+elev_all*snow_cat)
fsites_lm_all <- lm(cit_all~elev_all*cair_all)
summary(fsites_lm)
summary(fsites_lm_all)
AIC(fsites_lm)
AIC(fsites_lm_all)

#### Model with measurement error in both X and y. ####

# specify model
cat("
    model {
    ## Priors
    alpha[1] ~ dnorm(0, .001)
    alpha[2] ~ dnorm(0, .001)
    alpha[3] ~ dnorm(0, .001)
    alpha[4] ~ dnorm(0, .001)
    beta.elev[1] ~ dnorm(0, .001)
    beta.elev[2] ~ dnorm(0, .001)
    beta.elev[3] ~ dnorm(0, .001)
    beta.elev[4] ~ dnorm(0, .001)
    sdy ~ dnorm(pop_sdy,0.1)T(0.05,)
    tauy <- 1 / (sdy * sdy)
    sdx ~ dnorm(pop_sdx,0.1)T(0.01,)
    taux <- 1 / (sdx * sdx)
    
    ## Likelihood
    for (i in 1:n){
    xtrue[i] ~ dnorm(0,taux)
    x[i] ~ dnorm(xtrue[i], obs_taux)
    y[i] ~ dnorm(ytrue[i], obs_tauy)
    ytrue[i] ~ dnorm(mu[i],tauy)
    mu[i] <- alpha[desc[i]] + beta.elev[desc[i]]*xtrue[i]
    }

    ## Contrasts
    elev.cont24 <- beta.elev[2] - beta.elev[4]
    elev.cont23 <- beta.elev[2] - beta.elev[3]
    elev.cont34 <- beta.elev[3] - beta.elev[4]
    elev.cont14 <- beta.elev[1] - beta.elev[4]
    elev.cont12 <- beta.elev[1] - beta.elev[2]

    }
    ", fill=T, file="jagsmodel_CIT_xyerror.txt")

# bundle data
jags_d1 <- list(x = as.numeric(elev), y = as.numeric(cit), desc = cair_snow_num, obs_taux = obs_taux,
               obs_tauy = obs_tauy, pop_sdx = pop_sdx, 
               pop_sdy = pop_sdy, n = n)
mod1 <- jags.model("jagsmodel_CIT_xyerror.txt", data=jags_d1,
                   n.chains=3, n.adapt=1000)
out1 <- coda.samples(mod1, n.iter=30000, thin=30,
                     variable.names=c("alpha", "beta.elev","sdy","elev.cont24","elev.cont23","elev.cont34",
                                      "elev.cont14","elev.cont12"))

jags_d2 <- list(x = as.numeric(elev), y = as.numeric(cit_all), desc = cair_snow_num_all, obs_taux = obs_taux_all,
                obs_tauy = obs_tauy_all, pop_sdx = pop_sdx_all, 
                pop_sdy = pop_sdy_all, n = n)
mod2 <- jags.model("jagsmodel_CIT_xyerror.txt", data=jags_d2,
                   n.chains=3, n.adapt=1000)
out2 <- coda.samples(mod2, n.iter=30000, thin=30,
                     variable.names=c("alpha", "beta.elev","sdy","elev.cont24","elev.cont23","elev.cont34",
                                      "elev.cont14","elev.cont12"))

pdf("../results/CIT_cair_elev_snow2.pdf",width=5,height=5)
par(mfrow=c(1,1))
plot_cit_jags_out(fsites_CIT_all,out2,elev_all,cit_all,cair_snow_cat_all,
                  color_pal= c(rgb(0,0,0,0),"purple",rgb(0,0,0,0),"orange"))
dev.off()

####Looks at the relationship between estimated MAT and the CIT upper bound.####
pdf("../results/CIT_thresh_MAT.pdf",width=4.5,height=4.25)
par(mfrow=c(1,1),mar=c(4,4,2,2))
plot(fsites_CIT_all$MAT_c,fsites_CIT_all$CIT_upr.97.5.,
     xlab="Eestimated Mean Annual Temperature (C)",ylab="Community Upper Threshold MAT")
abline(0,1,lty=2)
dev.off()

##Calculates the Community Climate Resistance
fsites_CIT_all$CCR <- fsites_CIT_all$CIT_upr.97.5. - fsites_CIT_all$MAT_c

##Plots CCR by covariates
pdf("../results/CCR_cair_Franklin.pdf",width=4.5,height=4)
ggplot(data=fsites_CIT_all)+
  geom_boxplot(aes(x=cair_cat_all,y=CCR))+
#  facet_wrap(facets=~COMM_NAME)+
  theme_bw()
dev.off()

####Best model predicting CCR
ccr_mod <- lm(CCR~MORA_coldair_index+MORA_elev_3m+I(MORA_elev_3m^2)+
                utmy+MORA_snow_pred_2013_3m_parkbound_est,data=fsites_CIT_all)
summary(ccr_mod)

####Is this because climate niches are narrower or closer together?



