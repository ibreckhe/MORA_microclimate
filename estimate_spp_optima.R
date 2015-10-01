##Script to model temperature response curves for species on Mt. Rainier.
##Author: Ian Breckheimer
##Date: November 18th, 2014.

##Sets up workspace.
library(dplyr)
library(reshape2)
library(rjags)
library(ggmcmc)
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/")
source("~/code/MORA_microclimate/community_optimum_functions.R")

#### Load and Munge Data ####

##Brings in Data.
fspp <- read.csv("Franklindata_1978_sppmatrix.csv")
fspp$X <- gsub("\\s{2,}"," ",fspp$X)
fnam <- read.csv("Franklindata_sppnames.csv")
fspat <- read.csv("Franklindata_UTM.csv")
fcovar <- read.csv("franklindata_utm_covars.csv")
ecospp <- read.table("./Melanie/ecol.dat.txt",sep="\t",header=TRUE)
ecospp$plot.unique <- as.numeric(as.factor(paste(ecospp$PLOTNBR,ecospp$LAT,ecospp$LONG)))
ecospp$plot.visit <- as.factor(paste(ecospp$plot.unique,ecospp$YEAR,sep="."))
eco_coor <- unique(data.frame(ecospp$plot.visit,ecospp$LAT,ecospp$LONG,ecospp$ELEV))
colnames(eco_coor) <- c("PLOT","LAT","LONG")

##Subsets plot data to just western WA
eco_coor_wa <- eco_coor[eco_coor$LAT>45.5 & eco_coor$LONG< -119.6 & eco_coor$LAT < 50,]
plot(eco_coor_wa[,3],eco_coor_wa[,2],pch=".")
ecospp_wwa <- ecospp[ecospp$plot.visit %in% eco_coor_wa$PLOT,]
write.csv(ecospp_wwa,"./Melanie/FS_ecoshare_WWA_sppdata.csv",row.names=FALSE)
write.csv(eco_coor_wa,"./Melanie/FS_ecoshare_WWA_plotlocs.csv",row.names=FALSE)
eco_cor_wna <- read.csv("./Melanie/FS_ecoshare_WWA_plotlocs_climateWNA.csv")

##Replaces tree size codes with species codes.
fspp2 <- gsub(pattern="X+[0-9]",replacement="",colnames(fspp)[-1],perl=FALSE)
colnames(fspp) <- c("SITE",fspp2)
fnam$Code <- gsub(pattern="^[0-9]",replacement="",fnam$Code,perl=FALSE)
fnam2 <- unique(fnam)

##Melts dataset to merge columns with the same name.
fspp_mel <- melt(fspp)
fspp_mel$variable <- as.factor(as.character(fspp_mel$variable))
fspp_gr <- group_by(fspp_mel,SITE,variable)
fspp_sum <- dplyr::summarise(fspp_gr,tot_cover=sum(value))
fspp_mel2 <- melt(fspp_sum)[,-3]
colnames(fspp_mel2) <- c("SITE","SPECIES","COVER")
fspp_new <- dcast(fspp_mel2,formula=SITE~SPECIES)

##Finds MORA species represented in at least 50 plots.
numplots <- colSums(fspp_new[,-1]>0.1)
numplots_named <- merge(fnam2,data.frame(Code=names(numplots),n_plots=numplots))
numplots_ord <- numplots_named[order(numplots_named$n_plots,decreasing=TRUE),]
numplots_ord
spp_common <- numplots_ord[numplots_ord$n_plots>5,]

##Writes lists of species for TNRS taxonomy resolution (http://tnrs.iplantcollaborative.org/)
MORA_spp <- tolower(spp_common$Name)
write.csv(MORA_spp,"MORA_species.csv",row.names=FALSE)

ecospp_wwa$sp_lwr <- tolower(ecospp_wwa$SCIENTIFIC)
ECO_spp <- unique(ecospp_wwa$sp_lwr)
write.csv(ECO_spp,"ECO_species.csv",row.names=FALSE)

##Reads resolved species lists back in.
ECO_spp_tnrs <- read.csv("ECO_species_tnrs.csv")
MORA_spp_tnrs <- read.csv("MORA_species_tnrs.csv")

##Finds species that are in both lists.
joint_species <- ECO_spp_tnrs[ECO_spp_tnrs$Accepted_name %in% MORA_spp_tnrs$Accepted_name,]

##Merges resolved species lists with the data.
ecospp_wwa_tnrs <- merge(ecospp_wwa,ECO_spp_tnrs,by.x="sp_lwr",by.y="Name_submitted")
ecospp_moraspp <- ecospp_wwa_tnrs[ecospp_wwa_tnrs$Accepted_name %in% joint_species$Accepted_name,]
ecospp_reduced <- ecospp_moraspp[,c("plot.visit","Accepted_name","COVERAGE")]
ecospp_wide <- dcast(ecospp_reduced,formula=plot.visit~Accepted_name,fun.aggregate=sum)
ecospp_wide_YN <- data.frame(ecospp_wide[,1],ecospp_wide[,-1] > 0)
ecospp_numplots <- colSums(ecospp_wide[,-1]>0.1)

##Converts to presence/absence in long format.
ecospp_long <- melt(ecospp_wide,value.name="COVER")
colnames(ecospp_long) <- c("Plot.Visit","Species","Cover")
ecospp_long_coord <- merge(ecospp_long,eco_cor_wna,by.x="Plot.Visit",by.y="PLOT")
ecospp_long_coord$PresYN <- ecospp_long_coord$Cover > 0
Species <- unique(ecospp_long_coord$Species)

##Writes the joint Franklin data matrix to disk.
fnam2$sp_lwr <- tolower(fnam2$Name)
fspp_tnrs <- merge(fnam2,MORA_spp_tnrs,by.x="sp_lwr",by.y="Name_submitted")
fspp_ecospp <- fspp_tnrs[fspp_tnrs$Accepted_name %in% joint_species$Accepted_name,]
fspp_mel3 <- merge(fspp_mel2,fspp_ecospp,by.x="SPECIES",by.y="Code")
fspp_mel4 <- fspp_mel3[,c("Accepted_name","SITE","COVER")]
MORA_matrix_reduced <- dcast(fspp_mel4,formula=SITE~Accepted_name)
MORA_pres_abs_reduced <- MORA_matrix_reduced
MORA_pres_abs_reduced[,-1] <- as.numeric(MORA_pres_abs_reduced[,-1] > 0)
MORA_pres_abs_long <- melt(MORA_pres_abs_reduced)
colnames(MORA_pres_abs_long) <- c("Plot","Species","Presence")
write.csv(MORA_matrix_reduced,"MORA_matrix_reduced,csv",row.names=FALSE)
write.csv(MORA_pres_abs_reduced,"MORA_pres_abs_reduced.csv",row.names=FALSE)

##Writes the reduced Franklin data matrix to disk.
fnam2$sp_lwr <- tolower(fnam2$Name)
fspp_tnrs <- merge(fnam2,MORA_spp_tnrs,by.x="sp_lwr",by.y="Name_submitted")
fspp_mel5 <- merge(fspp_mel2,fspp_tnrs,by.x="SPECIES",by.y="Code")
fspp_mel6 <- fspp_mel5[,c("Accepted_name","SITE","COVER")]
MORA_matrix_common <- dcast(fspp_mel6,formula=SITE~Accepted_name)
MORA_pres_abs_common <- MORA_matrix_common
MORA_pres_abs_common[,-1] <- as.numeric(MORA_pres_abs_common[,-1] > 0)
MORA_pres_abs_common_long <- melt(MORA_pres_abs_common)
colnames(MORA_pres_abs_common_long) <- c("Plot","Species","Presence")
write.csv(MORA_matrix_common,"MORA_matrix_common.csv",row.names=FALSE)
write.csv(MORA_pres_abs_common,"MORA_pres_abs_common.csv",row.names=FALSE)

#### Plots species response curves fit using glm ####

# ##Models the probability of presence using glm.
# pdf("../results/fs_allspp_curves_glm.pdf",width=8,height=4)
# par(mfrow=c(1,2))
# newdata <- data.frame(MAT=seq(-10,20,by=0.01))
# plot(PresYN~MAT,data=subset(ecospp_long_coord,Species==Species[1]),
#      pch=".",col=rgb(0.25,0.25,0.25,0.8),ylim=c(0,1),main="",type="n",
#      xlab="Mean Annual Temperature (C)",ylab="Probability of Presence")
# for (i in 1:length(Species)){
#   spp_mod <- glm(PresYN~MAT+I(MAT^2),family="binomial",
#                  data=subset(ecospp_long_coord,Species==Species[i]))
#   spp_pred <- predict(spp_mod,newdata=newdata,type="response")
#   
#   ##Plots a few species vs. MAT.
#   points(newdata$MAT,spp_pred,type="l",col=Species[i])
# }
# 
# ##Models the probability of presence using glm.
# newdata <- data.frame(MAP=seq(0,5000,by=1))
# plot(PresYN~MAP,data=subset(ecospp_long_coord,Species==Species[1]),
#      pch=".",col=rgb(0.25,0.25,0.25,0.8),ylim=c(0,1),main="",type="n",
#      xlab="Mean Annual Precipitation (mm)",ylab="")
# for (i in 1:length(Species)){
#   spp_mod <- glm(PresYN~MAP+I(MAP^2),family="binomial",
#                  data=subset(ecospp_long_coord,Species==Species[i]))
#   spp_pred <- predict(spp_mod,newdata=newdata,type="response")
#   
#   ##Plots a few species vs. MAT.
#   points(newdata$MAP,spp_pred,type="l",col=Species[i])
# }
# dev.off()


# ####Bayesian model with species-level random effects.
# 
# ##Sets size of the analysis
# ntrain <- 2206 # Training plots
# ntest <- 2206 # Testing plots
# nsamples <- 5000 #MCMC samples
# 
# ##Splits data into a training and a testing set.
# all_plots <- unique(ecospp_long_coord$Plot.Visit)
# sample_plots <- sample(all_plots,ntrain+ntest,replace=FALSE)
# train_plots <- sample_plots[1:ntrain]
# test_plots <- sample_plots[(ntrain+1):(ntrain+ntest)]
# ecospp_train <- ecospp_long_coord[as.character(ecospp_long_coord$Plot.Visit)
#                                   %in% train_plots,]
# ecospp_test <- ecospp_long_coord[as.character(ecospp_long_coord$Plot.Visit)
#                                  %in% test_plots,]
# 
# ##Prepares training data for JAGS
# ecospp_train$Plot <- as.factor(paste("P.",ecospp_train$Plot.Visit,sep=""))
# ecospp_train$Plotnum <- as.numeric(ecospp_train$Plot)
# ecospp_train$X <- scale(ecospp_train$MAT)
# x_scale <- attr(ecospp_train$X,"scaled:scale")
# x_center <- attr(ecospp_train$X,"scaled:center")
# 
# ##Fits the model.
# jags.sample <- fit.jags.mixed.tree(y=ecospp_train$PresYN,x=ecospp_train$X,
#                                    species=ecospp_train$Species, thin=5,
#                                    plot=ecospp_train$Plotnum,nsamples=nsamples)
# jags.out <- jags.sample$out
# 
# ####Samples from the joint probability functions to create posterior samples of 
# 
# ####predicted CIT for the testing dataset.####
# test_reduced <- ecospp_test[,c("Plot.Visit","Species","Cover")]
# test_wide <- dcast(test_reduced,formula=Plot.Visit~Species)
# test_spp <- rowSums(test_wide[,-1] > 0)
# 
# ##Gets a vector of mean annual temperatures from the test data.
# test_MAT <- ecospp_test[,c("Plot.Visit","Species","MAT")]
# MAT_vec <- dcast(test_MAT,formula=Plot.Visit~Species)[,1:2]
# colnames(MAT_vec) <- c("Plot.Visit","MAT")
# 
# ##Makes the function matrix.
# mat1 <- make_fun_matrix(jags.out,n.samples=1000)
# 
# ##Sets the dimensions of the output.
# nsamples <-1000
# nplots <- ntest
# 
# ##Makes a matrix to hold outputs.
# CIT_samples <- matrix(NA,ncol=nsamples,nrow=nplots)
# 
# ##Loops through and calculates a CIT value for each sample.
# newx <- seq(-5,5,by=0.01)
# print("Computing CIT for sample")
# for(j in 1:nsamples){
#   print(j)
#   pred.fun <- make_combofun(mat1,iteration=j)
#   for(i in 1:nplots){
#     CIT_samples[i,j] <- sample_combo(newx,pred.fun,pres_vec=test_wide[i,-1])
#   }
# }
# rownames(CIT_samples) <- test_wide$Plot.Visit
# 
# ##Gets predicted values back on original scale.
# CIT_samples_unscaled <- CIT_samples * x_scale + x_center
# quant_fun <- function(x) quantile(x,probs = c(0.025,0.5,0.975),na.rm=FALSE)
# CIT_quants <- apply(CIT_samples_unscaled,FUN=quant_fun,MARGIN=1)
# CIT_test <- cbind(MAT_vec,t(CIT_quants),test_spp)
# CIT_fail <- CIT_test$MAT > CIT_test$'97.5%' | CIT_test$MAT < CIT_test$'2.5%'
# CIT_test$fail <- CIT_fail
# CIT_test$div_fail <- CIT_test$fail & CIT_test$test_spp > 2
# 
# ##Colors and jitter for plotting.
# cols <- c(rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1),
#           rep(rgb(0.1,0.1,0.1,1),12))
# line.cols <- c(rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1),rgb(0.1,0.1,0.1,0.1),
#                rep(rgb(0.2,0.2,0.2,0.6),12))
# xjit <- jitter(CIT_test$MAT,amount=0.1)
# 
# pdf("../results/CIT_testplots.pdf",width=5,height=5)
# plot(xjit,CIT_test$'50%',xlim=c(-2,11),ylim=c(-2,11),
#      cex=0.5,xlab="Mean Annual Temp. (C)",ylab="Inferred Mean Annual Temp. (C)",
#      col=cols[CIT_test$test_spp],type="n")
# arrows(x0=xjit,x1=xjit,y0=CIT_test$'2.5%',y1=CIT_test$'97.5%',
#        code=3,angle=90,length=0,col=line.cols[CIT_test$test_spp])
# points(xjit,CIT_test$'50%',pch=20,col=cols[CIT_test$test_spp],
#        cex=0.6)
# abline(0,1,lty=2)
# segments(x0=-1.5,x1=-1.5,y0=9.5,y1=10.5,lwd=2,col=rgb(0.1,0.1,0.1,0.1))
# points(x=-1.5,y=10,pch=20,cex=0.8,col=rgb(0.1,0.1,0.1,0.1))
# text(x=0,y=10,labels="<4 spp.",cex=0.7)
# segments(x0=-1.5,x1=-1.5,y0=7.5,y1=8.5,lwd=2,rgb(0.2,0.2,0.2,0.6))
# points(x=-1.5,y=8,pch=20,cex=0.8,col=1)
# text(x=0,y=8,labels="4+ spp.",cex=0.7)
# dev.off()
# 
# ##Predicts CIT for the reduced Franklin dataset.
# ##Sets the dimensions of the output.
# nsamples <-1000
# nplots <- dim(MORA_pres_abs_reduced)[1]
# 
# ##Makes a matrix to hold outputs.
# CIT_frank <- matrix(NA,ncol=nsamples,nrow=nplots)
# 
# ##Loops through and calculates a CIT value for each sample in the Franklin dataset.
# newx <- seq(-5,5,by=0.01)
# print("Computing CIT for sample")
# for(j in 1:nsamples){
#   print(j)
#   pred.fun <- make_combofun(mat1,iteration=j)
#   for(i in 1:nplots){
#     CIT_frank[i,j] <- sample_combo(newx,pred.fun,pres_vec=MORA_pres_abs_reduced[i,-1])
#   }
# }
# rownames(CIT_frank) <- MORA_pres_abs_reduced$SITE
# CIT_frank_unscaled <-  CIT_frank * x_scale + x_center
# write.csv(CIT_frank_unscaled,"Franklindata_CIT_samples.csv")
# 
# CIT_frank_means <- rowMeans(CIT_frank_unscaled,na.rm=TRUE)
# quantile.fun <- function(x) {quantile(x,probs=c(0.025,0.25,0.5,0.75,0.975))}
# CIT_frank_quants <- data.frame(t(apply(CIT_frank_unscaled,FUN=quantile.fun,MARGIN=1)))
# 
# ##Merges CIT values with plot metadata.
# fcovar$SITE <- paste(fcovar$SITE_NAME,fcovar$PLOT_NO,sep=" ")
# CIT_frank_quants$SITE <- rownames(CIT_frank_quants)
# fcovar_cit <- merge(fcovar,CIT_frank_quants)
# 
# ##Removes outlier solar radiation.
# fcovar_cit$MORA_srad_yearsum_9m[fcovar_cit$MORA_srad_yearsum_9m>5e08] <- NA
# fcovar_cit$ccov_trans <- logit(fcovar_cit$MORA_can_pct_focal81m)
# fcovar_cit <- fcovar_cit[complete.cases(fcovar_cit),]
# 
# cit_lm <- lm(X50.~MORA_elev_3m+MORA_coldair_index+MORA_srad_yearsum_9m+utmx*utmy,data=fcovar_cit)
# summary(cit_lm)
# AIC(cit_lm)
# 
# pdf("../results/Franklin_CIT_elev_rad.pdf",width=8,height=4)
# par(mfrow=c(1,2))
# plot(X50.~MORA_elev_3m,data=fcovar_cit,type="n",xlab="Elev. (m)",ylab="Community Inferred Temp. (C)")
# points(fcovar_cit$MORA_elev_3m,fcovar_cit$X50.,pch=20,cex=0.5)
# arrows(x0=fcovar_cit$MORA_elev_3m,x1=fcovar_cit$MORA_elev_3m,
#        y0=fcovar_cit$X25.,y1=fcovar_cit$X75.,
#        code=3,angle=90,length=0,col=rgb(0.25,0.25,0.25,0.5),lwd=0.5)
# 
# plot(X50.~MORA_srad_yearsum_9m,data=fcovar_cit,type="n",xlab="Annual Solar Radiation (Wh/M^2)",
#      ylab="",yaxt="t")
# points(fcovar_cit$MORA_srad_yearsum_9m,fcovar_cit$X50.,pch=20,cex=0.5)
# arrows(x0=fcovar_cit$MORA_srad_yearsum_9m,x1=fcovar_cit$MORA_srad_yearsum_9m,
#        y0=fcovar_cit$X25.,y1=fcovar_cit$X75.,
#        code=3,angle=90,length=0,col=rgb(0.25,0.25,0.25,0.5),lwd=0.5)
# dev.off()
# 
# ##Identifies and maps CIT outliers.
# fcovar_cit$CIT_resid <- cit_lm$residuals
# write.csv(fcovar_cit,"franklin_cit_covar_resid.csv",row.names=FALSE)

####Creates CIT based on all of the Franklin plots and species.####
MORA_spp_covar <- merge(MORA_pres_abs_common,fcovar[,c("SITE","MAT","MAT_c","MORA_elev_3m","MORA_dry_index_81m_precip")])
MORA_covar <- MORA_spp_covar[,c("SITE","MAT","MAT_c","MORA_elev_3m","MORA_dry_index_81m_precip")]
MORA_spp <- MORA_spp_covar[,-c(which(colnames(MORA_spp_covar) %in% c("MAT","MAT_c","MORA_elev_3m","MORA_dry_index_81m_precip")))]

##Separates training and testing plots.
set.seed(42)
rand <- runif(n=dim(MORA_covar)[1],0,1)
MORA_train_spp <- MORA_spp[rand > 0.9,]
MORA_train_spp_long <- melt(MORA_train_spp)
f_train <- merge(MORA_train_spp_long,MORA_covar,all.x = TRUE)
colnames(f_train) <- c("SITE","SPECIES","PRES_ABS","WNA_MAT","PRSM_MAT","elev","dry")

MORA_test_spp <- MORA_spp[rand <= 0.9,]
MORA_test_spp_long <- melt(MORA_test_spp)
f_test <- merge(MORA_test_spp_long,MORA_covar,all.x = TRUE)
colnames(f_test) <- c("SITE","SPECIES","PRES_ABS","WNA_MAT","PRSM_MAT","elev","dry")

MORA_spp_long <- melt(MORA_spp)
f_all <- merge(MORA_spp_long,MORA_covar,all.x=TRUE)
colnames(f_all) <- c("SITE","SPECIES","PRES_ABS","WNA_MAT","PRSM_MAT","elev","dry")

##Ensures training data are complete.
f_train <- f_train[complete.cases(f_train),]
f_all <- f_all[complete.cases(f_all),]

##Scales training data for JAGS
f_dry_train <- scale(f_train$dry)
f_dry_train_scale <- attr(f_dry,"scaled:scale")
f_dry_train_center <- attr(f_dry,"scaled:center")
f_mat_train <- scale(f_train$PRSM_MAT)
f_mat_train_scale <- attr(f_mat,"scaled:scale")
f_mat_train_center <- attr(f_mat,"scaled:center")
f_sppnum_train <- as.numeric(as.factor(f_train$SPECIES))
f_plotnum_train <- as.numeric(as.factor(f_train$SITE))
f_y_train <- f_train$PRES_ABS

##Scales full data for JAGS
f_dry <- scale(f_all$dry)
f_dry_scale <- attr(f_dry,"scaled:scale")
f_dry_center <- attr(f_dry,"scaled:center")
f_mat <- scale(f_all$PRSM_MAT)
f_mat_scale <- attr(f_mat,"scaled:scale")
f_mat_center <- attr(f_mat,"scaled:center")
f_sppnum <- as.numeric(as.factor(f_all$SPECIES))
f_plotnum <- as.numeric(as.factor(f_all$SITE))
f_y <- f_all$PRES_ABS

# ## Fits the univariate model
# nsamples <- 5000
# f.jags.sample <- fit.jags.mixed.tree(y=f_y,x=f_mat,
#                                    species=f_sppnum, thin=5,
#                                    plot=f_plotnum,nsamples=nsamples)
# save(f.jags.sample,file="CIT_franklindata_jagsoutput.R",compress=TRUE)
# f.jags.out <- f.jags.sample$out

##Fits a model with two interacting variables
nsamples <- 5000
f.jags.2var <- fit.jags.mixed.2var(y=f_y,x1=f_mat,x2=f_dry,
                                     species=f_sppnum, thin=5,
                                     plot=f_plotnum,nsamples=nsamples)
save(f.jags.2var,file="franklindata_jagsoutput_2var.R",compress=TRUE)
require(mcmcplots)
mcmcplot(f.jags.2var$out,regex=c("\\.mu","\\.sigma","int","width","opt.g.x2"))
load("franklindata_jagsoutput_2var.R")

## Makes the function matrix for the interaction model
f.jags.2var.out <- f.jags.2var$out
f.mat.2var <- make_fun_matrix_2var(f.jags.2var.out,n.samples=1000)
#save(f.mat.2var,file="franklindata_funmatrix_2var.R",compress=TRUE)

## Makes a function matrix for the median parameter estimates.
f.mat.2var.med <- make_med_fun_matrix_2var(f.jags.2var.out)

## Plots response curves for species.
require(fields)
pdf("../results/Temp_Dry_Spp_Curves.pdf",width=10,height=10)
breaks <- seq(0,1,length.out=1001)
cols <- rev(terrain.colors(length(breaks)-1))
set.panel()
par(mar=c(2,3,3,0),oma=c(3,3,0,10))
set.panel(4,4)
for(i in 1:16){
  x1_seq <- seq(-4,4,length.out=100)
  x2_seq <- seq(-4,4,length.out=100)
  x1_seq_un <- x1_seq * f_mat_scale + f_mat_center
  x2_seq_un <- x2_seq * f_dry_scale + f_dry_center
  x_grid <- expand.grid(x1=x1_seq,x2=x2_seq)
  x_grid$y <- f.mat.2var.med[1,i][[1]](x_grid$x1,x_grid$x2)
  image(x=x1_seq_un,y=x2_seq_un,z=matrix(x_grid$y,ncol=100),xlab="",ylab="",
        main=levels(as.factor(f_all$SPECIES))[i],breaks=breaks,col=cols)
}
set.panel()
mtext(text="Mean Annual Temperature",side = 1,outer = TRUE,line=1,adj=0.6)
mtext(text="Topographic Dryness Index",side = 2,outer = TRUE,line=1)
par(oma=c(0,10,0,1))
image.plot(x=x1_seq_un,y=x2_seq_un,z=matrix(runif(length(x_grid$y),0,1),ncol=100),xlab="",ylab="",
            main=levels(as.factor(f_train$SPECIES))[1],breaks=breaks,col=cols,legend.only=TRUE,
           legend.shrink=0.6,legend.args=list(text="Prob.",adj=-1))
dev.off()

##Sets the dimensions of the output.
nsamples <-1000
nplots <- dim(MORA_spp)[1]

##Makes an array to hold outputs.
f_CIT_samples <- matrix(NA,ncol=nsamples,nrow=nplots)
f_CID_samples <- matrix(NA,ncol=nsamples,nrow=nplots)

##Loops through and calculates a CIT value for each sample.
print("Computing CIT and CID for sample")
for(j in 1:nsamples){
  print(j)
  pred.fun <- make_combofun_2var(f.mat.2var,iteration=j)
  for(i in 1:nplots){
    newx1 <- runif(1000,-5,5)
    newx2 <- runif(1000,-5,5)
    accepted <- sample_combo_2var(newx1,newx2,pred.fun,pres_vec=MORA_spp[i,-1])
    f_CIT_samples[i,j] <- accepted[1]
    f_CID_samples[i,j] <- accepted[2]
  }
}
rownames(f_CIT_samples) <- MORA_spp$SITE
rownames(f_CID_samples) <- MORA_spp$SITE


##Gets predicted values back on the original scale.
f_CIT_samples_unscaled <- f_CIT_samples * f_mat_scale + f_mat_center
f_CID_samples_unscaled <- f_CID_samples * f_dry_scale + f_dry_center

##Writes samples to disk.
write.csv(f_CIT_samples_unscaled,"Franklindata_CIT_samples_all.csv")
write.csv(f_CID_samples_unscaled,"Franklindata_CID_samples_all.csv")

##Reads them back in
f_CIT_samples_unscaled <- read.csv("Franklindata_CIT_samples_all.csv")
rownames(f_CIT_samples_unscaled) <- f_CIT_samples_unscaled$X
f_CIT_samples_unscaled <- as.matrix(f_CIT_samples_unscaled[,-1])
f_CID_samples_unscaled <- read.csv("Franklindata_CID_samples_all.csv")
rownames(f_CID_samples_unscaled) <- f_CID_samples_unscaled$X
f_CID_samples_unscaled <- as.matrix(f_CID_samples_unscaled[,-1])

#Computes CIT and CID covariance.
require(MASS)
require(cluster)
slopes <- data.frame(SITE=MORA_spp$SITE,CIT_CID_slope=rep(NA,nplots),CIT_CID_covar=rep(NA,nplots))
for (i in 1:nplots){
  mcmc <- cbind(f_CIT_samples[i,],f_CID_samples[i,])
  samples <- cov.mve(mcmc,quantile.used=nrow(mcmc)*0.75)
  ellipse <- ellipsoidhull(mcmc[samples$best,])
  bound <- predict(ellipse)
  axis <- ellipse_majoraxis(ellipse)
  slopes$CIT_CID_slope[i] <- axis[2]
  slopes$CIT_CID_covar[i] <- samples$cov[1,2]
}

quant_fun <- function(x) quantile(x,probs = c(0.025,0.5,0.975),na.rm=FALSE)
f_CIT_quants <- apply(f_CIT_samples_unscaled,FUN=quant_fun,MARGIN=1)
f_CID_quants <- apply(f_CID_samples_unscaled,FUN=quant_fun,MARGIN=1)
f_all_covar <- merge(MORA_spp,MORA_covar,all.x = TRUE)
f_slopes <- merge(MORA_spp,slopes,all.x = TRUE)

##Computes number of species in each plot.
n_spp <- rowSums(MORA_spp[,-1] > 0)
f_CIT_all <- data.frame(cbind(f_all_covar$MAT,f_all_covar$MAT_c,f_all_covar$MORA_elev_3m,
                              f_all_covar$MORA_dry_index_81m_precip,
                              t(f_CIT_quants),t(f_CID_quants),f_slopes$CIT_CID_slope,
                              f_slopes$CIT_CID_covar,n_spp))
colnames(f_CIT_all) <- c("WNA_MAT","PRISM_MAT","Elev","Dry","CIT_lwr","CIT_med","CIT_upr",
                         "CID_lwr","CID_med","CID_upr","CIT_CID_slope","CIT_CID_covar","nspp")
f_CIT_all$DryFact <- cut(f_CIT_all$Dry,breaks=quantile(f_CIT_all$Dry,probs=c(0.01,0.25,0.5,0.75),na.rm=TRUE))
f_CIT_all$ColdFact <- cut(f_CIT_all$PRISM_MAT,breaks=quantile(f_CIT_all$PRISM_MAT,probs=c(0.01,0.25,0.5,0.75),na.rm=TRUE))

##Writes summary stats to disk.
write.csv(f_CIT_all,"franklindata_CIT_CID_summary.csv",row.names=FALSE)

##Plots distribution of covariances
pdf("../results/CIT_CID_covar_density.pdf",width=6,height=6)
ggplot(f_CIT_all)+
  geom_point(aes(x=Elev,y=CIT_CID_covar,col=ColdFact))+
  stat_smooth(aes(x=Elev,y=CIT_CID_covar),method="loess")+
  scale_x_continuous(limits=c(400,2100))+
  scale_y_continuous(limits=c(-0.1,0.05))+
  theme_bw()
dev.off()

##Plots CID and CIT
require(MASS)
require(cluster)
pdf("../results/CIT_CID_covar_lines.pdf",width=6,height=6)
par(mfrow=c(1,1),mar=c(4,4,0,0),oma=c(1,1,1,1))
cols <- rep(1,nplots)
plot(f_CIT_samples[1,],f_CID_samples[1,],xlim=c(-3,3),ylim=c(-3,3),type="n",
     xlab="Scaled Community MAT",ylab="Scaled Community Dryness")
for (i in 1:nplots){
  mcmc <- cbind(f_CIT_samples[i,],f_CID_samples[i,])
  samples <- cov.mve(mcmc,quantile.used=nrow(mcmc)*0.75)
  ellipse <- ellipsoidhull(mcmc[samples$best,])
  bound <- predict(ellipse)
  axis <- ellipse_majoraxis(ellipse)
  xbnds <- quantile(mcmc[,1],probs=c(0.25,0.75))
  yaxis <- axis[2] * xbnds + axis[1]
  #points(x=f_CIT_samples[i,],y=f_CID_samples[i,],
  #       col=adjustcolor(cols[i],alpha.f = 0.05),pch=20,cex=0.4)
  #polygon(bound,col=adjustcolor(cols[i],alpha.f = 0.3),border=NA)
  abline(axis[1],axis[2],lwd=0.3,col=adjustcolor(cols[i],alpha.f = 0.5))
  points(x=ellipse$loc[1],y=ellipse$loc[2],cex=1,pch=20,
         col=adjustcolor(cols[i],alpha.f = 0.8))
  #segments(x0=x0,x1=x1,y0=y0,y1=y1,lwd=2,col=adjustcolor(i,alpha.f = 0.3))
  #points(f_CIT_samples_unscaled[i,],f_CID_samples_unscaled[i,],col=i,pch=20,cex=0.2)
}
dev.off()

##Plots CIT and MAT for all plots.

##Colors and jitter for plotting.
cols <- rgb(0.1,0.1,0.1,1)
line.cols <- rgb(0.2,0.2,0.2,0.4)

pdf("../results/CIT_CID_MORA_allplots.pdf",width=10,height=5)
par(mfrow=c(1,2),mar=c(4,4,2,2))
xjit <- jitter(f_CIT_all$PRISM_MAT,amount=0.001)
plot(xjit,f_CIT_all$CIT_med,xlim=c(2,9),ylim=c(2,9),
     cex=0.5,xlab="Mean Annual Temp. (C)",ylab="Community Inferred Mean Annual Temp. (C)",
     col=cols,type="n")
arrows(x0=xjit,x1=xjit,y0=f_CIT_all$CIT_lwr,y1=f_CIT_all$CIT_upr,
       code=3,angle=90,length=0,col=line.cols)
points(xjit,f_CIT_all$CIT_med,pch=20,col=cols,
       cex=0.6)
abline(0,1,lty=2)
xjit <- jitter(f_CIT_all$Dry,amount=0.0001)
plot(xjit,f_CIT_all$CID_med,xlim=c(0.04,0.16),ylim=c(0.04,0.16),
     cex=0.5,xlab="Dryness Index",ylab="Community Inferred Dryness Index",
     col=cols,type="n")
arrows(x0=xjit,x1=xjit,y0=f_CIT_all$CID_lwr,y1=f_CIT_all$CID_upr,
       code=3,angle=90,length=0,col=line.cols)
points(xjit,f_CIT_all$CID_med,pch=20,col=cols,
       cex=0.6)
abline(0,1,lty=2)
dev.off()

pdf("../results/CIT_MORA_elev.pdf",width=5,height=5)
xjit <- jitter(f_CIT_all$Elev,amount=0.001)
par(mfrow=c(1,1))
plot(xjit,f_CIT_all$CIT_med,xlim=c(500,2000),ylim=c(2,9.5),
     cex=0.5,xlab="Elevation (m)",ylab="Community Inferred Mean Annual Temp. (C)",
     col=cols,type="n")
arrows(x0=xjit,x1=xjit,y0=f_CIT_all$CIT_lwr,y1=f_CIT_all$CIT_upr,
       code=3,angle=90,length=0,col=line.cols)
points(xjit,f_CIT_all$CIT_med,pch=20,col=cols,
       cex=0.6)
abline(lm(PRISM_MAT~Elev,data=f_CIT_all),lty=2)
dev.off()

pdf("../results/CID_MORA_elev.pdf",width=5,height=5)
xjit <- jitter(f_CIT_all$Elev,amount=0.001)
par(mfrow=c(1,1))
plot(xjit,f_CIT_all$CID_med,xlim=c(500,2000),ylim=c(0.04,0.16),
     cex=0.5,xlab="Elevation (m)",ylab="Community Inferred Dryness Index",
     col=cols,type="n")
arrows(x0=xjit,x1=xjit,y0=f_CIT_all$CID_lwr,y1=f_CIT_all$CID_upr,
       code=3,angle=90,length=0,col=line.cols)
points(xjit,f_CIT_all$CID_med,pch=20,col=cols,
       cex=0.6)
#abline(lm(PRISM_MAT~Elev,data=f_CIT_all),lty=2)
dev.off()
