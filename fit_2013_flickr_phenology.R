##Script to fit phenology curves to the 2013 flickr data.

####Setup####

##Sets up workspace
library(ggplot2)
library(ggmap)
library(lme4)
library(plyr)
setwd("~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Phenology/Flickr&iNaturalist/Flickr/")

##Sources my functions.
source("~/code/MORA_microclimate/phenology_functions.R")

##Brings in 2013 flickr data
fl <- read.csv("MORA_flickr_classified_all_2013_covars.csv")

####Format and filter data####

##Converts date taken to julian day.
fl$datetaken <- strptime(fl$datetaken,format="%m/%d/%y %H:%M")
fl$datetaken <- as.POSIXct(fl$datetaken) # convert to POSIXct to use in data frames / ddply
fl$datetaken_DOY <- as.numeric(strftime(fl$datetaken,format="%j"))

##Converts flowering NA's to NF.
sp <- as.character(fl$SPECIES)
sp[is.na(sp)] <- "NF"
sp <- as.factor(sp)
fl$SPECIES <- sp

##Creates a binary column for the presence of focal species.
fl$focalsp <- fl$SPECIES != "NF"

##Creates days since snow variable.
fl$dss <- fl$datetaken_DOY - fl$snow_diss_2013
hist(fl$dss)

##Filters data to bring in only data in the park bounds.
fl <- fl[!is.na(fl$can_pct),]

##Filters data to get only records with high spatial accuracy
fl <- fl[fl$accuracy==16,]

##Filters data to exclude photos not from the park.
fl <- fl[-c(which(fl$SPECIES=="NP")),]

##Filters data to exclude flowers on glacier.
fl <- fl[-c(which(fl$SPECIES!="NF" & fl$snow_melt==2)),]

##Filters data to exclude flowers on snow.
#fl <- fl[-c(which(fl$focalsp==TRUE & fl$dss < -1)),]

##Removes photos from user 88227046@N00
fl <- fl[-c(which(fl$owner=="88227046@N00")),]

##Creates a unique factor for photos taken on a particular 
##date by a particular person.
fl$ownerdate <- factor(paste(fl$owner,fl$datetaken_DOY,sep="_"))

##Puts all of the points on a map to check the spatial distribution.
##Quick map to check data.
bbox <- c(min(fl$long) - 0.05,
          min(fl$lat) - 0.05,
          max(fl$long) + 0.05,
          max(fl$lat) + 0.05)

gmap <- get_map(bbox)

map <- ggmap(gmap)+
  geom_point(aes(x=long,y=lat),
             data=fl,
             alpha=0.3,
             position=position_dodge(width = 0.01,height=0.01))+
  geom_point(aes(x=long,y=lat,color=SPECIES),
             data=subset(fl,PHEN_PHASE=="FLOWERING"),
             position=position_dodge(width = 0.01,height=0.01))
map

####Bins data for modeling####
# Creates the bin breaks for dss.
dss_breaks <- seq(-50,150,by=1)

# Loops through each species and counts the number of successes and trials
# in each bin

# Focal species
spp <- c("Erythronium montanum",
         "Erigeron peregrinus",
         "Valeriana sitchensis",
         "Polygonum bistortoides",
         "Anemone occidentalis",
         "Castilleja parviflora",
         "Lupinus arcticus") 

summed_all <- data.frame()
for (i in 1:length(spp)){
  classed <- class_obs(fl,"SPECIES",spp[i],"PHEN_PHASE","FLOWERING")
  binned <- bin_obs(classed,column="dss",breaks=dss_breaks,
                    new_colname="dss_bin")
  summed <- sum_obs(binned,"dss_bin",success_colname=spp[i],
                    trial_colname="All_Photos")
  summed <- summed[-c(dim(summed)[1]),] #Drops observations outside of the bin range.
  summed$Species <- spp[i]
  summed$dss_num <- as.numeric(as.character(summed$dss_bin))
  colnames(summed) <- c("dss_bin","Flower_Phot","All_Phot","Species")
  summed_all <- rbind(summed,summed_all)
}
colnames(summed_all) <- c("dss_bin","Flower_Phot","All_Phot","Species","dss_num")

# Creates the response variable: proportion of photos with a flower.
summed_all$prop_flower <- summed_all$Flower_Phot / summed_all$All_Phot

#### Models using glm.####
null_form <- formula(prop_flower~1)
quad_form <- formula(prop_flower~dss_num+I(dss_num^2))
log_form <- formula(prop_flower~log(dss_num)+I(log(dss_num)^2))

# Models with just a mean
null_mods <- list()
for (i in 1:length(spp)){
  null_mod <- glm(null_form,
                  weights=All_Phot,
                  data=subset(summed_all, Species==spp[i]),
                  family=binomial(link="logit"))
  null_mods[[i]] <- null_mod
}
names(null_mods) <- spp

# Models with a logit-quadratic (gaussian) trend
quad_mods <- list()
for (i in 1:length(spp)){
  quad_mod <- glm(quad_form,
                  weights=All_Phot,
                  data=subset(summed_all, Species==spp[i]),
                  family=binomial(link="logit"))
  quad_mods[[i]] <- quad_mod
}
names(quad_mods) <- spp

# Models with a log-gaussian trend
log_mods <- list()
for (i in 1:length(spp)){
  log_mod <- glm(log_form,
                 weights=All_Phot,
                 data=subset(summed_all, Species==spp[i]),
                 family=binomial(link="logit"))
  log_mods[[i]] <- log_mod
}
names(log_mods) <- spp

# Computes AIC values for all species for the three types of models.
aic_table <- data.frame(null_aic=aic_list(null_mods),
                        quad_aic=aic_list(quad_mods),
                        log_aic=aic_list(log_mods))
aic_table$delta <- aic_table$quad_aic - aic_table$log_aic
aic_table

# Finds the species for which the log model does better better.
log_table <- subset(aic_table,delta > -1)
gaus_table <- subset(aic_table,delta <= -1)

####Uses the predict function to generate predicted values####
dss_days <- seq(-50,150,by=0.1)

preds <- data.frame()
for (i in 1:length(spp)){
  pred_gauss <- predict(quad_mods[[i]],newdata=data.frame(dss_num=dss_days),type="response")
  pred_log <- predict(log_mods[[i]],newdata=data.frame(dss_num=dss_days),type="response")
  pred_gauss <- data.frame(type="gaussian",dss=dss_days,pred=pred_gauss)
  pred_log <- data.frame(type="log-gaussian",dss=dss_days,pred=pred_log)
  preddf <- data.frame(Species=spp[i],rbind(pred_log,pred_gauss))
  preds <- rbind(preds,preddf)
}


####Plots the predictions and binned data.
g <- ggplot(data=preds,groups=Species)+
  facet_wrap(~Species,ncol=4)+
  geom_path(aes(x=dss,y=pred,color=Species,linetype=type),data=preds)+
  geom_point(aes(x=dss_num,y=prop_flower,color=Species,size=All_Phot),
             data=summed_all,shape=21)+
  scale_linetype_discrete(name="Fit Type",guide='legend')+
  scale_size_continuous(name="Total # \nPhotos",guide='legend')+
  scale_color_discrete(name="Species",guide=FALSE)+
  coord_cartesian(xlim=c(-30,150),ylim=c(-0.03,0.2))+
  xlab("Days Since Snow")+
  ylab("Proportion of All Photos with Flowers")+
  theme_bw()+
  theme(legend.position=c(0.85,0.20),
        legend.key=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background=element_blank())
g

####Computes optimum and tolerance from model predictions####

# Gets rid of NaN predictions for observations less than zero.
preds$pred[is.nan(preds$pred)] <- NA

# Finds the maxima
optimums <- ddply(preds,~Species*type,.fun=summarise,
                  max=dss[which.max(pred)])

# Finds the upper and lower 68.5% cumulative probability intervals.
intervals <- ddply(preds,~Species*type,.fun=obs_intervals,
                   threshold=0.158,bin_width=1)

# Assembles all the relevant parameters.
intervals$opt <- optimums$max
intervals$std.dev <- (intervals$upr_bound - intervals$lwr_bound) / 2
intervals$AIC <- NA
intervals$AIC[intervals$type=="gaussian"] <- aic_list(quad_mods)
intervals$AIC[intervals$type=="log-gaussian"] <- aic_list(log_mods)
intervals

# Adds the intervals to the previous plot.
f <- g + geom_segment(data=subset(intervals,type == "log-gaussian"),
                      aes(x=lwr_bound,xend=upr_bound,y=-0.01,yend=-0.01,color=Species),
                      linetype=1,linewidth=3)+
  geom_segment(data=subset(intervals,type == "gaussian"),
               aes(x=lwr_bound,xend=upr_bound,y=-0.02,yend=-0.02,color=Species),
               linetype=2,linewidth=3)
f

# Adds points at the optima.
# Adds the intervals to the previous plot.
h <- f + geom_point(data=subset(optimums,type == "log-gaussian"),
               aes(x=max,y=-0.01,color=Species),size=3)+
  geom_point(data=subset(optimums,type == "gaussian"),
             aes(x=max,y=-0.02,color=Species),size=3,shape=21,fill="white")
h

#### Plots curves of the best model by species.

# Assembles the best models into a list.
best_mods <- list()
for (i in 1:length(rownames(log_table))){
  best_mods[[rownames(log_table)[i]]] <- log_mods[[rownames(log_table)[i]]]
}
for (i in 1:length(rownames(gaus_table))){
  best_mods[[rownames(gaus_table)[i]]] <- quad_mods[[rownames(gaus_table)[i]]]
}
names(best_mods)

# Gets predictions from the best model in a data frame.
best_preds <- data.frame()
for (i in 1:length(spp)){
  pred_best <- predict(best_mods[[i]],newdata=data.frame(dss_num=dss_days),type="response")
  preddf <- data.frame(Species=names(best_mods)[i],dss=dss_days,pred_best)
  best_preds <- rbind(best_preds,preddf)
}

####Plots the predictions for the best models.

# colors consistent with JHRL analysis
plotcol <- c("salmon","yellowgreen","red","lightpink",
             "lightblue","maroon","plum")

k <- ggplot(data=best_preds,groups=Species)+
  geom_path(aes(x=dss,y=pred_best,color=Species),data=best_preds)+
  coord_cartesian(xlim=c(-5,100))+
  scale_color_manual(values=plotcol)+
  xlab("Days Since Snow")+
  ylab("Proportion of All Photos with Flowers")+
  theme_bw()+
  theme(legend.position='right',
        legend.key=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background=element_blank())
k

####Computes optimum and tolerance from the best model predictions####

# Gets rid of NaN predictions for observations less than zero.
best_preds$pred_best[is.nan(best_preds$pred_best)] <- NA

# Finds dss at the maxima
best_max <- ddply(best_preds,~Species,.fun=summarise,
                  max=max(pred_best,na.rm=TRUE))

# Finds dss at the maxima
best_optimums <- ddply(best_preds,~Species,.fun=summarise,
                       optim=dss[which.max(pred_best)])

# Finds the upper and lower 68.5% cumulative probability intervals.
best_intervals <- ddply(best_preds,~Species,.fun=obs_intervals,
                        threshold=0.158,bin_width=1)
best_intervals$rangef <- (best_intervals$upr_bound - best_intervals$lwr_bound)/2

####Assembles predictions####
best_stats <- cbind(best_optimums,best_intervals$rangef,best_max$max)
colnames(best_stats) <- c("species","peakf","rangef","maxf")
format(best_stats,digits=4)

####Boostraps the best models to assess uncertainty in the fit.####
start <- Sys.time()
boots <- data.frame()
boot.opts <- data.frame()
boot.tols <- data.frame()
for (i in 1:length(spp)){
  flush.console()
  print(paste("Now bootstrapping the model for ",
        best_mods[[i]]$data$Species[1],
        ". Species ", i," of ", length(spp),sep=""))
  boot <- boot_preds2(best_mods[[i]],
                     data=fl,
                     dss_breaks=dss_breaks,
                     pred.data=dss_days,
                     clustervar="ownerdate",
                     n_replicates=5000)
  boots <- rbind(boots,boot[[1]])
  boot.opts <- rbind(boot.opts,data.frame(species=best_mods[[i]]$data$Species[1],
                                          opt_2.5=boot[[2]][1],
                                          opt_97.5=boot[[2]][2]))
  boot.tols <- rbind(boot.tols,data.frame(species=best_mods[[i]]$data$Species[1],
                                          tol_2.5=boot[[3]][1],
                                          tol_97.5=boot[[3]][2]))
}
elapsed <- Sys.time() - start
elapsed

##Appends the bootstrap intervals to the summary table.
best.stats.boot <- merge(x=best_stats,y=boot.opts,by="species")
best.stats.boot2 <- merge(x=best.stats.boot,y=boot.tols,by="species")
colnames(best.stats.boot2) <- c("Species","peak","range","max","peak_2.5","peak_97.5","range_2.5","range_97.5")
format(best.stats.boot2,digits=3)

##graphs the bootstrap intervals.
best.preds.boot <- cbind(best_preds,boots$lwr_bnd,boots$upr_bnd)
colnames(best.preds.boot) <- c('Species','dss','pred_best','lwr_bnd','upr_bnd')
best.preds.boot$range <- best.preds.boot$upr_bnd - best.preds.boot$lwr_bnd

####Plots the predictions and binned data.
l <- ggplot(data=best.preds.boot,groups=Species)+
  facet_wrap(~Species,ncol=4)+
  geom_ribbon(aes(x=dss,ymin=lwr_bnd,ymax=upr_bnd,fill=Species),data=best.preds.boot,alpha=0.3)+
  geom_path(aes(x=dss,y=pred_best,color=Species),data=best.preds.boot)+
  geom_point(aes(x=dss_num,y=prop_flower,color=Species,size=All_Phot),
             data=summed_all,shape=21)+
  geom_segment(data=best_intervals,
               aes(x=lwr_bound,xend=upr_bound,y=-0.018,yend=-0.018,color=Species),
               linetype=1,linewidth=3)+
  geom_point(data=best_optimums,
           aes(x=max,y=-0.018,color=Species),size=3)+
  scale_size_continuous(name="Total # \nPhotos",guide='legend')+
  scale_color_discrete(name="Species",guide=FALSE)+
  scale_fill_discrete(name="Species",guide=FALSE)+
  coord_cartesian(xlim=c(-30,150),ylim=c(-0.03,0.2))+
  xlab("Days Since Snow")+
  ylab("Proportion of All Photos with Flowers")+
  theme_bw()+
  theme(legend.position=c(0.85,0.20),
        legend.key=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background=element_blank())
l

setwd("../analysis")
svg("flickr_2013_bootstrap.svg",width=8,height=5)
l
dev.off()

