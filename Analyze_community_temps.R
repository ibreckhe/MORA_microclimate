##Sets up workspace
library(ggplot2)
library(dplyr)
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw")

##Reads in data

fd_summary <- read.csv("franklindata_CIT_CID_summary_microMAT.csv")
fd_sens <- read.csv("../results/franklin_sites_microclimate.csv")

##Merges with CIT data.
fd_summary$loc <- paste(fd_summary$UTM_X,fd_summary$UTM_Y,sep="by")
fd_sens$loc <- paste(fd_sens$utmx,fd_sens$utmy,sep="by")

fd_summary <- merge(fd_summary,data.frame(loc=fd_sens$loc,
                                          SSMOD_SENS=fd_sens$ssmod_TAVG_sens,
                                          SSMOD_MAT=fd_sens$ssmod_MAT,
                                          GAM_tmax_sens=fd_sens$tmax_sens_est,
                                          GAM_tmin_sens=fd_sens$tmin_sens_est,
                                          GAM_tavg_sens=fd_sens$tavg_sens_est),
                    by="loc",all.x=TRUE)

##Calculates breadth.
fd_summary$CIT_br <- fd_summary$CIT_max_med - fd_summary$CIT_min_med
fd_summary$CIT_res <- fd_summary$CIT_max_med - fd_summary$MicroMAT
fd_summary$CIT_diff <- fd_summary$CIT_med - fd_summary$MicroMAT

##Calculates regional change required to push local site out of CIT_br.
fd_summary$CIT_lcr <- fd_summary$CIT_res / fd_summary$GAM_tavg_sens

##Creates an indicator value for sensitivity.
fd_summary$GAM_tavg_sens_fact <-  cut(fd_summary$GAM_tavg_sens,
                                      breaks=c(0.7,0.95,1.05,1.3),
                                      labels=c("Buffered","~1","Exposed"))
lcr_mod <- lm(CIT_lcr~GAM_tavg_sens_fact,data=fd_summary)
summary(lcr_mod)

p1 <- ggplot(fd_summary)+
        geom_point(aes(y=CIT_res,x=GAM_tavg_sens,size=num_spp),shape=21)+
        geom_smooth(aes(CIT_res,x=GAM_tavg_sens),color="black",method="lm",se=TRUE)+
        geom_vline(aes(xintercept=1),linetype="dotted")+
        geom_hline(aes(yintercept=0),linetype="dotted")+
        scale_size_continuous("Spp. \nRichness",range=c(0.1,2.5))+
        ylim(c(-1.25,2.25))+
        xlab("Local Climate Exposure")+
        ylab("Community Climate Resistance (CCR)")+
        theme_bw()
p2 <- ggplot(fd_summary)+
        geom_boxplot(aes(y=CIT_lcr,x=GAM_tavg_sens_fact),outlier.size=0)+
        ylim(c(-1,3))+
        xlab("Local Climate Exposure")+
        ylab("Regional Climate Resistance (CCR/Sensitivity)")+
        theme_bw()

#Plots it
library(gridExtra)
pdf("../results/CIT_resistance_RCR_lm.pdf",width=8,height=4)
grid.arrange(p1,p2,ncol=2,widths=c(5,3))
dev.off()

##Best model to predict CCR
model <- lm(CIT_res~GAM_tavg_sens,data=fd_summary)
summary(model)

##Splits dataset by sensitivity
fd_summary$GAM_tavg_sens_fact <- cut(fd_summary$GAM_tavg_sens,
                                     breaks=c(0.8,1,1.2))
ggplot(fd_summary)+
  geom_boxplot(aes(x=GAM_tavg_sens_fact,y=CIT_lcr))+
  theme_bw()

ggplot(fd_summary)+
  geom_point(aes(x=GAM_tavg_sens,CIT_br,size=num_spp),fill="white")+
  theme_bw()

ggplot(fd_summary)+
  geom_point(aes(x=GAM_tavg_sens,CIT_res,size=(1/CIT_br)),fill="white")+
  theme_bw()


## Generalized Additive Model Predicting CCR
library(mgcv)
fd_summary$utmx_scale <- as.numeric(scale(fd_summary$UTM_X))
fd_summary$utmy_scale <- as.numeric(scale(fd_summary$UTM_Y))

ccr_gam <- gam(CIT_res~s(Cold,Elev,k=5)+s(Canvol,Elev,k=5)+utmx_scale+utmy_scale,data=fd_summary)
plot(ccr_gam)
summary(ccr_gam)

##Computes GAM predictions.
cold_ind_med <- median(fd_summary$Cold,na.rm=TRUE)
can_vol_med <- median(fd_summary$Canvol,na.rm=TRUE)

pred_data_cold <- expand.grid(type="Cold-Air Index",
                              Cold=seq(0.3,0.8,by=0.02),
                              Elev=seq(500,2000,by=50),
                              Canvol=can_vol_med,
                              utmx_scale=0,
                              utmy_scale=0)
pred_data_can <- expand.grid(type="Canopy Volume",
                             Cold=cold_ind_med,
                             Elev=seq(500,2000,by=50),
                             Canvol=seq(0,22000,by=1000),
                             utmx_scale=0,
                             utmy_scale=0)
pred_data <- rbind(pred_data_cold,pred_data_can)

pred_data$CCR_gampred <- predict(ccr_gam,newdata=pred_data)
pred_data$CCR_gampred_se <- predict(ccr_gam,newdata=pred_data,se.fit=TRUE)$se.fit

##Estimates if CCR is significantly different from 0.
pred_data$CCR_gampred_low <- pred_data$CCR_gampred - pred_data$CCR_gampred_se * 2
pred_data$CCR_gampred_high <- pred_data$CCR_gampred + pred_data$CCR_gampred_se * 2
pred_data$CCR_gampred_sig <- as.numeric(!(pred_data$CCR_gampred_low < 0 & pred_data$CCR_gampred_high > 0))

##Plots predictions and standard errors.
p1 <- ggplot(subset(pred_data,type=="Cold-Air Index"))+
  geom_raster(aes(y=Elev/1000,x=Cold,fill=CCR_gampred))+
  geom_point(aes(y=Elev/1000,
                 x=Cold,
                 shape=factor(CCR_gampred_sig),
                 alpha=factor(CCR_gampred_sig)),
             size=0.5)+
  scale_shape_manual(values=c(1,3))+
  scale_alpha_discrete(range=c(0,0.8))+
  scale_fill_gradient2(limits=c(-0.3,1.3),midpoint=0,low="red",high="blue")+
  ylab("Elevation (km)")+
  xlab("Cold-Air Index (Unitless)")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  #facet_wrap(facets=~type)+
  guides(fill = "none",shape="none",alpha="none")+
  theme_bw()
p2 <- ggplot(subset(pred_data,type=="Canopy Volume"))+
  geom_raster(aes(y=Elev/1000,x=Canvol,fill=CCR_gampred))+
  scale_fill_gradient2("Climate \nResistance",limits=c(-0.3,1.3),midpoint=0,low="red",high="blue")+
  geom_point(aes(y=Elev/1000,
                 x=Canvol,
                 shape=factor(CCR_gampred_sig),
                 alpha=factor(CCR_gampred_sig)),
             size=0.5)+
  scale_shape_manual(values=c(1,3))+
  scale_alpha_discrete(range=c(0,0.8))+
  ylab("")+
  xlab("Canopy Volume (m^3)")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  #facet_wrap(facets=~type)+
  guides(shape="none",alpha="none")+
  theme_bw()+
  theme(axis.text.y = element_blank())

##Plots predictions of mean CCR.
library(gridExtra)
pdf("../results/microMAT_CCR_gam_2d.pdf",width=8,height=4)
grid.arrange(p1, p2, ncol=2, nrow=1,widths=c(3.5,4.5))
dev.off()

## Median slope and CCR of NVC communities.
fd_NVC_grp <- group_by(fd_summary,NVC_CLASS)
fd_NVC_sum <- summarise(fd_NVC_grp,n=n(),
                        CIT_res_NVC_med=median(CIT_res,na.rm=TRUE),
                        CIT_lcr_NVC_med=median(CIT_lcr,na.rm=TRUE),
                        Med_elev=median(Elev,na.rm=TRUE),
                        Med_mat=median(MicroMAT,na.rm=TRUE),
                        Med_cold=median(Cold,na.rm=TRUE),
                        Med_dry=median(Dry,na.rm=TRUE),
                        Med_mat_width=median(mat_width_avg,na.rm=TRUE),
                        Med_opt_sd=median(mat_opt_sd,na.rm=TRUE),
                        Med_tavg_sens=median(GAM_tavg_sens,na.rm=TRUE))
                        
fd_NVC_sum$Med_elev_cat <- cut(fd_NVC_sum$Med_elev,3)
write.csv(fd_NVC_sum,"../results/NVC_res_microMAT.csv",row.names=FALSE)

##Grabs the most common communities.
fd_NVC_sum_common <- filter(fd_NVC_sum,n>14)

ggplot(fd_NVC_sum)+
  geom_point(aes(y=CIT_lcr_NVC_med,x=Med_tavg_sens,size=n,color=Med_mat_width))+
  stat_smooth(aes(y=CIT_lcr_NVC_med,x=Med_tavg_sens,weight=n),color="black",formula=y~x,
              method="lm",se=TRUE)+
  ylim(c(-0.5,1.5))+
  theme_bw()

##Joins with plot data.                                                          
fd_NVC_citmed <- left_join(fd_summary,fd_NVC_sum,by="NVC_CLASS")
fd_NVC_common <- filter(fd_NVC_citmed,n>14)
fd_NVC_common$NVC_CLASS <- factor(fd_NVC_common$NVC_CLASS,
                                  levels=fd_NVC_sum_common$NVC_CLASS[order(fd_NVC_sum_common$CIT_res_NVC_med)])

##Creates an indicator variable for NVC associations that are rare.
common_NVC <- as.character(unique(fd_NVC_common$NVC_CLASS))
fd_NVC_sum$NVC_common <- as.numeric(fd_NVC_sum$NVC_CLASS %in% common_NVC)
fd_NVC_sum$NVC_CLASS_com <- fd_NVC_sum$NVC_CLASS
fd_NVC_sum$NVC_CLASS_com[fd_NVC_sum$NVC_common==0] <- NA
fd_NVC_sum$NVC_CLASS_com <- factor(fd_NVC_sum$NVC_CLASS_com,
                                   levels=fd_NVC_sum_common$NVC_CLASS[order(fd_NVC_sum_common$CIT_res_NVC_med)])



##Plots distribution of CCR in common NVC classes.
p1 <- ggplot(data=fd_NVC_common)+
             geom_boxplot(aes(y=CIT_res,x=NVC_CLASS,fill=NVC_CLASS),
               notch=FALSE,
               outlier.size=0)+
            #geom_point(aes(y=CIT_res,x=NVC_CLASS,color="Elev"),size=0.5,alpha=0,
            #           position=position_jitter(width=0.2))+
            geom_hline(yintercept=0,linetype="dotted")+
            scale_fill_discrete("NVC")+
            coord_flip()+
            ylim(c(-0.7,1.5))+
            ylab("CCR (C)")+
            xlab("NVC Classification")+
            theme_bw()+
            theme(legend.position = "none",
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank())

p2 <- ggplot(fd_NVC_sum)+
              geom_point(aes(y=CIT_res_NVC_med,x=Med_mat,size=n,fill=NVC_CLASS_com),shape=21)+
              geom_smooth(aes(y=CIT_res_NVC_med,x=Med_mat),color="black",method="gam",se=TRUE,
                          formula=y~s(x,k=6))+
              geom_hline(yintercept=0,linetype="dotted")+
              scale_fill_discrete("NVC Class")+
              scale_size_continuous("Prevalence (Num. Plots)",range=c(0.5,3.5))+
              xlim(limits=c(3.0,8.2))+
              ylim(limits=c(-0.5,2.0))+
              xlab("Mean Annual Temperature (C)")+
              ylab("CCR (C)")+
              theme_bw()+
              theme(legend.position = "none",
                    legend.direction = "vertical",
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank())

library(gridExtra)
pdf("../results/CCR_NVC_MAT.pdf",width=12,height=4)
grid.arrange(p1,p2,ncol=2,widths=c(5,3))
dev.off()



##Plot of niche property correlations.
p3 <- ggplot(fd_summary)+
             geom_point(aes(x=mat_width_avg,y=mat_opt_sd,color=MicroMAT))+
             theme_bw()

p4 <- ggplot(fd_summary)+
             geom_point(aes(y=mat_width_avg,x=MicroMAT))+
             theme_bw()
p5 <- ggplot(fd_summary)+
  geom_point(aes(y=mat_width_avg,x=MicroMAT))+
  theme_bw()