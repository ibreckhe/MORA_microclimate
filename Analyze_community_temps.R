##Sets up workspace
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw")

##Reads in data
fd_summary <- read.csv("franklindata_CIT_CID_summary.csv")

##Calculates breadth.
elev <- fd_summary$Elev
dry <- fd_summary$DryYN
br <- fd_summary$CIT_upr - fd_summary$CIT_lwr

##Eliminates Outliers
elev_out <- elev[br<2.5 & elev<1600]
br_out <- br[br<2.5 & elev<1600]
dry_out <- dry[br<2.5 & elev<1600]

#Plots it
plot(elev_out,br_out, xlab= "Elevation (m)",ylab="Community Temperature Range (C)")

model <- lm(br_out~elev_out+I(elev_out^2))
model2 <- lm(br_out~elev_out)
model3 <- lm(br_out~dry_out)

xnew <- 500:1600
ynew <- coef(model)[1] + coef(model)[2] * xnew + coef(model)[3] * xnew^2
points(xnew,ynew,type="l",col="red",lwd=2)


