#SV Scarpino
#April 2020
#Crowding analysis of China and Italy

###########
#libraries#
###########
library(ggplot2)
library(dplyr)
library(tiff)
library(rgeos)
library(maptools)
library(raster)
library(proj4)
library(rgdal)
library(sp)
library(epiphy)
library(entropy)
library(RVAideMemoire)
library(glmulti)

#########
#Globals#
#########
time_stamp <- as.numeric(Sys.time())
save_new <- FALSE

###########
#Acc Funcs#
###########
#Crowding was calculated at various spatial resolutions as defined by arc minutes (for reasons outlined in earlier discussion): 30 Second ( ~ 1km), 2.5 minutes ( ~ 5km), 15 minutes ( ~ 30 km), 30 minutes ( ~ 55km), 60 minutes (~110km)

######
#Data#
######
global <- readRDS("../Data/380_cities_covariates_with_geometry_5_6_20.csv")

####################
#Run china analysis#
####################
source("china_crowding_cross_validation_4k.R")
source("italy_crowding_cross_validation_4k.R")

combined_reg_data <- rbind(italy_reg_data, china_reg_data)

################
#Building Model#
################
combined_model <- lm(log10(intensity) ~ log10(crowding) + log10(mobility) + log10(humidity), data = combined_reg_data)

############
#Predicting#
############
globe.pred <- data.frame(global$Population, global$crowding2point5m, 1/global$within_april_mean_over_december, global$tempF, global$spechum)
colnames(globe.pred) <- c("pop", "crowding", "mobility", "temp", "humidity")

pred.global <-predict(object = china_model, newdata = globe.pred, type = "response")

dat.out <- data.frame(pred.global, global[,-ncol(global)])
colnames(dat.out)[1] <- "intensity_pred_log10"
dat.out$intensity_pred_log10[which(dat.out$intensity_pred_log10 > 0)] <- NA #censoring noisy outbreaks

if(save_new == TRUE){
  file_name <- paste0(time_stamp, "prediction_380_cities_covariates_with_geometry_7_15_20.csv")
  write.csv(dat.out, file = file_name, row.names = FALSE)
}