#SV Scarpino
#April 2020
#Crowding intensity and attack rate analysis of China and Italy

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
run_glmulti <- FALSE
do_plot <- FALSE


###########
#Acc Funcs#
###########


#########################
#Build combined datasets#
#########################
source("italy_crowding_cross_validation.R")
source("china_crowding_cross_validation.R")

############
#Fit models#
############
quartz()
layout(matrix(1:4, nrow = 2, byrow = TRUE))
ints <- c()
vars <- c()
for(l in LETTERS[1:4]){
  i.all <- which(group == l & log10(results[[1]]$total_pops) < 6.5 & log10(results[[1]]$total_pops) > 6.4)
  i <- i.all[1]
  objectid.i <- china$OBJECTID[which(china$ADM2_PCODE %in% prefects[i])]
  use.int.i <- which(china_cases$OBJECTID == objectid.i)
  ord.int.i <- order(china_cases$date[use.int.i])
  new.cases.i <- china_cases$Cases[use.int.i][ord.int.i]
  p.i <- new.cases.i/sum(new.cases.i, na.rm = TRUE) #calculate fraction of total cases for each week
  ints <- c(ints, entropy(p.i)^-1)
  vars <- c(vars, var(diff(new.cases.i), na.rm = TRUE))
  plot(p.i, type = "b", main = paste0(l, ": Intensity = ", round(entropy(p.i)^-1,2), ", Variance = ", round(var(diff(new.cases.i), na.rm = TRUE),2)), ylim = c(0, 1), xlab = "Days since first case", col = cols_list[[l]], pch = 16, ylab = "Proportion")
  points(cumsum(p.i), type = "l", col = cols_list[[l]])
}
