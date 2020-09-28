#SV Scarpino
#April 2020
#Patchiness v Crowding analysis of China and Italy

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
china_model_patch <- lm(log10(intensity) ~ log10(crowding) + log10(pop) + log10(patchiness) + log10(mobility) + log10(temp) + log10(humidity), data = china_reg_data)

italy_model_patch <- lm(log10(intensity) ~ log10(crowding) + log10(pop) + log10(patchiness) + log10(mobility) + log10(temp) + log10(humidity), data = italy_reg_data)

plot(log10(china_reg_data$crowding), log10(china_reg_data$intensity), main = "Crowding v. intensity (China)", xlab = "log crowding", ylab = "log intensity")
plot(log10(china_reg_data$patchiness), log10(china_reg_data$intensity), main = "Crowding v. patchiness (China)", xlab = "log crowding", ylab = "log intensity")

plot(log10(china_reg_data$patchiness), log10(china_reg_data$crowding), main = "Crowding v. patchiness (China)", xlab = "log patchiness", ylab = "log crowding")

plot(log10(italy_reg_data$crowding), log10(italy_reg_data$intensity), main = "Crowding v. intensity (Italy)", xlab = "log crowding", ylab = "log intensity")
plot(log10(italy_reg_data$patchiness), log10(italy_reg_data$intensity), main = "Crowding v. patchiness (Italy)", xlab = "log patchiness", ylab = "log intensity")

#extended data figure 10
plot(log10(china_reg_data$patchiness)-log10(mean(china_reg_data$patchiness, na.rm = T)), log10(china_reg_data$crowding)-log10(mean(china_reg_data$crowding, na.rm = T)), main = "Crowding v. patchiness (China)", xlab = "log patchiness", ylab = "log crowding")

if(run_glmulti == TRUE){
  glmulti.china <-
    glmulti(log10(intensity) ~ log10(pop) + log10(patchiness) + log10(crowding) + log10(mobility) + log10(temp) + log10(humidity), data = italy_reg_data,
            level = 2,               # pairwise interaction considered
            method = "h",            # Exhaustive approach
            crit = "bic",            # BIC as criteria
            confsetsize = 1,        # Keep 10 best models
            plotty = F, report = F,  # Plot or interim reports
            fitfunction = glm)
}