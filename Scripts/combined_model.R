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
run_glmulti <- FALSE
do_plot <- FALSE


###########
#Acc Funcs#
###########
excess <- function(x){
  out <- log10(x) - mean(log10(x), na.rm = TRUE)
  return(out)
}

#########################
#Build combined datasets#
#########################
source("italy_crowding_cross_validation.R")
source("china_crowding_cross_validation.R")

############
#Fit models#
############
italy_reg_data$country <- rep("Italy", nrow(italy_reg_data))
china_reg_data$country <- rep("China", nrow(china_reg_data))

china_model <- lm(log10(intensity) ~ log10(contacts) + log10(mobility), data = china_reg_data)

italy_pred <- predict(object = china_model, newdata = italy_reg_data, type = "response")

spear <- spearman.ci(log10(italy_reg_data$intensity), italy_pred, nrep = 1000, conf.level = 0.95)
spear

mod_out_it_cha <- lm(log10(italy_reg_data$intensity)~italy_pred)
R2_out_it_cha <- round(summary(mod_out_it_cha)$adj.r.squared, 2)

quartz()
plot(italy_pred, log10(italy_reg_data$intensity), pch = 16, bty = "n", xlab = "Predicted from China-based model", ylab = "Observed", main = paste0("Obs. vs. pred. from a China-based model intensity in Italy: R2 = ", R2_out_it_cha), ylim = c(-2, 0), xlim = c(-2, 0))
abline(0, 1, col = "red", lty = 3, lwd = 3)
text(-0.5, -.68, "y = x line")

quartz(width = 10, height = 5)
layout(matrix(1:2, nrow = 1))
plot(excess(china_reg_data$crowding), excess(china_reg_data$intensity), pch = 16, col = "#00000090", bty = "n", main = "China", xlab = "Excess log10 Lloyd's mean crowding", ylab = "Excess log10 inverse Shannon entropy")
abline(lm(excess(china_reg_data$intensity) ~ excess(china_reg_data$crowding)), lty = 3, lwd = 3)

plot(excess(italy_reg_data$crowding), excess(italy_reg_data$intensity), pch = 16, col = "#00000090", bty = "n", main = "Italy", xlab = "Excess log10 Lloyd's mean crowding", ylab = "Excess log10 inverse Shannon entropy")
abline(lm(excess(italy_reg_data$intensity) ~ excess(italy_reg_data$crowding)), lty = 3, lwd = 3)


quartz(width = 10, height = 5)
layout(matrix(1:2, nrow = 1))
plot(excess(china_reg_data$patchiness), excess(china_reg_data$intensity), pch = 16, col = "#00000090", bty = "n", main = "China", xlab = "Excess log10 deviation from Poisson (patchiness)", ylab = "Excess log10 inverse Shannon entropy")
abline(lm(excess(china_reg_data$patchiness) ~ excess(china_reg_data$patchiness)), lty = 3, lwd = 3)

plot(excess(italy_reg_data$patchiness), excess(italy_reg_data$intensity), pch = 16, col = "#00000090", bty = "n", main = "Italy", xlab = "Excess log10 deviation from Poisson (patchiness)", ylab = "Excess log10 inverse Shannon entropy")
abline(lm(excess(italy_reg_data$intensity) ~ excess(italy_reg_data$patchiness)), lty = 3, lwd = 3)


quartz(width = 10, height = 10)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
plot(log10(china_reg_data$prop_peak), log10(china_reg_data$intensity), pch = 16, col = "#00000090", bty = "n", main = "China", xlab = "Proportion of cases (peak +/- 1 day) log10", ylab = "Inverse Shannon entropy (log10")

plot(log10(china_reg_data$range), log10(china_reg_data$intensity), pch = 16, col = "#00000090", bty = "n", main = "China", xlab = "Proportion of days above 50th percentile relative to peak height (log10)", ylab = "Inverse Shannon entropy (log10")

plot(log10(italy_reg_data$prop_peak), log10(italy_reg_data$intensity), pch = 16, col = "#00000090", bty = "n", main = "Italy", xlab = "Proportion of cases (peak +/- 1 day) log10", ylab = "Inverse Shannon entropy (log10")

plot(log10(italy_reg_data$range), log10(italy_reg_data$intensity), pch = 16, col = "#00000090", bty = "n", main = "Italy", xlab = "Proportion of days above 50th percentile relative to peak height (log10)", ylab = "Inverse Shannon entropy (log10")
