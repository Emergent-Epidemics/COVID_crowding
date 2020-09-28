#SV Scarpino
#May 2020
#Crowding intensity and attack rate analysis of China 

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
do_country <- "italy"

###########
#Acc Funcs#
###########


#########################
#Build combined datasets#
#########################
source(paste0(do_country, "_crowding_cross_validation.R"))

##########
#Analysis#
##########
var_diff[which(is.finite(var_diff) == FALSE)] <- NA
int[which(is.finite(int) == FALSE)] <- NA
int[which(int > 1/log(2))] <- NA
int[which(int == 0)] <- int[which(int == 0)] + min(int[which(int != 0)], na.rm = TRUE) * 0.98
int <- int/(1/log(2))
ranges[which(is.finite(ranges) == FALSE)] <- NA

resid_int <- log10(int) - mean(log10(int), na.rm = TRUE)
resid_var <- log10(var_diff + 1) - mean(log10(var_diff + 1), na.rm = TRUE)
resid_crowd <- log10(crowding) - mean(log10(crowding), na.rm = TRUE)
resid_patch <- log10(patchy) - mean(log10(patchy), na.rm = TRUE)
resid_peak <- lm(log10(peak_cases) ~ log10(total))$residuals
resid_range <- log10(ranges+1) - mean(log10(ranges+1), na.rm = TRUE)

quartz(width = 10, height = 5)
layout(matrix(1:3, nrow = 1, byrow = TRUE))
plot(resid_crowd, resid_int, pch = 16, col = "#00000095", bty = "n", xlab = "Residual Lloyd's mean crowding", ylab = "Residual inverse Shannon entropy")
abline(lm(resid_int ~ resid_crowd), lty = 3, lwd = 3)

plot(resid_var, resid_int, pch = 16, col = "#00000095", bty = "n", xlab = "Residual stationary variance", ylab = "Residual inverse Shannon entropy")
abline(lm(resid_int ~ resid_var), lty = 3, lwd = 3)

plot(resid_patch, resid_int, pch = 16, col = "#00000095", bty = "n", xlab = "Residual deviation from Poisson", ylab = "Residual inverse Shannon entropy")
abline(lm(resid_int ~ resid_patch), lty = 3, lwd = 3)

#Group analysis  
group <- rep(NA, length(resid_int))
group[which(resid_int < 0 & resid_var > 0)] <- "A"
group[which(resid_int > 0 & resid_var > 0)] <- "B"
group[which(resid_int < 0 & resid_var < 0)] <- "C"
group[which(resid_int > 0 & resid_var < 0)] <- "D"

cols <- c("#b35806", "#d73027", "#4575b4", "#542788")
cols_list <- list("A" = "#b35806", "B" = "#d73027", "C" = "#4575b4", "D" = "#542788")

dat.plot <- data.frame(group, resid_int, resid_var)
colnames(dat.plot) <- c("Curve_type", "Excess_intensity", "Excess_variance")
ggplot(dat.plot, aes(x = Excess_intensity, y =  Excess_variance, color = Curve_type)) + geom_point() + scale_color_manual(values =  cols, guide_legend(title = "Outbreak category")) + xlab("Excess log10 intensity") + ylab("Excess log10 variance") + theme(legend.position = c(0.8,0.8), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 15), panel.grid.minor = element_line(colour = "#00000050",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + scale_y_continuous(expand = c(0.01,0.01)) + geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed")

quartz()
layout(matrix(1:4, nrow = 2, byrow = TRUE))
ints <- c()
vars <- c()
for(l in LETTERS[1:4]){
  if(do_country == "italy"){
    i.all <- which(group == l & log10(results[[1]]$total_pops) < 5.5 & log10(results[[1]]$total_pops) > 5.4)
    i <- i.all[1]
    use.int.i <- which(italy_cases$denominazione_provincia == prefects[i])
    ord.int.i <- order(italy_cases$data[use.int.i])
    new.cases.i <- diff(italy_cases$totale_casi[use.int.i][ord.int.i])
  }
  
  if(do_country == "china"){
    i.all <- which(group == l & log10(results[[1]]$total_pops) < 6.5 & log10(results[[1]]$total_pops) > 6.4)
    i <- i.all[1]
    objectid.i <- china$OBJECTID[which(china$ADM2_PCODE %in% prefects[i])]
    use.int.i <- which(china_cases$OBJECTID == objectid.i)
    ord.int.i <- order(china_cases$date[use.int.i])
    new.cases.i <- china_cases$Cases[use.int.i][ord.int.i]
  }
  p.i <- new.cases.i/sum(new.cases.i, na.rm = TRUE) #calculate fraction of total cases for each week
  ints <- c(ints, entropy(p.i)^-1)
  vars <- c(vars, var(diff(new.cases.i), na.rm = TRUE))
  plot(new.cases.i, type = "b", main = paste0(l, ": Intensity = ", round(entropy(p.i)^-1,2), ", Variance = ", round(var(diff(new.cases.i), na.rm = TRUE),2)), xlab = "Days since first case", col = cols_list[[l]], pch = 16, ylab = "Proportion")
}

quartz()
layout(matrix(c(1:4), nrow = 2, byrow = TRUE))
boxplot(log10(total/total_pops) ~ group, main = "Total attack rate", col = cols, xlab = "", ylab = "Total attack rate (log10)")
boxplot(log10(peak_cases/total_pops) ~ group, main = "Peak (+/- 1 day) attack rate", col = cols, xlab = "", ylab = "Peak (+/- 1 day) attack rate (log10)")
boxplot(log10(crowding) ~ group, main = "Inverse Shannon entropy (crowding)", col = cols, xlab = "", ylab = "Inverse Shannon entropy (log10)")
boxplot(log10(total_pops) ~ group, main = "Population size", col = cols, xlab = "", ylab = "Population size (log10)")
