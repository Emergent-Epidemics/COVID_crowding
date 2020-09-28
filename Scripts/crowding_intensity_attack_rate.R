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
plot(log10(china_reg_data$intensity), log10(china_reg_data$total_cases/(china_reg_data$pop)), xlab = "Intensity (log10)", ylab = "Final attack rate (log10)", main = "China: Intensity and Attack Rate", pch = 16, bty = "n", col = "#00000090")
attack_china <- log10(china_reg_data$total_cases/(china_reg_data$pop))
mod_att_int_china <- lm(attack_china ~ log10(china_reg_data$intensity))
abline(mod_att_int_china, col = "#e41a1c", lty = 3, lwd = 2)

cor.test(attack_china, log10(china_reg_data$intensity))
spearman.ci(attack_china, log10(china_reg_data$intensity), nrep = 1000, conf.level = 0.95)

quartz()
layout(matrix(1:4, ncol = 2, byrow = TRUE))
plot(log10(italy_reg_data$intensity), log10(italy_reg_data$peak_cases/italy_reg_data$total_cases), xlab = "Intensity log10", ylab = "Prop cases peak (+/- 1 day) log10", main = "Italy")

plot(log10(italy_reg_data$intensity), log10(italy_reg_data$total_cases/italy_reg_data$pop), xlab = "Intensity log10", ylab = "Final attack rate log10", main = "Italy")

plot(log10(china_reg_data$intensity), log10(china_reg_data$peak_cases/china_reg_data$total_cases), xlab = "Intensity log10", ylab = "Prop cases peak (+/- 1 day) log10", main = "China")

plot(log10(china_reg_data$intensity), log10(china_reg_data$total_cases/china_reg_data$pop), xlab = "Intensity log10", ylab = "Final attack rate log10", main = "China")

#extended data figure 9
quartz()
layout(matrix(1:2, nrow = 1))
plot(log10(china_reg_data$crowding)-log10(mean(china_reg_data$crowding, na.rm = TRUE)), log10(china_reg_data$intensity)-log10(mean(china_reg_data$intensity, na.rm = TRUE)), ylab = "Excess log10 intensity", xlab = "Excess log10 crowding", main = "China", pch = 16)

data.c <- data.frame(log10(china_reg_data$crowding)-log10(mean(china_reg_data$crowding, na.rm = TRUE)), log10(china_reg_data$intensity)-log10(mean(china_reg_data$intensity, na.rm = TRUE)))
colnames(data.c) <- c("excess_crowding", "excess_int")
mod.c.excess <- lm(excess_int ~ excess_crowding, data = data.c)
abline(mod.c.excess)

plot(log10(italy_reg_data$crowding)-log10(mean(italy_reg_data$crowding, na.rm = TRUE)), log10(italy_reg_data$intensity)-log10(mean(italy_reg_data$intensity, na.rm = TRUE)), ylab = "Excess log10 intensity", xlab = "Excess log10 crowding", main = "Italy", pch = 16)

data.i <- data.frame(log10(italy_reg_data$crowding)-log10(mean(italy_reg_data$crowding, na.rm = TRUE)), log10(italy_reg_data$intensity)-log10(mean(italy_reg_data$intensity, na.rm = TRUE)))
colnames(data.i) <- c("excess_crowding", "excess_int")
mod.it.excess <- lm(excess_int ~ excess_crowding, data = data.i)
abline(mod.it.excess)


quartz()

colors <- colorRampPalette(c("#feb24c", "#fc4e2a", "#800026"))
cols <- c("gray","gray","gray",colors(4))

layout(matrix(1:2, nrow = 1))
plot(log10(china_reg_data$patchiness)-log10(china_reg_data$crowding), log10(china_reg_data$intensity)-log10(mean(china_reg_data$intensity, na.rm = TRUE)), xlab = "Patchiness log10 - crowding log10", ylab = "Excess log10 intensity", main = "China", pch = 16, col = cols[round(log10(china_reg_data$crowding),0)])

data.c <- data.frame(log10(china_reg_data$patchiness)-log10(china_reg_data$crowding), log10(china_reg_data$intensity)-log10(mean(china_reg_data$intensity, na.rm = TRUE)))
colnames(data.c) <- c("crowding_m_patch", "excess_int")
abline(lm(excess_int ~ crowding_m_patch, data = data.c))

plot(log10(italy_reg_data$patchiness)-log10(italy_reg_data$crowding), log10(italy_reg_data$intensity)-log10(mean(italy_reg_data$intensity, na.rm = TRUE)), xlab = "Patchiness log10 - crowding log10", ylab = "Excess log10 intensity", main = "Italy", pch = 16, col = cols[round(log10(italy_reg_data$crowding),0)])

data.i <- data.frame(log10(italy_reg_data$patchiness)-log10(italy_reg_data$crowding), log10(italy_reg_data$intensity)-log10(mean(italy_reg_data$intensity, na.rm = TRUE)))
colnames(data.i) <- c("crowding_m_patch", "excess_int")
abline(lm(excess_int ~ crowding_m_patch, data = data.i))

##########################
#Crowding and attack rate#
##########################
quartz(width = 10, height = 5)
layout(matrix(1:2, nrow = 1))
mod.attack.italy <- lm(log(total_cases/pop) ~ log(intensity), data = italy_reg_data)
R2_italy <- round(summary(mod.attack.italy)$adj.r.squared, 2)
plot(log(italy_reg_data$intensity), log(italy_reg_data$total_cases/italy_reg_data$pop), pch = 16, bty = "n", ylab = "Proportion infected", xlab = "Epidemic intensity", main = paste0("Italy: R2 ~ ", R2_italy))
abline(mod.attack.italy, lwd = 3, col = "red", lty = 3)

summary(lm(log10(total_cases/pop) ~ log10(pop) + log10(crowding) + log10(mobility) + log10(temp) + log10(humidity) + log10(intensity), data = italy_reg_data))

mod.attack.china <- lm(log(total_cases/(pop)) ~ log(intensity), data = china_reg_data)
R2_china <- round(summary(mod.attack.china)$adj.r.squared, 2)
plot(log(china_reg_data$intensity), log(china_reg_data$total_cases/(china_reg_data$pop)), pch = 16, bty = "n", ylab = "Proportion infected", xlab = "Epidemic intensity", main = paste0("China: R2 ~ ", R2_china))
abline(mod.attack.china, lwd = 3, col = "red", lty = 3)

summary(lm(log10(total_cases/(pop)) ~ log10(pop) + log10(crowding) + log10(mobility) + log10(temp) + log10(humidity) + log10(intensity), data = china_reg_data))

if(run_glmulti == TRUE){
  glmulti.italy.attack <-
    glmulti(log10(total_cases/pop) ~ log10(pop) + log10(crowding) + log10(mobility) + log10(temp) + log10(humidity) + log10(intensity), data = italy_reg_data,
            level = 2,               # pairwise interaction considered
            method = "h",            # Exhaustive approach
            crit = "bic",            # BIC as criteria
            confsetsize = 1,        # Keep 10 best models
            plotty = F, report = F,  # Plot or interim reports
            fitfunction = glm)
  
  glmulti.china.attack <-
    glmulti(log10(total_cases/(pop)) ~ log10(pop) + log10(crowding) + log10(mobility) + log10(temp) + log10(humidity), data = china_reg_data,
            level = 2,               # pairwise interaction considered
            method = "h",            # Exhaustive approach
            crit = "bic",            # BIC as criteria
            confsetsize = 1,        # Keep 10 best models
            plotty = F, report = F,  # Plot or interim reports
            fitfunction = glm)
}


###############################
#Crowding and PEAK attack rate#
###############################

quartz(width = 10, height = 5)
layout(matrix(1:2, nrow = 1))
#china
mod.attack.china <- lm(log10(peak_cases) ~ log10(total_cases), data = china_reg_data)
mod.attack.china.resid <- lm(mod.attack.china$residuals ~ log10(china_reg_data$intensity))
R2_china <- round(summary(mod.attack.china.resid)$adj.r.squared, 2)

plot(log10(china_reg_data$intensity), mod.attack.china$residuals, pch = 16, bty = "n", ylab = "Excess peak cases", xlab = "Inverse Shannon entropy (log10)", main = paste0("China: R2 ~ ", R2_china), col = "#00000095")
abline(mod.attack.china.resid, lwd = 3, col = "#b2182b", lty = 3)

summary(lm(log10(peak_cases) ~ log10(intensity)+log10(total_cases), data = china_reg_data))
summary(lm(log10(peak_cases) ~ log10(crowding)+log10(pop)+log10(mobility)+log10(humidity)+log10(temp), data = china_reg_data))

#italy
mod.attack.italy <- lm(log10(peak_cases) ~ log10(total_cases), data = italy_reg_data)
mod.attack.italy.resid <- lm(mod.attack.italy$residuals ~ log10(italy_reg_data$intensity[-mod.attack.italy$na.action]))
R2_italy <- round(summary(mod.attack.italy.resid)$adj.r.squared, 2)

plot(log10(italy_reg_data$intensity[-mod.attack.italy$na.action]), mod.attack.italy$residuals, pch = 16, bty = "n", ylab = "Excess peak cases", xlab = "Inverse Shannon entropy (log10)", main = paste0("Italy: R2 ~ ", R2_italy), col = "#00000095")
abline(mod.attack.italy.resid, lwd = 3, col = "#b2182b", lty = 3)
summary(lm(log10(peak_cases) ~ log10(intensity)+log10(total_cases), data = italy_reg_data))
summary(lm(log10(peak_cases) ~ log10(crowding)+log10(pop)+log10(mobility)+log10(humidity)+log10(temp), data = italy_reg_data))

if(run_glmulti == TRUE){
  glmulti.italy.peak.attack <-
    glmulti(log10(peak_cases) ~ log10(pop) + log10(crowding) + log10(mobility) + log10(temp) + log10(humidity), data = italy_reg_data,
            level = 1,               # pairwise interaction considered
            method = "h",            # Exhaustive approach
            crit = "bic",            # BIC as criteria
            confsetsize = 1,        # Keep 10 best models
            plotty = F, report = F,  # Plot or interim reports
            fitfunction = glm)
  
  glmulti.china.peak.attack <-
    glmulti(log10(peak_cases) ~ log10(pop) + log10(crowding) + log10(mobility) + log10(temp) + log10(humidity), data = china_reg_data,
            level = 1,               # pairwise interaction considered
            method = "h",            # Exhaustive approach
            crit = "bic",            # BIC as criteria
            confsetsize = 1,        # Keep 10 best models
            plotty = F, report = F,  # Plot or interim reports
            fitfunction = glm)
}
