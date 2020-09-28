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
library(jtools)
library(glmnet)
library(weathermetrics)

#########
#Globals#
#########
run_glmulti <- FALSE
do_plot <- FALSE
do_kelvin <- FALSE
do_z_score <- TRUE

models <- list("Model 1" = c("temp"), "Model 2" = c("humidity"), "Model 3" = c("contacts"), "Model 4" = c("pop"),  "Model 5" = c("mobility"), "Model 6" = c("contacts", "pop"), "Model 7" = c("contacts", "pop", "mobility"), "Model 8" = c("humidity", "contacts", "pop", "mobility"), "Model 9" = c("temp","humidity", "contacts", "pop", "mobility"), "Model 10" = c("temp","humidity", "crowding", "mobility"))

###########
#Acc Funcs#
###########


#########################
#Build combined datasets#
#########################
source("italy_crowding_cross_validation.R")
source("china_crowding_cross_validation.R")

if(do_kelvin == TRUE){
  china_reg_data$temp_f <- china_reg_data$temp
  italy_reg_data$temp_f <- italy_reg_data$temp
  
  china_reg_data$temp <- fahrenheit.to.kelvin(china_reg_data$temp)
  italy_reg_data$temp <- fahrenheit.to.kelvin(italy_reg_data$temp)
  
  file.china <- "model_table_china_kelvin.csv"
  file.italy <- "model_table_italy_kelvin.csv"
  file.combined <- "model_table_combined_kelvin.csv"
}else{
  file.china <- "model_table_china.csv"
  file.italy <- "model_table_italy.csv"
  file.combined <- "model_table_combined.csv"
}

if(do_z_score == TRUE){
  china_reg_data$temp <- (china_reg_data$temp - mean(china_reg_data$temp, na.rm = TRUE))/sd(china_reg_data$temp, na.rm = TRUE)
  italy_reg_data$temp <- (italy_reg_data$temp - mean(italy_reg_data$temp, na.rm = TRUE))/sd(italy_reg_data$temp, na.rm = TRUE)
  temp_combined <- c(china_reg_data$temp, italy_reg_data$temp)
  temp_combined <- (temp_combined - mean(temp_combined, na.rm = TRUE))/sd(temp_combined, na.rm = TRUE)
  file.china <- "model_table_china_zscore.csv"
  file.italy <- "model_table_italy_zscore.csv"
  file.combined <- "model_table_combined_zscore.csv"
}

for(i in 3:ncol(china_reg_data)){
  china_reg_data[,i] <- log10(china_reg_data[,i])
}

for(i in 3:ncol(italy_reg_data)){
  italy_reg_data[,i] <- log10(italy_reg_data[,i])
}

italy_reg_data$country <- rep("Italy", nrow(italy_reg_data))
china_reg_data$country <- rep("China", nrow(china_reg_data))

dat.combined <- rbind(china_reg_data, italy_reg_data)
dat.combined$temp <- log10(temp_combined)

############
#Fit models#
############
results_china <- matrix(NA, ncol = length(models), nrow = 25)
colnames(results_china) <- names(models)
rownames(results_china) <- c("(Intercept)", "(Intercept)-CI", "(Intercept)-p", "temp", "temp-CI", "temp-p", "humidity", "humidity-CI", "humidity-p", "contacts", "contacts-CI", "contacts-p", "pop", "pop-CI", "pop-p", "mobility", "mobility-CI", "mobility-p","crowding", "crowding-CI", "crowding-p", "N", "AIC", "BIC", "R2")
results_china <- as.data.frame(results_china)
models_china <- list()

results_italy <- results_china
models_italy <- list()

results_combined <- results_china
models_combined <- list()

Y.china <- china_reg_data$intensity
Y.italy <- italy_reg_data$intensity
Y.combined <- c(Y.china, Y.italy)

for(i in 1:length(models)){
  #china
  dat.i.china <- data.frame(china_reg_data[,models[[i]]])
  colnames(dat.i.china) <- models[[i]]
  mod.i.china <- lm(Y.china ~ . , data = dat.i.china)
  models_china[[i]] <- mod.i.china
  mt_coefs.i.china <- match(names(coefficients(mod.i.china)), rownames(results_china))
  results_china[mt_coefs.i.china,i] <- round(coefficients(mod.i.china), 3)
  
  for(k in 1:length(mt_coefs.i.china)){
    results_china[mt_coefs.i.china[k]+1,i] <- paste0("[", round(coefficients(mod.i.china)[k]-1.96*summary(mod.i.china)$coefficients[k,"Std. Error"], 3), ",", round(coefficients(mod.i.china)[k]+1.96*summary(mod.i.china)$coefficients[k,"Std. Error"], 3), "]")
  }
  
  results_china[mt_coefs.i.china+2,i] <- round(summary(mod.i.china)$coefficients[,"Pr(>|t|)"], 4)
  results_china["N",i] <- length(mod.i.china$residuals)
  results_china["AIC",i] <- round(aic(mod.i.china), 2)
  results_china["BIC",i] <- round(bic(mod.i.china), 2)
  results_china["R2",i] <- round(summary(mod.i.china)$adj.r.squared, 2)
  
  #italy
  dat.i.italy <- data.frame(italy_reg_data[,models[[i]]])
  colnames(dat.i.italy) <- models[[i]]
  mod.i.italy <- lm(Y.italy ~ . , data = dat.i.italy)
  models_italy[[i]] <- mod.i.italy
  mt_coefs.i.italy <- match(names(coefficients(mod.i.italy)), rownames(results_italy))
  results_italy[mt_coefs.i.italy,i] <- round(coefficients(mod.i.italy), 3)
  
  for(k in 1:length(mt_coefs.i.italy)){
    results_italy[mt_coefs.i.italy[k]+1,i] <- paste0("[", round(coefficients(mod.i.italy)[k]-1.96*summary(mod.i.italy)$coefficients[k,"Std. Error"], 3), ",", round(coefficients(mod.i.italy)[k]+1.96*summary(mod.i.italy)$coefficients[k,"Std. Error"], 3), "]")
  }
  
  results_italy[mt_coefs.i.italy+2,i] <- round(summary(mod.i.italy)$coefficients[,"Pr(>|t|)"], 4)
  results_italy["N",i] <- length(mod.i.italy$residuals)
  results_italy["AIC",i] <- round(aic(mod.i.italy), 2)
  results_italy["BIC",i] <- round(bic(mod.i.italy), 2)
  results_italy["R2",i] <- round(summary(mod.i.italy)$adj.r.squared, 3)
  
  #combined
  dat.i.combined <- data.frame(dat.combined[,models[[i]]])
  colnames(dat.i.combined) <- models[[i]]
  mod.i.combine <- lm(Y.combined ~ . , data = dat.i.combined)
  models_combined[[i]] <- mod.i.combine
  mt_coefs.i.combine <- match(names(coefficients(mod.i.combine)), rownames(results_combined))
  results_combined[mt_coefs.i.combine,i] <- round(coefficients(mod.i.combine), 3)
  
  for(k in 1:length(mt_coefs.i.combine)){
    results_combined[mt_coefs.i.combine[k]+1,i] <- paste0("[", round(coefficients(mod.i.combine)[k]-1.96*summary(mod.i.combine)$coefficients[k,"Std. Error"], 3), ",", round(coefficients(mod.i.combine)[k]+1.96*summary(mod.i.combine)$coefficients[k,"Std. Error"], 3), "]")
  }
  
  results_combined[mt_coefs.i.combine+2,i] <- round(summary(mod.i.combine)$coefficients[,"Pr(>|t|)"], 4)
  results_combined["N",i] <- length(mod.i.combine$residuals)
  results_combined["AIC",i] <- round(aic(mod.i.combine), 2)
  results_combined["BIC",i] <- round(bic(mod.i.combine), 2)
  results_combined["R2",i] <- round(summary(mod.i.combine)$adj.r.squared, 3)
}

write.csv(results_china, file = file.china)
write.csv(results_italy, file = file.italy)
write.csv(results_combined, file = file.combined)

#glmmulti
glmulti.combined <-
  glmulti(intensity ~ patchiness + mobility + contacts + pop + temp + humidity, data = dat.combined,
          level = 1,               # 2 = pairwise interaction considered
          method = "h",            # h = Exhaustive approach
          crit = "bic",            # bic = BIC as criteria
          confsetsize = 1,        # Keep # best models
          plotty = F, report = F,  # Plot or interim reports
          fitfunction = lm)
best.mod.combine <- slot(object = glmulti.combined, name = "objects")[[1]]

#glmnet
# A tuning parameter for the lasso/elastic-net model-selection algorithm
# alpha = 1 corresponds to l1/lasso
# alpha = 0.5 corresponds to elastic net
# alpha = 0 gives ridge regression, which will not do any model selection, merely shrinkage
alpha <- 0.5
D <- dat.combined[,models[[length(models)]]]
D <- cbind(Y.combined, D)
D <- as.matrix(D)
D <- na.omit(D)
Y_cv <- D[,"Y.combined"]
D <- D[,-which(colnames(D) == "Y.combined")]

#pick lambda by loo cross-validation
mycv <- cv.glmnet(D, Y_cv, nfolds=length(Y_cv), parallel=TRUE, family='gaussian', alpha=alpha, keep=TRUE, standardize=TRUE)
lambda0 <- mycv$lambda.min

# Fit the model
glm.gaus <- glmnet(D, Y_cv, family='gaussian', lambda = lambda0, alpha=alpha,  maxit=1000000, standardize=TRUE)
glm.gaus$beta
