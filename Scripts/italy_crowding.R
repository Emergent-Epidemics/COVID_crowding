#SV Scarpino
#May 2020
#Crowding analysis of Italy

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
library(geosphere)

#########
#Globals#
#########
do_plot <- FALSE

###########
#Acc Funcs#
###########
lloyd <- function(qs){
  num <- sum((qs-1)*qs, na.rm = TRUE)
  denom <- sum(qs, na.rm = TRUE)
  out <- num/denom
  return(out)
}

reporting_days <- function(x){
  start <- which(x != 0)[1]
  x.filt <- x[start:length(x)]
  stop <- rev(which(x.filt != 0))[1]
  x.filt.2 <- x.filt[1:stop]
  return(x.filt.2)
}

range_50 <- function(x){
  med <- median(x)
  min.med <- min(which(x > med))
  max.med <- max(which(x > med))
  range_out <- max.med - min.med
  return(range_out)
}

######
#Data#
######
italy <- readOGR("../Data/Italy_shapefile/gadm36_ITA_2.shp")
italy_cases <- read.csv("../Data/dpc-covid19-ita-province.csv")
italy_cases$data <- as.POSIXct(italy_cases$data)
italy_data <- readRDS("../Data/italy_admin2_mobility_reductions_5_6_20.RDS")
italy_data$NAME_2 <- as.character(italy_data$NAME_2)
italy_covariates <- read.csv("../Data/italy_admin2_covariates_5_7_20.csv")

cols <- colorRampPalette(colors = c("#b2182b", "#f7f7f7", "#2166ac"))
colors <- cols(1000)

##########
#Mod data#
##########
#lining up the shapefile and case file
italy_cases$denominazione_provincia[which(italy_cases$denominazione_provincia == "Firenze")]  <- "Florence"
italy_cases$denominazione_provincia[which(italy_cases$denominazione_provincia == "Forlì-Cesena")]  <- "Forli' - Cesena"
italy_cases$denominazione_provincia[which(italy_cases$denominazione_provincia == "Mantova")]  <- "Mantua"
italy_cases$denominazione_provincia[which(italy_cases$denominazione_provincia == "Monza e della Brianza")]  <- "Monza and Brianza"
italy_cases$denominazione_provincia[which(italy_cases$denominazione_provincia == "Pesaro e Urbino")]  <- "Pesaro E Urbino"
italy_cases$denominazione_provincia[which(italy_cases$denominazione_provincia == "Siracusa")]  <- "Syracuse"
italy_cases$denominazione_provincia[which(italy_cases$denominazione_provincia == "Padova")]  <- "Padua"
italy_cases$denominazione_provincia[which(italy_cases$denominazione_provincia == "Reggio di Calabria")]  <- "Reggio Di Calabria"
italy_cases$denominazione_provincia[which(italy_cases$denominazione_provincia == "Reggio nell'Emilia")]  <- "Reggio Nell'Emilia"

italy$NAME_2[which(italy$NAME_2 %in% c("Carbonia-Iglesias", "Medio Campidano"))] <- "Sud Sardegna"
italy_data$NAME_2[which(italy_data$NAME_2 %in% c("Carbonia-Iglesias", "Medio Campidano"))] <- "Sud Sardegna"

mt_name <- match(italy_covariates$GID_2, italy_data$GID_2)
italy_covariates$NAME_2 <- italy_data$NAME_2[mt_name]
###############
#Modify raster#
###############
raster_files <- list.files("../Data/Italy_layers", full.names = TRUE)
pops <- list()
for(i in 1:length(raster_files)){
  name.ras.i <- strsplit(x = as.character(raster_files[i]), "../Data/Italy_layers/")[[1]][2]
  pops[[name.ras.i]] <- raster(raster_files[i])
}

##########
#Analysis#
##########
prefects <- unique(italy$NAME_2)

results <- list()

pb <- txtProgressBar(1, length(pops), style=3)

for(j in 1:length(pops)){
  use.j <- names(pops)[j]
  crowding <- rep(NA, length(prefects))
  total_pops <- rep(NA, length(prefects))
  ent <- rep(NA, length(prefects))
  dists <- rep(NA, length(prefects))
  missed <- c()
  
  for(i in 1:length(prefects)){
    use.i <- which(italy$NAME_2 == prefects[i])
    poly_coords.i <- coordinates(italy[use.i,])
    int.i <- intersect(pops[[use.j]], italy[use.i,])
    pop_points.i <- rasterToPoints(int.i)
    
    if(nrow(pop_points.i) == 0){
      crowding[i] <- NA
      total_pops[i] <- NA
      ent[i] <- NA
      dists[i] <- NA
      missed <- c(missed, as.character(prefects[i]))
      next
    }
    
    pop_points.mat.i <- matrix(pop_points.i[,1:2], ncol = 2, byrow = FALSE)
    coords.i <- SpatialPoints(pop_points.mat.i)
    proj4string(coords.i) <- proj4string(pops[[use.j]])
    over.i <- over(coords.i, italy[use.i,])
    use.pop.i <- which(is.na(over.i[,1]) == FALSE)
    
    centroid.i <- gCentroid(coords.i)
    dists.i <- distm(c(12.4964, 41.9028), centroid.i) #distance from Rome
    if(coordinates(centroid.i)[2] < 41.9028){
      dists.i <- -1 * dists.i
    }
    dists[i] <- dists.i
    
    if(do_plot == TRUE){
      brk <- seq(1, log(650000), length.out = 10)
      arg <- list(at=brk, labels=round(exp(brk)))
      quartz()
      plot(log(int.i+1), col = rev(colors), interpolate = TRUE, main = italy$NAME_2[i], zlim=c(0,log(650000)), axis.args=arg)
      plot(italy[use.i,], add = TRUE)
    }
    
    #lloyd.i <- lloyd(qs = pop_points.i[use.pop.i,3])
    lloyd.i <- agg_index(pop_points.i[use.pop.i,3], method = "lloyd", type = "mean-crowding")
    if(length(use.pop.i) == 0){
      crowding[i] <- NA
      total_pops[i] <- NA
      ent[i] <- NA
      missed <- c(missed, as.character(prefects[i]))
    }else{
      crowding[i] <- lloyd.i[[1]]
      total_pops[i] <- sum(pop_points.i[use.pop.i,3], na.rm = TRUE)
      if(length(use.pop.i) == 1){
        ent[i] <- NA
      }else{
        ent[i] <- entropy(discretize(pop_points.i[use.pop.i,3], numBins = 5))
      }
    }
  }
  results[[use.j]] <- list("crowding" = crowding, "total_pops" = total_pops, "ent" = ent, "dists" = dists, "missed" = missed)
  setTxtProgressBar(pb, j)
}


int <- rep(NA, length(prefects))
total <- rep(NA, length(prefects))
missed <- c()
for(i in 1:length(prefects)){
  #intensity and total
  use.int.i <- which(italy_cases$denominazione_provincia == prefects[i])
  if(length(use.int.i) == 0){
    missed <- c(missed, prefects[i])
    next
  }
  ord.int.i <- order(italy_cases$data[use.int.i])
  new.cases.i <- diff(italy_cases$totale_casi[use.int.i][ord.int.i])
  rep_days.i <- reporting_days(new.cases.i)
  new.cases.i <- rep_days.i
  p.i <- new.cases.i/sum(new.cases.i, na.rm = TRUE) #calculate fraction of total cases for each week
  epidemic_intensity.i <- entropy(p.i)^-1 #calculate inverse Shannon entropy
  min_int.i <- entropy(rep(1/length(p.i), length(p.i)))^-1
  epidemic_intensity_norm.i <- epidemic_intensity.i - min_int.i
  int[i] <- epidemic_intensity_norm.i
  total[i] <- sum(new.cases.i, na.rm = TRUE)
}

mt_reg <- match(prefects, italy_covariates$NAME_2)


R2s <- rep(NA, length(results))
coefs <- rep(NA, length(results))
coef_crowd <- rep(NA, length(results))
p_crowd <- rep(NA, length(results))
rhos <- matrix(NA, ncol = 3, nrow = length(results))
colnames(rhos) <- c("rho", "lower_95", "upper_95")
rhos <- as.data.frame(rhos)

for(i in 1:length(results)){
  data.reg.i <- data.frame(prefects, results[[i]]$crowding, int, results[[i]]$total_pops)
  colnames(data.reg.i) <- c("location", "crowding", "intensity", "pop")
  
  data.reg.i$temp <- italy_covariates$temp[mt_reg]
  data.reg.i$humidity <- italy_covariates$spechum[mt_reg]
  data.reg.i$mobility <- italy_covariates$within_april_mean_over_december[mt_reg]
  
  data.reg.i$intensity[which(data.reg.i$intensity > 1/log(2))] <- NA
  data.reg.i$intensity[which(data.reg.i$intensity == 0)] <- data.reg.i$intensity[which(data.reg.i$intensity == 0)] + min(data.reg.i$intensity[which(data.reg.i$intensity != 0)], na.rm = TRUE) * 0.98
  data.reg.i$intensity <- data.reg.i$intensity/(1/log(2))
  
  if(nrow(na.omit(data.reg.i)) == 0){
    next
  }
  
  mod.i <- lm(log10(intensity) ~ log10(pop) + log10(crowding) + log10(mobility) + log10(humidity) + log10(temp), data = data.reg.i)
  r2.i <- summary(mod.i)$adj.r.squared
  slope.i <- coefficients(mod.i)["log10(crowding)"]
  p.i <- summary(mod.i)$coefficients["log10(crowding)","Pr(>|t|)"]
  spear.i <- spearman.ci(log10(data.reg.i$intensity), log10(data.reg.i$crowding), nrep = 1000, conf.level = 0.95)
  
  rhos$rho[i] <- spear.i$estimate
  rhos[i, c("lower_95", "upper_95")] <- spear.i$conf.int
  R2s[i] <- r2.i
  coef_crowd[i] <- slope.i
  p_crowd[i] <- p.i
}

scale <- unlist(strsplit(names(results), ".tif"))
scale <- unlist(strsplit(scale, "pop_italy_"))
scale <- unlist(strsplit(scale, "k"))
scale <- as.numeric(scale)

ord <- order(scale)

quartz(width = 14, height = 7)
layout(matrix(1:3, nrow = 1))
plot(scale[ord], R2s[ord], type = "b", main = "R2 of univariate regression and intensity: Italy", ylab = "R2", xlab = "Spatial aggregation", bty = "n", pch = 16)
plot(scale[ord], coef_crowd[ord], type = "b", main = "Coefficient of crowding: Italy", ylab = "Coefficient of crowding", xlab = "Spatial aggregation", bty = "n", pch = 16, ylim = c(min(coef_crowd[ord], na.rm = TRUE), 0))
abline(h = 0, lty = 3, col = "red")
plot(scale[ord], log10(p_crowd[ord]), type = "b", main = "p value crowding: Italy", ylab = "p value", xlab = "Spatial aggregation", bty = "n", pch = 16, yaxt = "n", ylim = c(-5, 0))
at <- c(-5, -4, -3, -2, -1, 0)
axis(2, at = at, labels = round(10^at, 4))
abline(h = log10(0.05), lty = 3, col = "red")

rownames(rhos) <- scale
rhos <- rhos[ord,]