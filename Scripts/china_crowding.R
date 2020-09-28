#SV Scarpino
#April 2020
#Crowding analysis of China

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

######
#Data#
######
china <- readOGR("../Data/China_shapefile/chn_admbnda_adm2_ocha.shp")
data <- readRDS("../Data/admin2_crowding_intensity_4_11.RDS")
curves <- readRDS("../Data/epi_curves_guangzhou_quzhou.RDS")
mobility <- read.csv("../Data/Google COVID-19 Aggregated Mobility Research Dataset.csv")
china_cases <- readRDS("../Data/chinese_prefactures_cases_5_7_20.RDS")
china_cases$date <- as.POSIXct(china_cases$date)
cols <- colorRampPalette(colors = c("#b2182b", "#f7f7f7", "#2166ac"))
colors <- cols(1000)

###############
#Modify raster#
###############
raster_files <- list.files("../Data/China_layers", full.names = TRUE)
pops <- list()
for(i in 1:length(raster_files)){
  name.ras.i <- strsplit(x = as.character(raster_files[i]), "../Data/China_layers/")[[1]][2]
  pops[[name.ras.i]] <- raster(raster_files[i])
}

##########
#Analysis#
##########
prefects <- unique(china$ADM2_PCODE)

results <- list()

pb <- txtProgressBar(1, length(pops), style=3)

for(j in 1:length(pops)){
  use.j <- names(pops)[j]
  crowding <- rep(NA, length(prefects))
  total_pops <- rep(NA, length(prefects))
  ent <- rep(NA, length(prefects))
  missed <- c()
  for(i in 1:length(prefects)){
    use.i <- which(china$ADM2_PCODE == prefects[i])
    poly_coords.i <- coordinates(china[use.i,])
    int.i <- intersect(pops[[use.j]], china[use.i,])
    pop_points.i <- rasterToPoints(int.i)
    
    if(nrow(pop_points.i) == 0){
      crowding[i] <- NA
      total_pops[i] <- NA
      ent[i] <- NA
      missed <- c(missed, as.character(prefects[i]))
      next
    }
    
    pop_points.mat.i <- matrix(pop_points.i[,1:2], ncol = 2, byrow = FALSE)
    coords.i <- SpatialPoints(pop_points.mat.i)
    proj4string(coords.i) <- proj4string(pops[[use.j]])
    over.i <- over(coords.i, china[use.i,])
    use.pop.i <- which(is.na(over.i[,1]) == FALSE)
    
    if(do_plot == TRUE){
      #i <- 291 #Quzhou
      #i <- 57 #Guangzhou
      brk <- seq(1, log(650000), length.out = 10)
      arg <- list(at=brk, labels=round(exp(brk)))
      quartz()
      plot(log(int.i+1), col = rev(colors), interpolate = TRUE, main = china$ADM2_EN[i], zlim=c(0,log(650000)), axis.args=arg)
      plot(china[use.i,], add = TRUE)
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
  results[[use.j]] <- list("crowding" = crowding, "total_pops" = total_pops, "ent" = ent, "missed" = missed)
  setTxtProgressBar(pb, j)
}

int <- rep(NA, length(prefects))
total <- rep(NA, length(prefects))
missed <- c()
for(i in 1:length(prefects)){
  #intensity and total
  objectid.i <- china$OBJECTID[which(china$ADM2_PCODE %in% prefects[i])]
  use.int.i <- which(china_cases$OBJECTID == objectid.i)
  if(length(use.int.i) == 0){
    missed <- c(missed, prefects[i])
    next
  }
  ord.int.i <- order(china_cases$date[use.int.i])
  if(length(which(diff(china_cases$date[use.int.i]) != 1)) > 0){
    stop("dates not sequential by days")
  }
  new.cases.i <- china_cases$Cases[use.int.i][ord.int.i]
  p.i <- new.cases.i/sum(new.cases.i, na.rm = TRUE) #calculate fraction of total cases for each week
  epidemic_intensity.i <- entropy(p.i)^-1 #calculate inverse Shannon entropy
  min_int.i <- entropy(rep(1/length(p.i), length(p.i)))^-1
  epidemic_intensity_norm.i <- epidemic_intensity.i - min_int.i
  int[i] <- epidemic_intensity_norm.i
  total[i] <- sum(new.cases.i, na.rm = TRUE)
}

mt_pre <- match(data$ADM2_PCODE, prefects)
mt_mob <- match(gsub(pattern = "CN", replacement = "", as.character(data$ADM2_PCODE)), mobility$CODE)


R2s <- rep(NA, length(results))
coef_crowd <- rep(NA, length(results))
p_crowd <- rep(NA, length(results))
rhos <- matrix(NA, ncol = 3, nrow = length(results))
colnames(rhos) <- c("rho", "lower_95", "upper_95")
rhos <- as.data.frame(rhos)

for(i in 1:length(results)){
  data.reg.i <- data.frame(prefects[mt_pre], as.character(data$ADM2_EN), results[[i]]$crowding[mt_pre], int[mt_pre], data$temp, data$`Population Size (Millions)`, data$meanhum, as.numeric(as.character(mobility$reduction.in.mobility.during.the.outbreak.period[mt_mob])), total[mt_pre])
  colnames(data.reg.i) <- c("ID","location", "crowding", "intensity", "temp", "pop", "humidity", "mobility", "total_cases")
  
  data.reg.i$intensity[which(data.reg.i$intensity > 1/log(2))] <- NA
  data.reg.i$intensity[which(data.reg.i$intensity == 0)] <- data.reg.i$intensity[which(data.reg.i$intensity == 0)] + min(data.reg.i$intensity[which(data.reg.i$intensity != 0)], na.rm = TRUE) * 0.98
  data.reg.i$intensity <- data.reg.i$intensity/(1/log(2))
  data.reg$pop <- data.reg$pop * 1e6
  
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
scale <- unlist(strsplit(scale, "pop_china_"))
scale <- unlist(strsplit(scale, "k"))
scale <- as.numeric(scale)

ord <- order(scale)

quartz(width = 14, height = 7)
layout(matrix(1:3, nrow = 1))
plot(scale[ord], R2s[ord], type = "b", main = "R2 of univariate regression and intensity: China", ylab = "R2", xlab = "Spatial aggregation", bty = "n", pch = 16)
plot(scale[ord], coef_crowd[ord], type = "b", main = "Coefficient of crowding: China", ylab = "Coefficient of crowding", xlab = "Spatial aggregation", bty = "n", pch = 16, ylim = c(min(coef_crowd[ord], na.rm = TRUE), 0))
abline(h = 0, lty = 3, col = "red")
plot(scale[ord], log10(p_crowd[ord]), type = "b", main = "p value crowding: China", ylab = "p value", xlab = "Spatial aggregation", bty = "n", pch = 16, yaxt = "n")
at <- c(-5, -4, -3, -2, -1, 0)
axis(2, at = at, labels = round(10^at, 4))
abline(h = log10(0.05), lty = 3, col = "red")

rownames(rhos) <- scale
rhos <- rhos[ord,]
