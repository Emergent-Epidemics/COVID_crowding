#SV Scarpino
#May 2020
#Crowding analysis of Italy (cross validation)

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
library(glmulti)
library(igraph)
library(zoo)

#########
#Globals#
#########
do_plot <- FALSE
run_glmulti <- FALSE
peak_window <- 1 #number of days around the peak for peak intensity
res <- "pop_italy_10k.tif"

###########
#Acc Funcs#
###########
range_cut <- function(x, prob1 = 0.5, prob2 = 1){
  out <- length(which(x >= quantile(x = x, probs = prob1)))/quantile(x = x, probs = prob2)
  return(out)
}

reporting_days <- function(x){
  start <- which(x != 0)[1]
  x.filt <- x[start:length(x)]
  stop <- rev(which(x.filt != 0))[1]
  x.filt.2 <- x.filt[1:stop]
  return(x.filt.2)
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

run_pops <- which(names(pops) %in% res)

for(j in run_pops){
  use.j <- names(pops)[j]
  crowding <- rep(NA, length(prefects))
  total_pops <- rep(NA, length(prefects))
  ent <- rep(NA, length(prefects))
  dists <- rep(NA, length(prefects))
  patchy <- rep(NA, length(prefects))
  contacts <- rep(NA, length(prefects))
  missed <- c()
  pb <- txtProgressBar(1, length(prefects), style=3)
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
    
    #lloyd.i <- lloyd_crowding(qs = pop_points.i[use.pop.i,3])
    lloyd.i <- agg_index(pop_points.i[use.pop.i,3], method = "lloyd", type = "mean-crowding")
    patchy.i <- var(pop_points.i[use.pop.i,3], na.rm = TRUE)/mean(pop_points.i[use.pop.i,3], na.rm = TRUE)
    contacts.i <- sum(pop_points.i[use.pop.i,3] * (pop_points.i[use.pop.i,3]-1), na.rm = TRUE)
    if(length(use.pop.i) == 0){
      crowding[i] <- NA
      total_pops[i] <- NA
      contacts[i] <- NA
      ent[i] <- NA
      patchy[i] <- NA
      missed <- c(missed, as.character(prefects[i]))
    }else{
      crowding[i] <- lloyd.i[[1]]
      patchy[i] <- patchy.i[[1]]
      contacts[i] <- contacts.i
      total_pops[i] <- sum(pop_points.i[use.pop.i,3], na.rm = TRUE)
      if(length(use.pop.i) == 1){
        ent[i] <- NA
      }else{
        ent[i] <- entropy(discretize(pop_points.i[use.pop.i,3], numBins = 5))
      }
    }
    setTxtProgressBar(pb, i)
  }
  results[[use.j]] <- list("crowding" = crowding, "total_pops" = total_pops, "ent" = ent, "dists" = dists, "missed" = missed, "patchy" = patchy, "contacts" = contacts)
  
}


int <- rep(NA, length(prefects))
total <- rep(NA, length(prefects))
peak_cases <- rep(NA, length(prefects))
rep_days <- rep(NA, length(prefects))
var_diff <- rep(NA, length(prefects))
ranges <- rep(NA, length(prefects))
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
  rep_days[i] <- length(rep_days.i)
  new.cases.i <- rep_days.i
  p.i <- new.cases.i/sum(new.cases.i, na.rm = TRUE) #calculate fraction of total cases for each week
  epidemic_intensity.i <- entropy(p.i)^-1 #calculate inverse Shannon entropy
  min_int.i <- entropy(rep(1/length(p.i), length(p.i)))^-1
  epidemic_intensity_norm.i <- epidemic_intensity.i - min_int.i
  if((which.max(new.cases.i)-peak_window) < 0){
    epi_range.i <- (which.max(new.cases.i)-peak_window):(which.max(new.cases.i)+peak_window) + 1 + abs(which.max(new.cases.i)-peak_window)
  }else{
    epi_range.i <- (which.max(new.cases.i)-peak_window):(which.max(new.cases.i)+peak_window)
  }
  peak_cases.i <- sum(new.cases.i[epi_range.i], na.rm = TRUE)
  peak_cases[i] <- peak_cases.i
  int[i] <- epidemic_intensity_norm.i
  total[i] <- sum(new.cases.i, na.rm = TRUE)
  var_diff[i] <- var(diff(new.cases.i), na.rm = TRUE)
  ranges[i] <- range_cut(new.cases.i)
}

data.reg <- data.frame(italy$GID_2[match(prefects, italy$NAME_2)], prefects, results[[res]]$crowding, int, results[[res]]$total_pops, total, peak_cases, peak_cases/total, rep_days, results[[res]]$patchy, results[[res]]$ent, var_diff, ranges, results[[res]]$contacts)
colnames(data.reg) <- c("ID","location", "crowding", "intensity", "pop", "total_cases", "peak_cases", "prop_peak", "reporting_days", "patchiness", "pop_entropy", "var_difference", "range", "contacts")

data.reg$intensity[which(data.reg$intensity > 1/log(2))] <- NA
data.reg$intensity[which(data.reg$intensity == 0)] <- data.reg$intensity[which(data.reg$intensity == 0)] + min(data.reg$intensity[which(data.reg$intensity != 0)], na.rm = TRUE) * 0.98

data.reg$raw_intensity <- data.reg$intensity
data.reg$intensity <- data.reg$intensity/(1/log(2))

#data.reg$patchiness <- (data.reg$patchiness - min(data.reg$patchiness, na.rm = TRUE))
#data.reg$patchiness <- data.reg$patchiness/max(data.reg$patchiness, na.rm = TRUE)
#data.reg$patchiness <- data.reg$patchiness + min(data.reg$patchiness[which(data.reg$patchiness > 0)], na.rm = TRUE)

#data.reg$crowding <- (data.reg$crowding - min(data.reg$crowding, na.rm = TRUE))
#data.reg$crowding <- data.reg$crowding/max(data.reg$crowding, na.rm = TRUE)
#data.reg$crowding <- data.reg$crowding + min(data.reg$crowding[which(data.reg$crowding > 0)], na.rm = TRUE)

mt_reg <- match(data.reg$location, italy_covariates$NAME_2)
data.reg$temp <- italy_covariates$temp[mt_reg]
data.reg$humidity <- italy_covariates$spechum[mt_reg]
data.reg$mobility <- 1/italy_covariates$within_april_mean_over_december[mt_reg]
oos <- rep(NA, nrow(data.reg))
missed_pred <- c()
for(i in 1:nrow(data.reg)){
  in.i <- data.reg[-i,]
  out.i <- data.reg[i,]
  mod.i <- lm(log10(intensity) ~ log10(crowding)+log10(pop)+log10(humidity)+log10(temp)+log10(mobility), data = in.i)
  pred.i <- predict(mod.i, newdata = out.i, type = "response")
  if(length(pred.i) == 0){
    missed_pred <- c(missed_pred, as.character(prefects[-i]))
    next
  }
  oos[i] <- pred.i
}

mod_in <- lm(log10(intensity) ~ log10(pop) + log10(crowding) + log10(mobility) + log10(humidity) + log10(temp), data = data.reg)
R2_in <- round(summary(mod_in)$adj.r.squared, 2)

mod_out <- lm(log10(data.reg$intensity)~oos)
R2_out <- round(summary(mod_out)$adj.r.squared, 2)

layout(matrix(1:2, nrow = 1))
plot(log10(data.reg$intensity)[-mod_in$na.action], mod_in$fitted.values, pch = 16, bty = "n", xlab = "Observed log intensity", ylab = "In-sample estimated log intensity", main = paste0("In Sample R2 ~ ", R2_in))
abline(a = 0, b = 1, lty = 3, lwd = 3, col = "gray")
plot(log10(data.reg$intensity), oos, pch = 16, bty = "n", xlab = "Observed log intensity", ylab = "Out-of-sample estimated log intensity", main = paste0("Out-of-Sample R2 ~ ", R2_out))
abline(a = 0, b = 1, lty = 3, lwd = 3, col = "gray")

spear_italy <- spearman.ci(log10(data.reg$intensity)[-mod_in$na.action], mod_in$fitted.values, nrep = 1000, conf.level = 0.95)

italy_reg_data <- data.reg

italy_reg_data$var_difference <- italy_reg_data$var_difference + min(italy_reg_data$var_difference[which(italy_reg_data$var_difference > 0)], na.rm = TRUE)

if(run_glmulti == TRUE){
  glmulti.italy <-
    glmulti(log10(intensity) ~ log10(patchiness) + log10(pop) + log10(mobility) + log10(contacts), data = italy_reg_data,
            level = 1,               # 2 = pairwise interaction considered
            method = "h",            # h = Exhaustive approach
            crit = "bic",            # bic = BIC as criteria
            confsetsize = 1,        # Keep # best models
            plotty = F, report = F,  # Plot or interim reports
            fitfunction = lm)
  best.mod.italy <- slot(object = glmulti.italy, name = "objects")[[1]]
}

write.csv(italy_reg_data, file = "../Data/italy_reg_data.csv", row.names = FALSE, quote = FALSE)
