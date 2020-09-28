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
library(glmulti)

#########
#Globals#
#########
do_plot <- FALSE
run_glmulti <- FALSE
peak_window <- 1 #number of days around the peak for peak intensity
res <- "pop_china_10k.tif"

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
china <- readOGR("../Data/China_shapefile/chn_admbnda_adm2_ocha.shp")
china_cases <- readRDS("../Data/chinese_prefactures_cases_5_7_20.RDS")
china_cases$date <- as.POSIXct(china_cases$date)
china_cases$OBJECTID <- as.character(china_cases$OBJECTID) 
data <- readRDS("../Data/admin2_crowding_intensity_4_11.RDS")
curves <- readRDS("../Data/epi_curves_guangzhou_quzhou.RDS")
mobility <- read.csv("../Data/Google COVID-19 Aggregated Mobility Research Dataset.csv")
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

run_pops <- which(names(pops) %in% res)

for(j in run_pops){
  use.j <- names(pops)[j]
  crowding <- rep(NA, length(prefects))
  total_pops <- rep(NA, length(prefects))
  contacts <- rep(NA, length(prefects))
  ent <- rep(NA, length(prefects))
  patchy <- rep(NA, length(prefects))
  missed <- c()
  pb <- txtProgressBar(1, length(prefects), style=3)
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
    
    lloyd.i <- agg_index(pop_points.i[use.pop.i,3], method = "lloyd", type = "mean-crowding")
    contacts.i <- sum(pop_points.i[use.pop.i,3] * (pop_points.i[use.pop.i,3]-1), na.rm = TRUE)
    patchy.i <- var(pop_points.i[use.pop.i,3], na.rm = TRUE)/mean(pop_points.i[use.pop.i,3], na.rm = TRUE)
    if(length(use.pop.i) == 0){
      crowding[i] <- NA
      total_pops[i] <- NA
      ent[i] <- NA
      patchy[i] <- NA
      contacts[i] <- NA
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
  results[[use.j]] <- list("crowding" = crowding, "total_pops" = total_pops, "ent" = ent, "missed" = missed, "patchy" = patchy)
}

int <- rep(NA, length(prefects))
var_diff <- rep(NA, length(prefects))
total <- rep(NA, length(prefects))
peak_cases <- rep(NA, length(prefects))
ranges <- rep(NA, length(prefects))
days <- rep(NA, length(prefects))
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
  days[i] <- length(p.i)
}

mt_pre <- match(data$ADM2_PCODE, prefects)
mt_mob <- match(gsub(pattern = "CN", replacement = "", as.character(data$ADM2_PCODE)), mobility$CODE)

data.reg <- data.frame(prefects[mt_pre], as.character(data$ADM2_EN), results[[res]]$crowding[mt_pre], int[mt_pre], data$temp, total_pops[mt_pre], data$meanhum, as.numeric(as.character(mobility$reduction.in.mobility.during.the.outbreak.period[mt_mob])), total[mt_pre], peak_cases[mt_pre], (peak_cases[mt_pre]/total[mt_pre]), data$repday, results[[res]]$patchy[mt_pre], results[[res]]$ent[mt_pre], var_diff[mt_pre], ranges[mt_pre], contacts[mt_pre])
colnames(data.reg) <- c("ID","location", "crowding", "intensity", "temp", "pop", "humidity", "mobility", "total_cases", "peak_cases", "prop_peak", "reporting_days", "patchiness", "pop_entropy", "var_difference", "range", "contacts")

data.reg$intensity[which(data.reg$intensity > 1/log(2))] <- NA
data.reg$intensity[which(data.reg$intensity == 0)] <- data.reg$intensity[which(data.reg$intensity == 0)] + min(data.reg$intensity[which(data.reg$intensity != 0)], na.rm = TRUE) * 0.98

data.reg$raw_intensity <- data.reg$intensity
data.reg$intensity <- data.reg$intensity/(1/log(2))

oos <- rep(NA, nrow(data.reg))
missed_pred <- c()
for(i in 1:nrow(data.reg)){
  in.i <- data.reg[-i,]
  out.i <- data.reg[i,]
  mod.i <- lm(log10(intensity) ~ log10(crowding)+log10(pop)+log10(temp)+log10(humidity)+log10(mobility), data = in.i)
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
plot(log10(data.reg$intensity)[-mod_in$na.action], mod_in$fitted.values, pch = 16, bty = "n", xlab = "Observed log intensity", ylab = "In-sample estimated log intensity",main = paste0("In Sample R2 ~ ", R2_in))
abline(a = 0, b = 1, lty = 3, lwd = 3, col = "gray")
plot(log10(data.reg$intensity), oos, pch = 16, bty = "n", xlab = "Observed log intensity", ylab = "Out-of-sample estimated log intensity", main = paste0("Out-of-Sample R2 ~ ", R2_out))
abline(a = 0, b = 1, lty = 3, lwd = 3, col = "gray")

spear_china <- spearman.ci(log10(data.reg$intensity)[-mod_in$na.action], mod_in$fitted.values, nrep = 1000, conf.level = 0.95)

china_reg_data <- data.reg

china_reg_data$var_difference <- china_reg_data$var_difference + min(china_reg_data$var_difference[which(china_reg_data$var_difference > 0)], na.rm = TRUE)

if(run_glmulti == TRUE){
  glmulti.china <-
    glmulti(log10(intensity) ~ log10(patchiness) + log10(pop) + log10(crowding) + log10(mobility) + log10(temp) + log10(humidity) + log10(contacts), data = china_reg_data,
            level = 1,               # pairwise interaction considered
            method = "h",            # Exhaustive approach
            crit = "bic",            # BIC as criteria
            confsetsize = 1,        # Keep 10 best models
            plotty = F, report = F,  # Plot or interim reports
            fitfunction = glm)
}

cor.test(log10(china_reg_data$intensity), log10(china_reg_data$total_cases/(china_reg_data$pop)))

china_model <- lm(log10(intensity) ~ log10(crowding) + log10(pop) + log10(humidity) + log10(mobility) + log10(temp), data = china_reg_data)

if(run_glmulti == TRUE){
  glmulti.china <-
    glmulti(log10(intensity) ~ log10(patchiness) + log10(pop) + log10(mobility) + log10(contacts), data = china_reg_data,
            level = 1,               # 2 = pairwise interaction considered
            method = "h",            # h = Exhaustive approach
            crit = "bic",            # bic = BIC as criteria
            confsetsize = 1,        # Keep # best models
            plotty = F, report = F,  # Plot or interim reports
            fitfunction = lm)
  best.mod.china <- slot(object = glmulti.china, name = "objects")[[1]]
}

write.csv(china_reg_data, file = "../Data/china_reg_data.csv", row.names = FALSE, quote = FALSE)
