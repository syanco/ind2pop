###############################################
####                                       ####
####    Individual to Population Niches    ####
####         Home Range Estimation         ####
####          Scott Yanco, PhD             ####
####        scott.yanco@yale.edu           ####
####                                       ####
###############################################

# This script estimates home ranges for each individual which are later used to 
# set the spatial extend for climate vuinerability analyses.

#----   Init   ----#
library(tidyverse)
library(glue)
library(ctmm)
library(raster)

#----   Load Data   ----#
message("Loading data...")
load("out/anno.Rdata")



#----   Perform Analyses   ----#
inds <- unique(gad_anno$individual.local.identifier)

# initialize list to collect results
out <- list()
for(i in 1:length(inds)){
 
  message(glue("Starting ctmm process for Individual {inds[i]}..."))
   #format data as `telemetry` object for ctmm
  evt_tmp <- gad_anno %>% 
    filter(individual.local.identifier == inds[i])
  
  if(nrow(evt_tmp) == 0){
    message("Insufficient records during season, moving to next...")
  } else{
    message(glue("{nrow(evt_tmp)} records found..."))
    evt_telem <- as.telemetry(evt_tmp)
  } 
  
  # get initial acf guess
  guess <- ctmm.guess(evt_telem, interactive = F)
  
  # acf model selection using guess as init
  fit <- ctmm.select(evt_telem, guess)
  
  # fit 
  akde <- akde(evt_telem, fit)
  
  message(glue("Gathering output for Individual {inds[i]}..."))
  # write the akde object to list
  tmp_out <- list()
  tmp_out[["Individual"]] <- as.character(inds[i])
  tmp_out[["aKDE Object"]] <- akde
  tmp_out[["ctmm Fit"]] <- fit
  tmp_out[["PDF Raster"]] <- raster(akde, DF = "PDF")
  tmp_out[["SF Polygon"]] <- as.sf(akde, level.UD=.95)
  tmp_out[["Clipped Raster"]] <- mask(tmp_out[["PDF Raster"]], tmp_out[["SF Polygon"]])
  
  out[[i]] <- tmp_out
  
  message(glue("Individual {inds[i]} complete..."))
}

message("Saving aKDEs to out/akde_out.rdata...")
save(out, file = "out/akde_out.rdata")