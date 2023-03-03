###############################################
####                                       ####
####    Individual to Population Niches    ####
####       Calculate Vulnerabilities       ####
####          Scott Yanco, PhD             ####
####        scott.yanco@yale.edu           ####
####                                       ####
###############################################

# This script estimates individual climate change vulnerabilities


#----   Init   ----#
library(tidyverse)
library(glue)
library(ctmm)
library(raster)
library(assertthat)
library(terra)
library(sf)
library(glue)

source("src/funs/niche_funs.r")

warm <- 4

#----   Load Data   ----#
message("Loading data...")

# Movement Data
load("out/anno.Rdata")
# convert K to C
gad_anno <- gad_anno %>% 
  mutate(value = value*.02-273.15)

# aKDEs
load("out/akde_out.rdata")
akdes <- out

# Environmental Raster
lst <- raster("out/modis_lst_mean.tif") 
values(lst) <- values(lst)*.02-273.15 # convert to C

#----   Perform Analyses   ----#
inds <- unique(gad_anno$individual.local.identifier)

# initialize list to collect results
out <- list()
#i <- 1
for(i in 1:length(inds)){
  
  message(glue("Calculating vulnerability for Individual {inds[i]}..."))
  
  dat_ind <- gad_anno %>% 
    filter(individual.local.identifier == inds[i])
  
  akde <- akdes[[i]]
  hr_poly <- akde$`SF Polygon`[glue("{inds[i]} 95% est"),] %>% 
    st_transform(crs = projection(lst))
  
  # Check that we got the right ind...
  assert_that(akde[["Individual"]] == inds[i])
  
  # Make LST raster for individual home range
  ind_lst <- lst %>% terra::crop(hr_poly)
  vals <- values(ind_lst)
  
  cur <- sum(sapply(vals, getUDVal, dat = dat_ind$value), na.rm = T)
  fut <- sum(sapply(values(ind_lst)+warm, getUDVal, dat = dat_ind$value), na.rm = T)
  
  w <- fut/cur
  
  out[[i]] <- data.frame(individual = inds[i],
                         fut_weight = w)
  
}  

df_out <- do.call("rbind", out)

saveRDS(df_out, "out/fut-weights.rds")
