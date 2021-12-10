######################
#    Annotate HPC    #
######################

# Quick script to annotate larger datasets on HPC (rather than local) using `rstoat`


#----   Inits   ----#

# Libraries
library(tidyverse)
library(rstoat)
library(lubridate)

# Custom Functions
source("src/funs/get_mol_annos.r")
source("src/funs/big_stoat.r")


#----   Load Data   ----#

#-- Individual Data 

message("Loading movement data...")

# load Kruger Elephant data from file, modify vars to work with STOAT
elephants <- read.csv("data/African elephants in Etosha National Park (data from Tsalyuk et al. 2018).csv") %>% 
  rename(lng = location.long,
         lat = location.lat) %>% 
  mutate(date = as_date(timestamp))%>% 
  filter(!is.na(lng),
         !is.na(lat))

storks <- read.csv("data/MPIO white stork lifetime tracking data (2013-2014)-gps.csv") %>% 
  rename(lng = location.long,
         lat = location.lat) %>% 
  mutate(date = as_date(timestamp)) %>% 
  filter(!is.na(lng),
         !is.na(lat))


#-- GBIF/MOL Data 

message("Loading MOL data...")
# pull pre-annotations from MOL, filter to movement bounding box + 1 degree
elephant_mol <- get_species_scenario('Loxodonta africana', 'landsat8', 'evi', 100, 16) %>% 
  # filter(latitude >= (bb_ele$minlat-1) & latitude <= (bb_ele$maxlat+1),
  # longitude >= (bb_ele$minlon-1) & longitude <= (bb_ele$maxlon+1)) %>% 
  rename(lng = longitude,
         lat = latitude,
         lsat_evi = value) %>% 
  mutate(date = as_date(eventdate),
         week = week(date)) 

# make sure occurrences are not too temporally clustered within weeks
# get 0.75 quantile sample size within weeks
samp_targ <- elephant_mol %>% 
  group_by(week) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  summarise(quantile(n, 0.75)) %>% 
  as.numeric()

# set 0,75 quantile as the upper boundary for weekly samples
elephant_mol <- elephant_mol %>% 
  group_by(week) %>% 
  slice(1:samp_targ)

# pull pre-annotations from MOL, filter to movement bounding box + 1 degree
stork_mol <- get_species_scenario('Ciconia ciconia', 'landsat8', 'evi', 100, 16) %>% 
  # filter(latitude >= (bb_ele$minlat-1) & latitude <= (bb_ele$maxlat+1),
  # longitude >= (bb_ele$minlon-1) & longitude <= (bb_ele$maxlon+1)) %>% 
  rename(lng = longitude,
         lat = latitude,
         lsat_evi = value) %>% 
  mutate(date = as_date(eventdate),
         week = week(date)) 

# make sure occurrences are not too temporally clustered within weeks
# get 0.75 quantile sample size within weeks
samp_targ <- stork_mol %>% 
  group_by(week) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  summarise(quantile(n, 0.75)) %>% 
  as.numeric()

# set 0,75 quantile as the upper boundary for weekly samples
stork_mol <- stork_mol %>% 
  group_by(week) %>% 
  slice(1:samp_targ)


#----   Annotate Data   ----#

#declare varibales in Google Earth Engine format
vars <- list(list(id = "MODIS/MCD43A4_006_NDVI",
                  static = FALSE,
                  reducers = list("mean"),
                  s_buff = 1000,
                  t_buff = 1,
                  bands = list("NDVI")))

#-- Movement Data

#- Elephants
message("Annotating elephant movements...")
elephants_anno <- big_stoat(data = elephants, vars = vars)

# write to file to avoid re-annotating
write.csv(elephants_anno, file = "data/elephants_annotated.csv")


# Storks
message("Annotating stork movements...")
storks_anno <- big_stoat(data = storks, vars = vars)

# write to file to avoid re-annotating
write.csv(storks_anno, file = "data/storks_annotated.csv")


#-- MOL data

#- Elephants
message("Annotating elephants from MOL...")
elephants_mol_anno <- big_stoat(data = elephant_mol, vars =  vars)

# write to file to avoid re-annotating
write.csv(elephants_mol_anno, file = "data/elephants_mol_annotated.csv")

#- Storks
message("Annotating storks from MOL...")
stork_mol_anno <- big_stoat(data = stork_mol, vars =  vars)

# write to file to avoid re-annotating
write.csv(stork_mol_anno, file = "data/stork_mol_annotated.csv")

