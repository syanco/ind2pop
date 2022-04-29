###########################
#                         #
#      Annotate Data      #
#     Scott Yanco, PhD    #
#  scott.yancoe@yale.edu  #
#                         #
###########################

# Annotate datasets using `rstoat`.  Currently uses two species as exmaples: 
# African Elephant and Gadwall.  Movement dat come from movebank.org whereas
# static data are pulled directly from mol.org Annotated with max temp from the 
# CMIP5 models (including past climate hindecasts).
# 
# TODO: Add future annotation code below
# TODO: Consider making more flexible - right now this is fully hardcoded for 
#         each sp.



#-------------------#
#----   Inits   ----#
#-------------------#

# Libraries
library(tidyverse)
library(rstoat)
library(lubridate)
library(sf)

# Custom Functions
source("src/funs/get_mol_annos.r")
source("src/funs/big_stoat.r")

# Declare annotation variables
vars <- list(list(id = "NASA/NEX-GDDP",
                  static = FALSE,
                  reducers = list("mean"),
                  s_buff = 27830,
                  t_buff = 1,
                  bands = list("tasmax")))



#--------------------------------#
#--- Individual Movement Data ---#
#--------------------------------#

#-- Elephants --#
# 
# message("Loading elephant data...")
# elephants <- read.csv("data/African elephants in Etosha National Park (data from Tsalyuk et al. 2018).csv") %>%
#   rename(lng = location.long,
#          lat = location.lat) %>%
#   mutate(date = as_date(timestamp))%>%
#   filter(!is.na(lng),
#          !is.na(lat))
# 
# message("Annotating elephant movements...")
# elephants_anno <- big_stoat(data = elephants, vars = vars)

# write to file to avoid re-annotating
# write.csv(elephants_anno, file = "data/elephants_annotated.csv")

#-- Gadwall --#

message("Loading gadwall data...")
gadwall <- read.csv("data/MPIAB Gadwall.csv") %>% 
  rename(lng = location.long,
         lat = location.lat) %>% 
  mutate(date = as_date(timestamp)) %>% 
  filter(!is.na(lng),
         !is.na(lat))

message("Annotating gadwall movements...")
gad_anno <- big_stoat(data = gadwall, vars = vars)

# write to file to avoid re-annotating
# write.csv(gad_anno, file = "data/gadwall_annotated.csv")



#---------------------#
#--- GBIF/MOL Data ---#
#---------------------#


message("Loading MOL data...")

#-- Elephants --#

# # Download daily MODIS LST data for elephants (1km spatial res)
# # Just using the pre-annotations to get the data, we'll annotate to CMIP5 at
# # next step...
# elephant_mol <- get_species_scenario('Loxodonta africana', 'modis', 'lst_day', 
#                                      1000, 1) %>%
#   rename(lng = longitude,
#          lat = latitude,
#          lst_day = value) %>%
#   mutate(date = as_date(eventdate),
#          week = week(date))
# # Annotate
# message("Annotating elephants from MOL...")
# elephants_mol_anno <- big_stoat(data = elephant_mol, vars =  vars)
# 
# # Temporal Filter - ensure occurrences are temporally clustered within weeks 
# 
# # get 0.75 quantile sample size within weeks
# samp_targ <- elephants_mol_anno %>%
#   group_by(week) %>%
#   summarise(n=n()) %>%
#   ungroup() %>%
#   summarise(quantile(n, 0.75)) %>%
#   as.numeric()
# 
# # This dataset to max  0.75 quantile as the upper boundary for weekly samples
# elephants_mol_anno <- elephants_mol_anno %>%
#   group_by(week) %>%
#   slice(1:samp_targ) %>% 
#   ungroup()
# 
# # Spatial Filter - filter by grid cell to reduce spatial clustering of samples
# 
# # Create 10km by 10km grid across bounding box of data
# ele_grd <- st_as_sf(elephants_mol_anno, coords = c("lng", "lat"), crs = 4326) %>% 
#   st_transform(crs = 6933) %>% #transform to an equal area proj
#   st_make_grid(cellsize = c(10000, 10000)) %>% # make grid
#   st_sf() %>% # convert geometry sdet to sf
#   mutate(ID = 1:nrow(.)) # add identifier column for each cell
# 
# # create `sf` version of data for spatial join
# ele_mol_sf <- st_as_sf(elephants_mol_anno, coords = c("lng", "lat"), crs = 4326) %>% 
#   st_transform(crs = 6933) %>% #transform to an equal area proj
#   st_join(ele_grd, join = st_within) #join with grid polygons
# 
# # Get target max number of smaples per cell (set at 0.75 quantile of sample size) 
# spat_targ_ele <- count(as_tibble(ele_mol_sf), ID) %>%  # get num. samples in each cell 
#   summarise(quantile(n, 0.75)) %>% # calc 0.75 quantile
#   as.numeric() # convert to numeric
# 
# # set 0,75 quantile as the upper boundary for weekly samples w/i ea. cell
# ele_mol_sf <- ele_mol_sf %>%
#   group_by(ID) %>%
#   slice(1:spat_targ_ele)
# 
# # convert back to df
# elephant_mol <- ele_mol_sf %>% 
#   as.data.frame()
# 
# # write to file to avoid re-annotating
# write.csv(elephant_mol, file = "data/elephants_mol_annotated.csv")


#-- Gadwall --#

# TODO: move gadwall anotation up here like I did for elephants

# pull pre-annotations from MOL, filter to movement bounding box + 1 degree
gad_mol <- get_species_scenario('Anas strepera', 'modis', 'lst_day', 1000, 1) %>%
  rename(lng = longitude,
         lat = latitude,
         lsat_evi = value) %>%
  mutate(date = as_date(eventdate),
         week = week(date))

# Annotate
message("Annotating elephants from MOL...")
gadwall_mol_anno <- big_stoat(data = gad_mol, vars =  vars)

# Temporal Filter - ensure occurrences are temporally clustered within weeks 

# get 0.75 quantile sample size within weeks
samp_targ <- gadwall_mol_anno %>%
  group_by(week) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  summarise(quantile(n, 0.75)) %>%
  as.numeric()

# This dataset to max  0.75 quantile as the upper boundary for weekly samples
gadwall_mol_anno <- gadwall_mol_anno %>%
  group_by(week) %>%
  slice(1:samp_targ) %>% 
  ungroup()

# Spatial Filter - filter by grid cell to reduce spatial clustering of samples

# Create 10km by 10km grid across bounding box of data
gad_grd <- st_as_sf(gadwall_mol_anno, coords = c("lng", "lat"), crs = 4326) %>% 
  st_transform(crs = 6933) %>% #transform to an equal area proj
  st_make_grid(cellsize = c(10000, 10000)) %>% # make grid
  st_sf() %>% # convert geometry sdet to sf
  mutate(ID = 1:nrow(.)) # add identifier column for each cell

# create `sf` version of data for spatial join
gad_mol_sf <- st_as_sf(gadwall_mol_anno, coords = c("lng", "lat"), crs = 4326) %>% 
  st_transform(crs = 6933) %>% #transform to an equal area proj
  st_join(ele_grd, join = st_within) #join with grid polygons

# Get target max number of smaples per cell (set at 0.75 quantile of sample size) 
spat_targ_gad <- count(as_tibble(gad_mol_sf), ID) %>%  # get num. samples in each cell 
  summarise(quantile(n, 0.75)) %>% # calc 0.75 quantile
  as.numeric() # convert to numeric

# set 0,75 quantile as the upper boundary for weekly samples w/i ea. cell
gad_mol_sf <- gad_mol_sf %>%
  group_by(ID) %>%
  slice(1:spat_targ_gad)

# convert back to df
gadwall_mol <- gad_mol_sf %>% 
  as.data.frame()

# write to file to avoid re-annotating
write.csv(gadwall_mol, file = "data/gadwall_mol_annotated.csv")



#------------------------------#
#--- Get Future Niche Annos ---#
#------------------------------#

#-- Elephants --#

# Create new timestamp in movement data
elephant_fut <- elephants %>% 
  mutate(ts_fut = date,
         ts_fut = `year<-`(date, 2075),
         data = ts_fut)

# Annotate
message("Annotating elephant movements... in the *future*...")
ele_fut_anno <- big_stoat(data = elephant_fut, vars = vars)


# Write csv
write.csv(ele_fut_anno, file = "data/elephant_fut_annotated.csv")


#-- Gadwall --#

# Create new timestamp in movement data
gadwall_fut <- gadwall %>% 
  mutate(ts_fut = date,
         ts_fut = `year<-`(date, 2075),
         date = ts_fut)

# Annotate
message("Annotating gadwall movements... in the *future*...")
gad_fut_anno <- big_stoat(data = gadwall_fut, vars = vars)

# Write csv
write.csv(gad_fut_anno, file = "data/gadwall_fut_annotated.csv")



#-------------------------#
#--- Save Rdata Object ---#
#-------------------------#

save(
  # elephants_mol_anno, 
     gadwall_mol_anno, 
     # elephants_anno, 
     gad_anno, 
     # ele_fut_anno, 
     gad_fut_anno,
     file = "out/anno.Rdata")


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
\
#------------------#
#--- DEPRECATED ---#
#------------------#

# #-- Background Sampling
# #
# # Sample background points within movement data vicinity to compare distribution 
# # of values to annotations.
# 
# #- Elephants
# 
# lat_range_ele <-c(min(elephants$lat), max(elephants$lat))
# lng_range_ele <-c(min(elephants$lng), max(elephants$lng))
# date_range_ele <- c(min(elephants$date), max(elephants$date))
# 
# ele_bg <- data.frame(lat = runif(min=lat_range_ele[1], max=lat_range_ele[2], 
#                                  n=10000),
#                      lng = runif(min=lng_range_ele[1], max=lng_range_ele[2], 
#                                  n=10000),
#                      date = as_date(runif(date_range_ele[1], 
#                                           max=date_range_ele[2], 
#                                           n=10000)))
# 
# message("Annotating elephant background...")
# ele_bg_anno <- big_stoat(data = ele_bg, vars = vars)
# 
# # write to file to avoid re-annotating
# write.csv(ele_bg_anno, file = "data/elephants_bg.csv")
# 
# #- Gadwall
# 
# lat_range_gad <-c(min(gadwall$lat), max(gadwall$lat))
# lng_range_gad <-c(min(gadwall$lng), max(gadwall$lng))
# date_range_gad <- c(min(gadwall$date), max(gadwall$date))
# 
# gad_bg <- data.frame(lat = runif(min=lat_range_gad[1], max=lat_range_gad[2], 
#                                  n=10000),
#                      lng = runif(min=lng_range_gad[1], max=lng_range_gad[2], 
#                                  n=10000),
#                      date = as_date(runif(date_range_gad[1], 
#                                           max=date_range_gad[2], 
#                                           n=10000)))
# 
# message("Annotating elephant background...")
# gad_bg_anno <- big_stoat(data = gad_bg, vars = vars)
# 
# # write to file to avoid re-annotating
# write.csv(gad_bg_anno, file = "data/gadwall_bg.csv")
# 
# 
# elephants_anno <- read_csv("data/elephants_annotated.csv")
# ele_bg <- read_csv("data/elephants_bg.csv")
# ele_mol_anno <- read_csv("data/elephants_mol_annotated.csv")
# 
# gad_anno <- read_csv("data/gadwall_annotated.csv")
# gad_bg <- read_csv("data/gadwall_bg.csv")
# gad_mol_anno <- read_csv("data/gadwall_mol_annotated.csv")

