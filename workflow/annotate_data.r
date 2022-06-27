###########################
#                         #
#      Annotate Data      #
#     Scott Yanco, PhD    #
#  scott.yancoe@yale.edu  #
#                         #
###########################

# Annotate datasets using `rstoat`.  Currently uses Gadwall.  Movement data come from movebank.org, 
# static data are from GBIF. Annotated with LST day temp from MODIS.


#-------------------#
#----   Inits   ----#
#-------------------#

# Libraries
library(tidyverse)
library(rstoat)
library(lubridate)
library(sf)
library(glue)

# Custom Functions
# source("src/funs/get_mol_annos.r")
source("src/funs/big_stoat.r")

# Declare annotation variables
vars <- list(list(id = "MODIS/061/MOD11A1",
                  static = FALSE,
                  reducers = list("mean"),
                  s_buff = 1000,
                  t_buff = 1,
                  bands = list("LST_Day_1km")))

# declare start and stop dates for gadwall seasons
gadstart <- 121
gadstop <- 243
gaddur <- 30



#--------------------------------#
#--- Individual Movement Data ---#
#--------------------------------#


#-- Gadwall --#

message("Loading gadwall data...")
gadwall <- read.csv("data/MPIAB Gadwall.csv") %>% 
  rename(lng = location.long,
         lat = location.lat) %>% 
  mutate(date = as_date(timestamp)) %>% 
  filter(!is.na(lng),
         !is.na(lat)) %>% 
  mutate(doy = yday(date)) %>% 
  filter(doy > gadstart & doy < gadstop)

message("Annotating gadwall movements...")
gad_anno <- big_stoat(data = gadwall, vars = vars) %>% 
  mutate(date = date(as.character(date))) # this might not be necessary

# write to file to avoid re-annotating
write.csv(gad_anno, file = "data/gadwall_annotated.csv")
# gad_anno <- read.csv(file = "data/gadwall_annotated.csv")



#---------------------#
#--- GBIF/MOL Data ---#
#---------------------#


#-- Gadwall --#
message("Loading MOL data...")


gad_mol <- st_read("data/gbif_anas2/gbif_anas2.shp") %>%
  mutate(week = week(eventDt),
         doy = yday(eventDt)) %>%
  filter(!duplicated(molid),
         doy > gadstart & doy < gadstop) #looks like the dataset from Yani is duplicated???

# Temporal Filter - ensure occurrences are temporally clustered within weeks

# get 0.75 quantile sample size within weeks
samp_targ <- gad_mol %>%
  group_by(week) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  summarise(targ = quantile(n, 0.75)) %>%
  pull(targ) %>%
  as.numeric()

# This dataset to max  0.75 quantile as the upper boundary for weekly samples
gadwall_mol_filt <- gad_mol %>%
  group_by(week) %>%
  slice(1:samp_targ) %>%
  ungroup()

# Spatial Filter - filter by grid cell to reduce spatial clustering of samples

# Create 10km by 10km grid across bounding box of data
gad_grd <- st_as_sf(gadwall_mol_filt, coords = c("lng", "lat"), crs = 4326) %>%
  st_transform(crs = 6933) %>% #transform to an equal area proj
  st_make_grid(cellsize = c(10000, 10000)) %>% # make grid
  st_sf() %>% # convert geometry sdet to sf
  mutate(ID = 1:nrow(.)) # add identifier column for each cell

# create `sf` version of data for spatial join
gad_mol_sf <- gadwall_mol_filt %>%
  st_transform(crs = 6933) %>% #transform to an equal area proj
  st_join(gad_grd, join = st_within) #join with grid polygons

# Get target max number of smaples per cell (set at 0.75 quantile of sample size)
spat_targ_gad <- count(as_tibble(gad_mol_sf), ID) %>%  # get num. samples in each cell
  summarise(quantile(n, 0.75)) %>% # calc 0.75 quantile
  as.numeric() # convert to numeric

# set 0.75 quantile as the upper boundary for weekly samples w/i ea. cell
gad_mol_filt2 <- gad_mol_sf %>%
  group_by(ID) %>%
  slice(1:spat_targ_gad) %>%
  ungroup() %>%
  st_transform(4326)

# convert back to df and add stoat columns for annotation
gadwall_mol <- gad_mol_filt2 %>%
  mutate(lon = unlist(map(gad_mol_filt2$geometry,1)),
         lat = unlist(map(gad_mol_filt2$geometry,2))) %>%
  st_drop_geometry() %>%
  mutate(lat_stoat = round(lat, digits = 3),
         lon_stoat = round(lon, digits = 3),
         stoatId = glue("{lat_stoat}_{lon_stoat}"))


# Annotate
message("Annotating gadwall from MOL...")

# Shrink dataset for faster computation given spatial buffer
gad_stoat <- gadwall_mol %>%
  select(eventDt, lat_stoat, lon_stoat, stoatId) %>%
  rename(date = eventDt, lat = lat_stoat, lng = lon_stoat) %>%
  distinct(stoatId, date, .keep_all = T)

# join back to original
gadwall_mol_anno <- big_stoat(data = gad_stoat, vars =  vars) %>%
  left_join(gadwall_mol, by = c("stoatId" = "stoatId", "lat" = "lat_stoat", "lng" = "lon_stoat", "date" = "eventDt"))


# write to file to avoid re-annotating
write.csv(gadwall_mol_anno, file = "data/gadwall_mol_annotated.csv")



# #------------------------------#
# #--- Get Future Niche Annos ---#
# #------------------------------#
# 
# 
# #-- Gadwall --#
# 
# # Create new timestamp in movement data
# gadwall_fut <- gadwall %>% 
#   mutate(ts_fut = date,
#          ts_fut = `year<-`(date, 2075),
#          date = ts_fut)
# 
# # Annotate
# message("Annotating gadwall movements... in the *future*...")
# gad_fut_anno <- big_stoat(data = gadwall_fut, vars = vars)
# 
# # Write csv
# write.csv(gad_fut_anno, file = "data/gadwall_fut_annotated.csv")
# # gat_fut_anno <- read.csv("data/gadwall_fut_annotated.csv")



#-------------------------#
#--- Save Rdata Object ---#
#-------------------------#
message("Saving output...")
save(
    gad_anno,
    # gad_fut_anno,
    # gadwall_mol_anno, 
    gadstart,
    gadstop,
    gaddur,
  file = "out/anno.Rdata")


