###############################################
####                                       ####
####    Individual to Population Niches    ####
####          Scott Yanco, PhD             ####
####        scott.yanco@yale.edu           ####
####                                       ####
###############################################

# This script demonstrates the basic method of estimating population-level 
# niches from inidvidual-scale data. It also shows how the components of that
# estimation can be decomposed to estimate individual contributions to the 
# population level distribution.  The method relies on simple closed form 
# solutions for mixture distributions. This analysis is being prepared as a part
# of a larger manuscript effort led by Muyang Lu.

# Use "niche_mix" conda env (found in niche_mix.yml)

# TODO:
# + Swap datasets - too few individuals

#----   Inits   ----#

# Libraries
library(tidyverse)
library(rstoat)
library(lubridate)
library(lme4)
library(ggplot2)
library(viridisLite)

# Custom Functions

source("src/funs/niche_funs.r")
source("src/funs/get_mol_annos.r")
source("src/funs/big_stoat.r")

#----   Load Data   ----#

#-- Individual Data 

# load Kruger Elephant data from file, modify vars to work with STOAT
elephants <- read.csv("data/African elephants in Etosha National Park (data from Tsalyuk et al. 2018).csv") %>% 
  rename(lng = location.long,
         lat = location.lat) %>% 
  mutate(date = as_date(timestamp))

# # Create bounding box to restrict GBIF records to region of interest
# bb_ele <- data.frame(maxlat = max(elephants$lat),
#                      minlat = min(elephants$lat),
#                      maxlon = max(elephants$lng),
#                      minlon = min(elephants$lng))

storks <- read.csv("data/MPIO white stork lifetime tracking data (2013-2014)-gps.csv") %>% 
  rename(lng = location.long,
         lat = location.lat) %>% 
  mutate(date = as_date(timestamp))


# load Common Crane data from file, modify vars to work with STOAT
# cranes <- read.csv("data/Common Crane Lithuania GPS, 2015-2016.csv") %>% 
#   rename(lng = location.long,
#          lat = location.lat) %>% 
#   mutate(date = as_date(timestamp))

# Create bounding box to restrict GBIF records to region of interest
# bb_crane <- data.frame(maxlat = max(cranes$lat),
#                        minlat = min(cranes$lat),
#                        maxlon = max(cranes$lng),
#                        minlon = min(cranes$lng))

# tortoises <- read.csv("data/Galapagos Tortoise Movement Ecology Programme_2009-2018.csv") %>% 
#   filter(individual.taxon.canonical.name == "Chelonoides hoodensis")


#-- GBIF/MOL Data 

# scen <- get_scenarios()

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

# cranes_mol <- get_species_scenario('Grus grus', 'landsat8', 'evi', 100, 16) %>% 
#   filter(latitude >= (bb_crane$minlat-1) & latitude <= (bb_crane$maxlat+1),
#          longitude >= (bb_crane$minlon-1) & longitude <= (bb_crane$maxlon+1)) %>% 
#   rename(lng = longitude,
#          lat = latitude,
#          lsat_evi = value) %>% 
#   mutate(date = as_date(eventdate))

# tortoise_mol <- get_species_scenario('Chelonoidis hoodensis', 'landsat8', 'evi', 100, 16)


#----   Annotate Data   ----#

#List stoat native products
# stoat_prod <- get_products()

# declare varibale(s) for annotation
# vars <- "srtm-elevation-1000-1"

#declare varibales in Google Earth Engine format
vars <- list(list(id = "MODIS/MCD43A4_006_NDVI",
                  static = FALSE,
                  reducers = list("mean"),
                  s_buff = 1000,
                  t_buff = 1,
                  bands = list("NDVI")))


#-- Movement Data

# annotate elephants
# elephants_anno <- big_stoat(data = elephants, vars = vars)

# write to file to avoid re-annotating
# write.csv(elephants_anno, file = "data/elephants_annotated.csv")

# load from file
elephants_anno <- read.csv("data/elephants_annotated.csv")

# annotate storks
# storks_anno <- big_stoat(data = storks, vars = vars)

# write to file to avoid re-annotating
# write.csv(storks_anno, file = "data/storks_annotated.csv")

# load from file
storks_anno <- read.csv("data/storks_annotated.csv")

# annotate cranes
# cranes_anno <- big_stoat(data = cranes, vars =  vars)

# write to file to avoid re-annotating
# write.csv(cranes_anno, file = "data/cranes_annotated.csv")

# load from file
# cranes_anno <- read.csv("data/cranes_annotated.csv")


#-- MOL data

# Elephants

# elephants_mol_anno <- big_stoat(data = elephant_mol, vars =  vars)

# write to file to avoid re-annotating
# write.csv(elephants_mol_anno, file = "data/elephants_mol_annotated.csv")

# load from file
elephants_mol_anno <- read.csv("data/cranes_mol_annotated.csv")

# Storks

# stork_mol_anno <- big_stoat(data = stork_mol, vars =  vars)

# write to file to avoid re-annotating
# write.csv(stork_mol_anno, file = "data/stork_mol_annotated.csv")

# load from file
stork_mol_anno <- read.csv("data/stork_mol_annotated.csv")

# Cranes

# cranes_mol_anno <- big_stoat(data = cranes_mol, vars =  vars)

# write to file to avoid re-annotating
# write.csv(cranes_mol_anno, file = "data/cranes_mol_annotated.csv")

# load from file
# cranes_mol_anno <- read.csv("data/cranes_mol_annotated.csv")

#----   Pop Level Compare   ----#

#-- Population Level

#- Elephants

(pop_mean_ele <- mean(na.omit(elephants_mol_anno$value)))
(pop_var_ele <- var(na.omit(elephants_mol_anno$value)))

#- Storks

(pop_mean_stork <- mean(na.omit(stork_mol_anno$value)))
(pop_var_stork <- var(na.omit(stork_mol_anno$value)))

#- Cranes

# (pop_mean_cranes <- mean(na.omit(cranes_mol$value)))
# (pop_var_cranes <- var(na.omit(cranes_mol$value)))

#-- Rand Effect Model

#- Elephants

fm_ele <- lmer(value ~ 1 + (1|individual.local.identifier), data = elephants_anno)
vc_ele <- as.data.frame(VarCorr(fm_ele))
vc_var_ele <- sum(vc_ele$vcov)
vc_mean_ele <- fixef(fm_ele)

#- Storkss

fm_stork <- lmer(value ~ 1 + (1|individual.local.identifier), data = elephants_anno)
vc_stork <- as.data.frame(VarCorr(fm_ele))
vc_var_stork <- sum(vc_ele$vcov)
vc_mean_stork <- fixef(fm_ele)

#- Cranes

fm_crane <- lmer(value ~ 1 + (1|individual.local.identifier), data = cranes_anno)
vc_crane <- as.data.frame(VarCorr(fm_crane))
vc_var_crane <- sum(vc_crane$vcov)
vc_mean_crane <- fixef(fm_crane)


#-- New Method

#- Elephants

#TODO: tons of NAs in the elehpant annotations
# Calc individual means and vars
ind_sum_ele <- elephants_anno %>%
  # mutate(ind_f = as.factor(individual_id)) %>% 
  group_by(individual.local.identifier) %>% 
  summarise(mu = mean(na.omit(value)),
            var = var(na.omit(value)))

# make population estimates
(mix_mean_ele <-mean(na.omit(ind_sum_ele$mu))) 
(mix_var_ele <- estPopVar(var = na.omit(ind_sum_ele$var), 
                          means = na.omit(ind_sum_ele$mu),
                          pop_mean = mix_mean_ele))

#- Cranes

# Calc individual means and vars
ind_sum_cranes <- cranes_anno %>%
  # mutate(ind_f = as.factor(individual_id)) %>% 
  group_by(individual.local.identifier) %>% 
  summarise(mu = mean(na.omit(value)),
            var = var(na.omit(value)))

# make population estimates
(mix_mean_cranes <-mean(na.omit(ind_sum_cranes$mu))) 
(mix_var_cranes <- estPopVar(var= na.omit(ind_sum_cranes$var), 
                             means = na.omit(ind_sum_cranes$mu),
                             pop_mean = mix_mean_cranes))


#----   Individual Contributions    ----#

#- Elephants
contrib_ele <- c()
mean_con_ele <- c()
for(i in 1:nrow(ind_sum_ele)){
  contrib_ele[i] <- indContrib(ind_sum_ele$mu[i], 
                               ind_sum_ele$var[i], 
                               mu_pop = mix_mean_ele, 
                               n = nrow(ind_sum_ele))
  mean_con_ele[i] <- muContrib(ind_sum_ele$mu[i], 
                               mu_pop = mix_mean_ele, 
                               n = nrow(ind_sum_ele))
}
ind_sum_ele$tot_contrib <- contrib_ele
ind_sum_ele$mean_contrib <- mean_con_ele
ind_sum_ele$taxa <- "Loxodonta africana"

# TODO: remove this after test
#check that the variance contribution components
round(((1/nrow(ind_sum_ele))*ind_sum_ele$var)+ind_sum_ele$mean_contrib, digits = 17) == round(ind_sum_ele$tot_contrib, digits = 17)
# they match to within 17 decimals...

#- Cranes
contrib_cranes <- c()
mean_con_cranes <- c()
for(i in 1:nrow(ind_sum_cranes)){
  contrib_cranes[i] <- indContrib(ind_sum_cranes$mu[i], 
                                  ind_sum_cranes$var[i], 
                                  mu_pop = pop_mean_cranes, 
                                  n = nrow(ind_sum_cranes))
  mean_con_cranes[i] <- muContrib(ind_sum_cranes$mu[i], 
                                  mu_pop = pop_mean_cranes, 
                                  n = nrow(ind_sum_cranes))
}
ind_sum_cranes$tot_contrib <- contrib_cranes
ind_sum_cranes$mean_contrib <- mean_con_cranes
ind_sum_cranes$taxa <- "Grus grus"


#----   Output    ----#

#-- Method Comparison Tables

#- Elephants
df_ele <- data.frame(mu = c(mix_mean_ele, vc_mean_ele, pop_mean_ele),
                     var = c(mix_var_ele, vc_var_ele, pop_var_ele))
row.names(df_ele) <- c("Mixture Dist Est", "RE Var Components", "Pop Est")
df_ele

#- Cranes
df_cranes <- data.frame(mu = c(mix_mean_cranes, vc_mean_crane, pop_mean_cranes),
                        var = c(mix_var_cranes, vc_var_crane, pop_var_cranes))
row.names(df_cranes) <- c("Mixture Dist Est", "RE Var Components", "Pop Est")
df_cranes

#-- Individual Contributions

ind_contributions <- rbind(ind_sum_ele, ind_sum_cranes)

#- Individual Plots
ind_sum_ele %>% 
  # ind_contributions %>% 
  mutate(ind_f = factor(1:n()),
         ind_fro = fct_reorder(.f=ind_f, .x=tot_contrib, .fun=min)) %>% 
  ggplot()+  geom_point(aes(x=tot_contrib, y = ind_fro, color = tot_contrib))+
  scale_color_viridis_c() +
  geom_vline(data=NULL, aes(xintercept=0))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none") +
  # xlim(c(minx[s], maxx[s]))+
  xlab("")


ind_sum_ele %>% 
  # ind_contributions %>% 
  mutate(ind_f = factor(1:n()),
         ind_fro = fct_reorder(.f=ind_f, .x=var, .fun=min)) %>% 
  ggplot()+  geom_point(aes(x=var, y = ind_fro, color = var))+
  scale_color_viridis_c() +
  geom_vline(data=NULL, aes(xintercept=0))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none") +
  # xlim(c(minx[s], maxx[s]))+
  xlab("")


ind_sum_ele %>% 
  # ind_contributions %>% 
  mutate(ind_f = factor(1:n()),
         ind_fro = fct_reorder(.f=ind_f, .x=mean_contrib, .fun=min)) %>% 
  ggplot()+  geom_point(aes(x=mean_contrib, y = ind_fro, color = mean_contrib))+
  scale_color_viridis_c() +
  geom_vline(data=NULL, aes(xintercept=0))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none") +
  # xlim(c(minx[s], maxx[s]))+
  xlab("")
