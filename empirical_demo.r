###############################################
####                                       ####
####    Individual to Population Niches    ####
####          Scott Yanco, PhD             ####
####        scott.yanco@yale.edu           ####
####                                       ####
###############################################

# This script demonstrates the basic method of estimating population-level 
# niches from individual-scale data. It also shows how the components of that
# estimation can be decomposed to estimate individual contributions to the 
# population level distribution.  The method relies on simple closed form 
# solutions for mixture distributions. This analysis is being prepared as a part
# of a larger manuscript effort led by Muyang Lu.

# Use "niche_mix" conda env (found in niche_mix.yml)

# Script assumes that data have been previously annotated

# TODO:



#----   Inits   ----#

# Libraries
library(tidyverse)
# library(rstoat)
library(lubridate)
library(lme4)
library(glue)
library(ggplot2)
library(viridisLite)

# Custom Functions

source("src/funs/niche_funs.r")
# source("src/funs/get_mol_annos.r")
# source("src/funs/big_stoat.r")

#------------------------------------------------------------------------------#



#----   Load Data   ----#

#-- Individual Data 

elephants_anno <- read.csv("data/data/elephants_annotated.csv")

# storks_anno <- read.csv("data/data/storks_annotated.csv")

gad_anno <- read.csv("data/gadwall_annotated.csv")

# Create bounding box to restrict GBIF records to region of interest
bb_ele <- data.frame(maxlat = max(elephants_anno$lat),
                     minlat = min(elephants_anno$lat),
                     maxlon = max(elephants_anno$lng),
                     minlon = min(elephants_anno$lng))

# bb_sto <- data.frame(maxlat = max(storks_anno$lat),
#                      minlat = min(storks_anno$lat),
#                      maxlon = max(storks_anno$lng),
#                      minlon = min(storks_anno$lng))

bb_gad <- data.frame(maxlat = max(gad_anno$lat),
                     minlat = min(gad_anno$lat),
                     maxlon = max(gad_anno$lng),
                     minlon = min(gad_anno$lng))

#-- GBIF/MOL Data 

elephants_mol_anno <- read.csv("data/data/elephants_mol_annotated.csv") 
# %>% 
#   filter(lng < bb_ele$maxlon & lng > bb_ele$minlon,
#          lat < bb_ele$maxlat & lat > bb_ele$minlat)

# stork_mol_anno <- read.csv("data/data/stork_mol_annotated.csv") %>% 
#   filter(lng < bb_sto$maxlon & lng > bb_sto$minlon,
#          lat < bb_sto$maxlat & lat > bb_sto$minlat)

gadwall_mol_anno <- read.csv("data/gadwall_mol_annotated.csv") 
# %>%
# filter(lng < bb_gad$maxlon & lng > bb_gad$minlon,
#        lat < bb_gad$maxlat & lat > bb_gad$minlat)

#------------------------------------------------------------------------------#



#----   Check for overlapping records   ----#

# Check for records in common between MOL and MoveBank data

which(elephants_anno$lng %in% elephants_mol_anno$lng) %in% which(elephants_anno$lat %in% elephants_mol_anno$lat)
which(storks_anno$lng %in% stork_mol_anno$lng) %in% which(storks_anno$lat %in% stork_mol_anno$lat)
which(gad_anno$lng %in% gadwall_mol_anno$lng) %in% which(gad_anno$lat %in% gadwall_mol_anno$lat)

# no records with matching lat and long, i.e. no overlap between MOL and MoveBank

#------------------------------------------------------------------------------#




#----   Pop Level Compare   ----#


#-- Population Level

#- Elephants

(pop_mean_ele <- mean(na.omit(elephants_mol_anno$value)))
(pop_var_ele <- var(na.omit(elephants_mol_anno$value)))

# #- Storks
# 
# (pop_mean_stork <- mean(na.omit(stork_mol_anno$value)))
# (pop_var_stork <- var(na.omit(stork_mol_anno$value)))

#- Storks

(pop_mean_gad <- mean(na.omit(gadwall_mol_anno$value)))
(pop_var_gad <- var(na.omit(gadwall_mol_anno$value)))

#-- Rand Effect Model

#- Elephants

fm_ele <- lmer(value ~ 1 + (1|individual.local.identifier), data = elephants_anno)
vc_ele <- as.data.frame(VarCorr(fm_ele))
vc_var_ele <- sum(vc_ele$vcov)
vc_mean_ele <- fixef(fm_ele)

# #- Storks
# 
# fm_stork <- lmer(value ~ 1 + (1|individual.local.identifier), data = storks_anno)
# vc_stork <- as.data.frame(VarCorr(fm_stork))
# vc_var_stork <- sum(vc_stork$vcov)
# vc_mean_stork <- fixef(fm_stork)

#- Gadwall

fm_gad <- lmer(value ~ 1 + (1|individual.local.identifier), data = gad_anno)
vc_gad <- as.data.frame(VarCorr(fm_gad))
vc_var_gad <- sum(vc_gad$vcov)
vc_mean_gad <- fixef(fm_gad)

#-- New Method

#- Elephants

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

# #- Storks
# 
# # Calc individual means and vars
# ind_sum_storks <- storks_anno %>%
#   # mutate(ind_f = as.factor(individual_id)) %>% 
#   group_by(individual.local.identifier) %>% 
#   summarise(mu = mean(na.omit(value)),
#             var = var(na.omit(value)))
# 
# # make population estimates
# (mix_mean_storks <-mean(na.omit(ind_sum_storks$mu))) 
# (mix_var_storks <- estPopVar(var= na.omit(ind_sum_storks$var), 
#                              means = na.omit(ind_sum_storks$mu),
#                              pop_mean = mix_mean_storks))

#- Gadwall

# Calc individual means and vars
ind_sum_gad <- gad_anno %>%
  # mutate(ind_f = as.factor(individual_id)) %>% 
  group_by(individual.local.identifier) %>% 
  summarise(mu = mean(na.omit(value)),
            var = var(na.omit(value)))

# make population estimates
(mix_mean_gad <-mean(na.omit(ind_sum_gad$mu))) 
(mix_var_gad <- estPopVar(var= na.omit(ind_sum_gad$var), 
                          means = na.omit(ind_sum_gad$mu),
                          pop_mean = mix_mean_gad))

#------------------------------------------------------------------------------#



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


#- Storks

contrib_storks <- c()
mean_con_storks <- c()
for(i in 1:nrow(ind_sum_storks)){
  contrib_storks[i] <- indContrib(ind_sum_storks$mu[i], 
                                  ind_sum_storks$var[i], 
                                  mu_pop = pop_mean_stork, 
                                  n = nrow(ind_sum_storks))
  mean_con_storks[i] <- muContrib(ind_sum_storks$mu[i], 
                                  mu_pop = pop_mean_stork, 
                                  n = nrow(ind_sum_storks))
}
ind_sum_storks$tot_contrib <- contrib_storks
ind_sum_storks$mean_contrib <- mean_con_storks
ind_sum_storks$taxa <- "Ciconia ciconia"

#- Gadwall

contrib_gad <- c()
mean_con_gad <- c()
for(i in 1:nrow(ind_sum_gad)){
  contrib_gad[i] <- indContrib(ind_sum_gad$mu[i], 
                               ind_sum_gad$var[i], 
                               mu_pop = pop_mean_gad, 
                               n = nrow(ind_sum_gad))
  mean_con_gad[i] <- muContrib(ind_sum_gad$mu[i], 
                               mu_pop = pop_mean_gad, 
                               n = nrow(ind_sum_gad))
}
ind_sum_gad$tot_contrib <- contrib_gad
ind_sum_gad$mean_contrib <- mean_con_gad
ind_sum_gad$taxa <- "Anas strepera"


save.image(glue("out/out.Rdata"))


