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

# Script assumes that data have been previously annotated

# TODO:



#----   Inits   ----#

# Libraries
library(tidyverse)
# library(rstoat)
library(lubridate)
library(lme4)
library(ggplot2)
library(viridisLite)

# Custom Functions

source("src/funs/niche_funs.r")
# source("src/funs/get_mol_annos.r")
# source("src/funs/big_stoat.r")

#------------------------------------------------------------------------------#



#----   Load Data   ----#

#-- Individual Data 

elephants_anno <- read.csv("data/elephants_annotated.csv")

storks_anno <- read.csv("data/storks_annotated.csv")

# # Create bounding box to restrict GBIF records to region of interest
# bb_ele <- data.frame(maxlat = max(elephants$lat),
#                      minlat = min(elephants$lat),
#                      maxlon = max(elephants$lng),
#                      minlon = min(elephants$lng))

#-- GBIF/MOL Data 

elephants_mol_anno <- read.csv("data/cranes_mol_annotated.csv")

stork_mol_anno <- read.csv("data/stork_mol_annotated.csv")

#------------------------------------------------------------------------------#



#----   Pop Level Compare   ----#


#-- Population Level

#- Elephants

(pop_mean_ele <- mean(na.omit(elephants_mol_anno$value)))
(pop_var_ele <- var(na.omit(elephants_mol_anno$value)))

#- Storks

(pop_mean_stork <- mean(na.omit(stork_mol_anno$value)))
(pop_var_stork <- var(na.omit(stork_mol_anno$value)))


#-- Rand Effect Model

#- Elephants

fm_ele <- lmer(value ~ 1 + (1|individual.local.identifier), data = elephants_anno)
vc_ele <- as.data.frame(VarCorr(fm_ele))
vc_var_ele <- sum(vc_ele$vcov)
vc_mean_ele <- fixef(fm_ele)

#- Storks

fm_stork <- lmer(value ~ 1 + (1|individual.local.identifier), data = elephants_anno)
vc_stork <- as.data.frame(VarCorr(fm_ele))
vc_var_stork <- sum(vc_ele$vcov)
vc_mean_stork <- fixef(fm_ele)


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

#- Storks

# Calc individual means and vars
ind_sum_storks <- storks_anno %>%
  # mutate(ind_f = as.factor(individual_id)) %>% 
  group_by(individual.local.identifier) %>% 
  summarise(mu = mean(na.omit(value)),
            var = var(na.omit(value)))

# make population estimates
(mix_mean_storks <-mean(na.omit(ind_sum_storks$mu))) 
(mix_var_cranes <- estPopVar(var= na.omit(ind_sum_storks$var), 
                             means = na.omit(ind_sum_storks$mu),
                             pop_mean = mix_mean_storks))

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
                                  mu_pop = pop_mean_storks, 
                                  n = nrow(ind_sum_storks))
  mean_con_storks[i] <- muContrib(ind_sum_storks$mu[i], 
                                  mu_pop = pop_mean_storks, 
                                  n = nrow(ind_sum_storks))
}
ind_sum_storks$tot_contrib <- contrib_storks
ind_sum_storks$mean_contrib <- mean_con_storks
ind_sum_storks$taxa <- "Ciconia ciconia"


#----   Output    ----#

#-- Method Comparison Tables

#- Elephants
df_ele <- data.frame(mu = c(mix_mean_ele, vc_mean_ele, pop_mean_ele),
                     var = c(mix_var_ele, vc_var_ele, pop_var_ele))
row.names(df_ele) <- c("Mixture Dist Est", "RE Var Components", "Pop Est")
df_ele

#- Storks
df_storks <- data.frame(mu = c(mix_mean_storks, vc_mean_storks, pop_mean_storks),
                        var = c(mix_var_storks, vc_var_storks, pop_var_storks))
row.names(df_storks) <- c("Mixture Dist Est", "RE Var Components", "Pop Est")
df_storks

#-- Individual Contributions

ind_contributions <- rbind(ind_sum_ele, ind_sum_storks)

#-- Individual Plots

#- Total Ind Contribution
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

#- Ind Contribution via Variance
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

#- Ind Contribution via Mean
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
