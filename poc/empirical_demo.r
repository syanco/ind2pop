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

#------------------------------------------------------------------------------#



#----   Load Data   ----#

load("out/anno.Rdata")

# convert K to C
gad_anno <- gad_anno %>% 
  mutate(value = value-273.15)

gad_fut_anno <- gad_fut_anno %>% 
  mutate(value = value-273.15)

#-- Individual Data 

# elephants_anno <- read.csv("data/elephants_annotated.csv")
# # ele_bg <- read.csv("data/elephants_bg.csv")
# 
# gad_anno <- read.csv("data/gadwall_annotated.csv")
# # gad_bg <- read.csv("data/gadwall_bg.csv")
# 
# # # Create bounding box to restrict GBIF records to region of interest
# # bb_ele <- data.frame(maxlat = max(elephants_anno$lat),
# #                      minlat = min(elephants_anno$lat),
# #                      maxlon = max(elephants_anno$lng),
# #                      minlon = min(elephants_anno$lng))
# # 
# TODO: move this to the annotate script to speed up annos...
# bb_gad <- data.frame(maxlat = max(gad_anno$lat)+0,
#                      minlat = min(gad_anno$lat)-0,
#                      maxlon = max(gad_anno$lng)+0,
#                      minlon = min(gad_anno$lng)-0)
# 
# #-- GBIF/MOL Data 
# 
# elephants_mol_anno <- read.csv("data/elephants_mol_annotated.csv") 
# # %>% 
# #   filter(lng < bb_ele$maxlon & lng > bb_ele$minlon,
# #          lat < bb_ele$maxlat & lat > bb_ele$minlat)
# 
# # stork_mol_anno <- read.csv("data/data/stork_mol_annotated.csv") %>% 
# #   filter(lng < bb_sto$maxlon & lng > bb_sto$minlon,
# #          lat < bb_sto$maxlat & lat > bb_sto$minlat)
# 
# gadwall_mol_anno_box <- gadwall_mol_anno %>%
#   filter(lng < bb_gad$maxlon & lng > bb_gad$minlon,
#          lat < bb_gad$maxlat & lat > bb_gad$minlat,
#          year == 2009 | year == 2010) #match tracking years


#------------------------------------------------------------------------------#



#----   Check for overlapping records   ----#


# gadwall_pop <- anti_join(gadwall_mol_anno_box, gad_anno, by = c("lng", "lat", "date"))

# no records with matching lat and long, i.e. no overlap between MOL and MoveBank

#------------------------------------------------------------------------------#



# #-- Population Level --##
# 
# # N.B. rather than calculate the t-based intervals for sample mean I'm just 
# # using lm to shortcut this... results should be roughly identical
# 
# #- Gadwall
# 
# #
# # Mean + CIs
# gad_mod <- lm(value ~ 1, gadwall_pop)
# (pop_mean_gad <- as.numeric(gad_mod$coefficients))
# (pop_mean_gad_CI <- confint(gad_mod, level=0.95))
# 
# # Var + CIs
# (pop_var_gad <- var(na.omit(gadwall_pop$value)))
# (pop_var_gad_CI <- var_CI_chi(var = pop_var_gad, n = length(na.omit(gadwall_pop$value)), alpha = 0.05))


# #-- Rand Effect Model --##
# 
# #- Gadwall
# 
# # Filter Data, summer only, and min
# 
# #create duration summary
# gad_anno_sum <- gad_anno %>% 
#   # mutate(ts = date(as.character(date))) %>% 
#   group_by(individual.local.identifier) %>% 
#   summarise(mints = min(date),
#             maxts = max(date),
#             diffts = difftime(maxts, mints)) %>% 
#   right_join(gad_anno)
# 
# gad_anno2 <- gad_anno_sum %>%
#   filter(diffts > gaddur)
# 
# # fit random effects model to get within/between group variance
# fm_gad <- lmer(value ~ 1 + (1|individual.local.identifier), data = gad_anno2)
# 
# # extarct mean estimate and CIs
# (vc_mean_gad <- fixef(fm_gad))
# (vc_mean_CI_gad <- confint(fm_gad, parm = c("(Intercept)")))
# 
# 
# # Get pooled variance and CIs
# vc_gad <- as.data.frame(VarCorr(fm_gad))
# (vc_var_gad <- sum(vc_gad$vcov))
# (vc_gad_ci <- var_CI_chi(vc_var_gad, n = nrow(fm_gad@frame), 0.05))


#-- New Method --##


#- Gadwall

# Calc individual means and vars
ind_sum_gad <- gad_anno %>%
  mutate(ts = date(as.character(timestamp)),
         doy = yday(ts)) %>%
  filter(doy > 121 & doy < 243) %>%
  # mutate(ind_f = as.factor(individual_id)) %>% 
  group_by(individual.local.identifier) %>% 
  summarise(mu = mean(na.omit(value)),
            var = var(na.omit(value)),
            n = n(),
            m_lat = mean(lat, na.rm = T),
            m_lng = mean(lng, na.rm = T),
            mints = min(ts),
            maxts = max(ts),
            diffts = difftime(maxts, mints),
            mindoy = min(doy),
            maxdoy = max(doy)) %>% 
filter(diffts > 30)

# make population estimates
(mix_mean_gad <-mean(na.omit(ind_sum_gad$mu))) 
(mix_var_gad <- estPopVar(var= na.omit(ind_sum_gad$var), 
                          means = na.omit(ind_sum_gad$mu),
                          pop_mean = mix_mean_gad))
(gad_mean_ci <- mean_CIs(ind_sum_gad))
(gad_var_ci <- var_CIs(ind_sum_gad))


#------------------------------------------------------------------------------#



#----   Individual Contributions    ----#

#- Gadwall

contrib_gad <- c()
mean_con_gad <- c()
for(i in 1:nrow(ind_sum_gad)){
  contrib_gad[i] <- indContrib(ind_sum_gad$mu[i], 
                               ind_sum_gad$var[i], 
                               mu_pop = mix_mean_gad, 
                               n = nrow(ind_sum_gad))
  mean_con_gad[i] <- muContrib(ind_sum_gad$mu[i], 
                               mu_pop = mix_mean_gad, 
                               n = nrow(ind_sum_gad))
}
ind_sum_gad$tot_contrib <- contrib_gad
ind_sum_gad$mean_contrib <- mean_con_gad
ind_sum_gad$taxa <- "Anas strepera"


####----  Climate Vulnerability  ----####

# ind_vec <- unique(gad_anno$individual.local.identifier)

##-- Gadwall

gad_fut <- gad_fut_anno %>% 
  group_by(individual.local.identifier) %>% 
  summarise(fut_mu = mean(value),
            fut_var = var(value)) %>% 
  inner_join(ind_sum_gad) %>% 
  mutate(fut_d = dnorm(fut_mu, mean = mu, sd = sqrt(var)),
         cur_d = dnorm(mu, mean = mu, sd = sqrt(var)),
         fut_w = fut_d/cur_d
  )

(mix_mean_gad_fut <- gad_fut %>% 
    mutate(w_mu = fut_w*mu) %>% 
    summarise(w_sum = sum(w_mu, na.rm = T),
              pop_mu_fut = w_sum/sum(fut_w))
)
(mix_var_gad_fut <- estPopVar(var= na.omit(gad_fut$var), 
                              means = na.omit(gad_fut$mu),
                              pop_mean = mix_mean_gad_fut$pop_mu_fut))
(gad_mean_ci <- mean_CIs(ind_sum_gad))
(gad_var_ci <- var_CIs(ind_sum_gad))




####----  Save Results to File  ----####

save.image(glue("out/out.Rdata"))


