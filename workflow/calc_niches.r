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
library(moments)

# Custom Functions

source("src/funs/niche_funs.r")

#------------------------------------------------------------------------------#



#----   Load Data   ----#
message("Loading data...")
load("out/anno.Rdata")

# convert K to C
gad_anno <- gad_anno %>% 
  mutate(value = value*.02-273.15)

# gad_fut_anno <- gad_fut_anno %>% 
#   mutate(value = value-273.15)

#-- Individual Data 


# Calc individual means and vars
message("Calculating individual summaries...")
ind_sum_gad <- gad_anno %>%
  mutate(ts = date(as.character(timestamp)),
         doy = yday(ts)) %>%
  filter(doy > gadstart & doy < gadstop) %>% # this should be redundant 
  # mutate(ind_f = as.factor(individual_id)) %>% 
  group_by(individual.local.identifier) %>% 
  summarise(mu = mean(na.omit(value)),
            med = median(na.omit(value)),
            sigma = sd(na.omit(value)),
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

skews <- gad_anno %>% 
  group_by(individual.local.identifier) %>% 
  summarize(sk = skewness(value, na.rm = T))

tot <- individual_contribution(x=ind_sum_gad) %>% 
  left_join(ind_sum_gad, by = c("ID" = "individual.local.identifier")) %>% 
  left_join(skews, by = c("ID" = "individual.local.identifier"))



# make population estimates
message("Calculating population/mixture distribution metrics...")
(mix_mean_gad <-mean(na.omit(ind_sum_gad$mu))) 
(mix_var_gad <- estPopVar(var= na.omit(ind_sum_gad$var), 
                          means = na.omit(ind_sum_gad$mu),
                          pop_mean = mix_mean_gad))
(gad_mean_ci <- mean_CIs(ind_sum_gad))
(gad_var_ci <- var_CIs(ind_sum_gad))


#------------------------------------------------------------------------------#



# #----   Individual Contributions    ----#
# 
# #- Gadwall
# 
# contrib_gad <- c()
# mean_con_gad <- c()
# for(i in 1:nrow(ind_sum_gad)){
#   contrib_gad[i] <- indContrib(ind_sum_gad$mu[i], 
#                                ind_sum_gad$var[i], 
#                                mu_pop = mix_mean_gad, 
#                                n = nrow(ind_sum_gad))
#   mean_con_gad[i] <- muContrib(ind_sum_gad$mu[i], 
#                                mu_pop = mix_mean_gad, 
#                                n = nrow(ind_sum_gad))
# }
# ind_sum_gad$tot_contrib <- contrib_gad
# ind_sum_gad$mean_contrib <- mean_con_gad
# ind_sum_gad$taxa <- "Anas strepera"


####----  Climate Vulnerability  ----####
message("Calculating future vulnerability...")

# declare ind vector for loop
ind_vec <- unique(gad_anno$individual.local.identifier)

# declare future warming offset
warm <- 4

out <- list()

for(i in 1:length(ind_vec)){
  
  #filter to focal ind
  dat <- gad_anno %>% 
    filter(individual.local.identifier == ind_vec[i])
  
  # get minimum obs temp
  mtmp <- min(dat$value, na.rm = T)
  
  # get new min temp
  new_min <- mtmp + warm
  
  # Get empirical kernel density estimate
  d <- density(na.omit(dat$value),
               from = min(na.omit(gad_anno$value)),
               to = max(na.omit(gad_anno$value)),
               bw = 2.5)
  
  # Numerical integration to get AUC left of new min
  xx <- d$x  ## 512 evenly spaced points on [min(x) - 3 * d$bw, max(x) + 3 * d$bw]
  dx <- xx[2L] - xx[1L]  ## spacing / bin size
  yy <- d$y  ## 512 density values for `xx`
  C <- sum(yy) * dx  ## sum(yy * dx)
  p.unscaled <- sum(yy[xx < new_min]) * dx
  p.scaled <- p.unscaled / C
  
  gad_fut <- dat %>% 
    summarise(
      fut_mu = mean(value, na.rm = T) + warm,
      fut_med = median(value, na.rm = T) + warm,
      # fut_var = var(value, na.rm = T) + warm,
      # fut_sigma = sd(value, na.rm = T) + warm,
      ID = individual.local.identifier[1]) %>% 
    inner_join(tot, by = "ID") %>% 
    mutate(
      # fut_d = dnorm(fut_med, mean = mu, sd = fut_sigma),
      # cur_d = dnorm(med, mean = mu, sd = fut_sigma),
      # cur_d2 = max(d$y), #set the current pdens as the max dens in the kernal smooth
      # #then get the fut p by offseting from the x value corresponding to the current max p dens
      # x_max_p = d$x[d$y == cur_d2],
      # fut_d2 = approx(d$x, d$y, xout = c(x_max_p+warm))$y, 
      # fut_w = fut_d/cur_d,
      # fut_w2 = fut_d2/cur_d2
      vuln = p.scaled,
      fut_w = (1-vuln)/1
    )
  
  out[[i]] <- gad_fut
}

gad_fut <- do.call("rbind", out)


# 
# gad_fut2 <- gad_fut %>% 
#   select(fut_mu, fut_sigma, ID, fut_w) %>% 
#   rename(mu = fut_mu, sigma = fut_sigma, ID = ID, 
#          w = fut_w)

# tot_fut <- individual_contribution(gad_fut2, ID = "ID", w = gad_fut2$w)

(mix_mean_gad_fut <- gad_fut %>% 
    mutate(w_mu = fut_w*mu) %>% 
    summarise(w_sum = sum(w_mu, na.rm = T),
              pop_mu_fut = w_sum/sum(fut_w))
)

(mix_var_gad_fut <- estPopVar(var= na.omit(gad_fut$var), 
                              means = na.omit(gad_fut$mu),
                              pop_mean = mix_mean_gad_fut$pop_mu_fut,
                              w = na.omit(gad_fut$fut_w)))
(gad_mean_ci <- mean_CIs(ind_sum_gad))
(gad_var_ci <- var_CIs(ind_sum_gad))




####----  Save Results to File  ----####
message("Saving output to /out/out.Rdata...")
save.image(glue("out/out.Rdata"))


