###############################################
####                                       ####
####    Individual to Population Niches    ####
####          Scott Yanco, PhD             ####
####        scott.yanco@yale.edu           ####
####                                       ####
###############################################

# This script produces plots and table (except for maps) associated with the 
# analysis. Specifically produces individual contributions to niche breadth 
# (decomposed to variance components), pop-level estimates compared with ind2pop 
# method and lme model, bivariate niche position v breadth plots, and niche 
# vulnerability plots.
# 
# Use "plots" conda env (found in plots.yml)

# Script assumes that data have been previously annotated and empirical analysis 
# completed and saved to rdata file

# TODO: Bivariate plots (niche position X niche breadth - match concpetual fig)
# TODO: Niche vulnerability Plots


#----   Output    ----#

#-- Libraries
library(tidyverse)
# library(rstoat)
library(lubridate)
library(lme4)
library(glue)
library(ggplot2)
library(viridisLite)
library(lme4)
library(RColorBrewer)


#-- Init



#-- Load workspace
load("out/out.Rdata")

# #-- Method Comparison Tables
# 
# #- Gadwall
# df_gad <- data.frame(mu = c(mix_mean_gad, vc_mean_gad, pop_mean_gad),
#                      var = c(mix_var_gad, vc_var_gad, pop_var_gad),
#                      mu_ci_l = c(gad_mean_ci[1], vc_mean_CI_gad[1], pop_mean_gad_CI[1]),
#                      mu_ci_h = c(gad_mean_ci[2], vc_mean_CI_gad[2], pop_mean_gad_CI[2]),
#                      var_ci_l = c(gad_var_ci[1], vc_gad_ci[1], pop_var_gad_CI[1]), 
#                      var_ci_h = c(gad_var_ci[2], vc_gad_ci[2], pop_var_gad_CI[2]))
# row.names(df_gad) <- c("Mixture Dist Est", "RE Var Components", "Pop Est")
# df_gad
# 
# (gad_mu_comp <-df_gad %>% 
#     rownames_to_column(var = "Method") %>% 
#     ggplot() +
#     geom_point(aes(y=mu, x = Method), size = 1.5) +
#     geom_errorbar(aes(x = Method, ymin = mu_ci_l, ymax = mu_ci_h), width = 0.1) +
#     ylim(c(285, 300)) +
#     ylab("Mean") +
#     theme_linedraw(base_size = 12) +
#     coord_flip())
# ggsave(filename = "out/figs/gad_mu_com.png", gad_mu_comp, width = 2.6, height = 1.5, units = "in")
# 
# (gad_var_comp <-df_gad %>% 
#     rownames_to_column(var = "Method") %>% 
#     ggplot() +
#     geom_point(aes(y=var, x = Method), size = 1.5) +
#     geom_errorbar(aes(x = Method, ymin = var_ci_l, ymax = var_ci_h), width = 0.1) +
#     ylim(c(0, 75)) +
#     ylab("Variance") +
#     
#     theme_linedraw(base_size = 12) +
#     theme(axis.title.y=element_blank(),
#           axis.text.y=element_blank(),
#           axis.ticks.y=element_blank())+
#     coord_flip())
# ggsave(filename = "out/figs/gad_var_com.png", gad_var_comp, width = 1.5, height = 1.5, units = "in")

# #-- Individual Contributions
# 
# ind_contributions <- rbind(ind_sum_ele, ind_sum_gad)
# 
# #-- Individual Plots
# 
# #- Total Ind Contribution
# (tot_ele <- ind_sum_ele %>% 
#     # ind_contributions %>% 
#     mutate(ind_f = factor(1:n()),
#            ind_fro = fct_reorder(.f=ind_f, .x=tot_contrib, .fun=min)) %>% 
#     ggplot()+  
#     geom_point(aes(x=tot_contrib, y = ind_fro, color = individual.local.identifier), size = 3)+
#     scale_fill_brewer(palette = "RdYlGn") +
#     # scale_color_viridis_c() +
#     geom_vline(data=NULL, aes(xintercept=0))+
#     theme(axis.title.y=element_blank(),
#           axis.text.y=element_blank(),
#           axis.ticks.y=element_blank(),
#           legend.position = "none") +
#     # xlim(c(minx[s], maxx[s]))+
#     theme_minimal(base_size = 13)+
#     theme(legend.position = "none",
#           axis.text.y = element_blank())+
#     scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+
#     xlab("Total")+
#     ylab(""))
# 
# ggsave(filename = "out/figs/ele_tot_ind.png", tot_ele, width = 2.4, height = 2.4, units = "in")

# #- Ind Contribution via Variance
# var_ele <- ind_sum_ele %>% 
#   # ind_contributions %>% 
#   mutate(ind_f = factor(1:n()),
#          ind_fro = fct_reorder(.f=ind_f, .x=var, .fun=min)) %>% 
#   ggplot()+  
#   geom_point(aes(x=var, y = ind_fro, color = individual.local.identifier), size = 3)+
#   scale_fill_brewer(palette = "RdYlGn") +
#   geom_vline(data=NULL, aes(xintercept=pop_var_ele))+
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         legend.position = "none") +
#   # xlim(c(minx[s], maxx[s]))+
#   theme_minimal()+
#   theme(legend.position = "none",
#         axis.text.y = element_blank())+
#   scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+
#   xlab("Variance")+
#   ylab("")
# 
# ggsave(filename = "out/figs/ele_var_ind.png", var_ele, width = 2.4, height = 2.4, units = "in")

# #- Ind Contribution via Mean
# mu_ele <- ind_sum_ele %>% 
#   # ind_contributions %>% 
#   mutate(ind_f = factor(1:n()),
#          ind_fro = fct_reorder(.f=ind_f, .x=mean_contrib, .fun=min)) %>% 
#   ggplot()+  
#   geom_point(aes(x=mean_contrib, y = ind_fro, color = individual.local.identifier), size = 3)+
#   scale_fill_brewer(palette = "RdYlGn") +
#   geom_vline(data=NULL, aes(xintercept=0))+
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         legend.position = "none") +
#   # xlim(c(minx[s], maxx[s]))+
#   theme_minimal()+
#   theme(legend.position = "none",
#         axis.text.y = element_blank())+
#   scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+
#   xlab("Mean")+
#   ylab("")
# 
# ggsave(filename = "out/figs/ele_mu_ind.png", mu_ele, width = 2.4, height = 2.4, units = "in")

#-- Gadwall

# #- Total Ind Contribution
# (tot_gad <- ind_sum_gad %>% 
#     # ind_contributions %>% 
#     mutate(ind_f = factor(1:n()),
#            ind_fro = fct_reorder(.f=ind_f, .x=tot_contrib, .fun=min)) %>% 
#     ggplot()+  
#     geom_point(aes(x=tot_contrib, y = ind_fro, color = individual.local.identifier), size = 3)+
#     scale_fill_brewer(palette = "RdYlGn") +
#     # scale_color_viridis_c() +
#     geom_vline(data=NULL, aes(xintercept=0))+
#     theme(axis.title.y=element_blank(),
#           axis.text.y=element_blank(),
#           axis.ticks.y=element_blank(),
#           legend.position = "none") +
#     # xlim(c(minx[s], maxx[s]))+
#     theme_minimal(base_size = 13)+
#     theme(legend.position = "none",
#           axis.text.y = element_blank())+
#     scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+
#     xlab("Total")+
#     ylab(""))
# 
# ggsave(filename = "out/figs/gad_tot_ind.png", tot_gad, width = 2.4, height = 2.4, units = "in")

# #- Ind Contribution via Variance
# var_gad <- ind_sum_gad %>% 
#   # ind_contributions %>% 
#   mutate(ind_f = factor(1:n()),
#          ind_fro = fct_reorder(.f=ind_f, .x=var, .fun=min)) %>% 
#   ggplot()+  
#   geom_point(aes(x=var, y = ind_fro, color = individual.local.identifier), size = 3)+
#   scale_fill_brewer(palette = "RdYlGn") +
#   geom_vline(data=NULL, aes(xintercept=pop_var_ele))+
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         legend.position = "none") +
#   # xlim(c(minx[s], maxx[s]))+
#   theme_minimal()+
#   theme(legend.position = "none",
#         axis.text.y = element_blank())+
#   scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+
#   xlab("Variance")+
#   ylab("")
# 
# ggsave(filename = "out/figs/gad_var_ind.png", var_gad, width = 2.4, height = 2.4, units = "in")
# 
# #- Ind Contribution via Mean
# mu_gad <- ind_sum_gad %>% 
#   # ind_contributions %>% 
#   mutate(ind_f = factor(1:n()),
#          ind_fro = fct_reorder(.f=ind_f, .x=mean_contrib, .fun=min)) %>% 
#   ggplot()+  
#   geom_point(aes(x=mean_contrib, y = ind_fro, color = individual.local.identifier), size = 3)+
#   scale_fill_brewer(palette = "RdYlGn") +
#   geom_vline(data=NULL, aes(xintercept=0))+
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         legend.position = "none") +
#   # xlim(c(minx[s], maxx[s]))+
#   theme_minimal()+
#   theme(legend.position = "none",
#         axis.text.y = element_blank())+
#   scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+
#   xlab("Mean")+
#   ylab("")
# 
# ggsave(filename = "out/figs/gad_mu_ind.png", mu_gad, width = 2.4, height = 2.4, units = "in")


# #----  Compare background
# 
# #- Elephants
# 
# elcomb <- elephants_anno %>% 
#   mutate(trt = "mv") %>% 
#   bind_rows(ele_bg %>% mutate(trt = "bg")) %>% 
#   bind_rows(elephants_mol_anno %>% mutate(trt = "mol",
#                                           value = lst_day))
# 
# ggplot(elcomb) + 
#   geom_density(aes(x=value, color = trt))
# 
# ggplot(elcomb)+
#   geom_point(aes(x = date, y = value, color = trt), alpha = 0.2)
# 
# #- Gadwall
# 
# gadcomb <- gad_anno %>% 
#   mutate(trt = "mv") %>% 
#   bind_rows(gad_bg %>% mutate(trt = "bg")) %>% 
#   bind_rows(gadwall_mol_anno %>% mutate(trt = "mol"))
# 
# ggplot(gadcomb) + 
#   geom_density(aes(x=value, color = trt))
# 
# ggplot(gadcomb)+
#   geom_point(aes(x = date, y = value, color = trt), alpha = 0.2)


#---- Density Plots

# ggplot()+
#   geom_density(data = elephants_anno, aes(x = value, color = individual.local.identifier)) +
#   # geom_density(data = elephants_mol_anno, aes(x = value), size = 2) +
#   theme_minimal()+
#   theme(legend.position = "none")

# ggplot()+
#   geom_density(data = gad_anno, aes(x = value, color = individual.local.identifier)) +
#   # geom_density(data = elephants_mol_anno, aes(x = value), size = 2) +
#   theme_minimal()+
#   theme(legend.position = "none")


#---- Bivariate plots ----#


# #-- Elephants
# 
# (ele_bivar <- ggplot(ind_sum_ele) +
#     geom_point(aes(x=mu, y = var, color = individual.local.identifier), size = 1.5)+
#     xlab("Niche Position") + 
#     ylab("Niche Breadth") +
#     theme_linedraw() +
#     theme(legend.position = "none")
# )
# 
# ggsave(filename = "out/figs/ele_bivar.png", ele_bivar, width = 2.4, 
#        height = 2.4, units = "in")

#-- Gadwall


(gad_bivar <- ggplot(gad_fut) +
    geom_point(aes(x=mu, y = var, color = fut_w), size = 1.5)+
    scale_color_viridis_c(direction = -1) +
    xlab("Niche Position") + 
    ylab("Niche Breadth") +
    theme_linedraw() +
    theme(legend.position = "none")
)

ggsave(filename = "out/figs/gad_bivar.png", gad_bivar, width = 2.4, 
       height = 2.4, units = "in")



#---- Climate Vulnerability ----#

# #-- Elephants
# 
# #- Ind vulnerability density plot
# #
# # Create df with sim values for plotting
# dens_sims_ele <- ind_sum_ele %>% 
#   group_by(individual.local.identifier) %>% 
#   group_modify(~data.frame(sims = rnorm(1000000, mean = .$mu, sd = sqrt(.$var)))) %>% 
#   full_join(ele_fut)
# 
# ggplot(dens_sims_ele) +
#   geom_density(aes(x = sims, group = individual.local.identifier, color = fut_w))

#- Future population niche denisty plot


#-- Gadwall

#- Ind vulnerability density plot
#
# Create df with sim values for plotting
dens_sims_gad <- ind_sum_gad %>% 
  group_by(individual.local.identifier) %>% 
  group_modify(~data.frame(sims = rnorm(100000, mean = .$mu, sd = sqrt(.$var)))) %>% 
  full_join(gad_fut)

(vul_plot_gad <- ggplot(dens_sims_gad) +
    geom_density(aes(x = sims, group = individual.local.identifier, color = fut_w)) +
    scale_color_viridis_c(direction = -1) +
    theme_linedraw() +
    theme(legend.position = "none") +
    xlab("Temperature")
)

ggsave(filename = "out/figs/gad_vuln.png", vul_plot_gad, width = 2.4, 
       height = 2.4, units = "in")

#- Future population niche density plot

fut_gad_pop_sim <- data.frame(type = c(rep("Future", 10000), 
                                       rep("Current", 10000)),
                              sims = c(rnorm(10000, mean = mix_mean_gad_fut$pop_mu_fut, 
                                             sd = sqrt(mix_var_gad_fut)),
                                       rnorm(10000, mean = mix_mean_gad, 
                                             sd = sqrt(mix_var_gad)))
)

(fut_gad_plot <- ggplot(fut_gad_pop_sim) +
    geom_density(aes(x = sims, color = type))+
    scale_color_manual(values = c("black", "red")) +
    theme_linedraw() +
    theme(legend.position = "none") +
    xlab("Temperature")
)

ggsave(filename = "out/figs/gad_fut.png", fut_gad_plot, width = 2.4, 
       height = 2.4, units = "in")
