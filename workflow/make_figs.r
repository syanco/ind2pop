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
# Use "plots" conda env

# Script assumes that data have been previously annotated and empirical analysis 
# completed and saved to rdata file



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
library(rnaturalearth)
library(raster)
library(sf)
library(ggmap)
library(patchwork)
library(ggpubr)

#-- Init
fig_w <- 3.2
fig_h <- 3.2

# Separate file to register static google maps api key not included in repo)
source("src/init/google_reg.r")


#-- Load workspace
message("Loading data...")
load("out/out.Rdata")

fut_weights <- readRDS("out/fut-weights.rds")


#---- Density Plot ----#

message("Creating density plot...")

# Join observed data with projected future data
tot2 <- tot %>%   
  inner_join(fut_weights, by = c("ID" = "individual"))

# get canonical list of ind from gad_anno...
inds_keep <- unique(tot2$ID)

gad_anno2 <- gad_anno %>% 
  inner_join(fut_weights, by = c("individual.local.identifier" = "individual")) %>% 
  filter(individual.local.identifier %in% inds_keep)


# Create plot
(gad_dens <- ggplot()+
    geom_density(data = gad_anno2, aes(x = value,
                                       group = individual.local.identifier,
                                       color = fut_weight),
                 alpha = 0.7,
                 bw = 2) +
    # scale_color_gradient2(low = "red", mid = "white", high = "blue")+
    scale_color_viridis_c(direction = 1) +
    xlab("Temperature (c)")+
    theme_linedraw() +
    theme(legend.position = "none"))

# Save plot
ggsave(filename = "out/figs/gad_dens.png", gad_dens, width = fig_w, 
       height = fig_h, units = "in")



#---- Bivariate Plot ----#

message("Creating bivariate plot...")

# dev.new(width=fig_w, height=fig_h)
# Create plot
(gad_bivar <- ggplot(tot2) +
    geom_point(aes(x=marginality_sigma2, y = specialization_sigma2,
                   color = fut_weight), alpha = 1)+
    scale_size_continuous(name = "Skewness") +
    # scale_color_gradient2(low = "red", mid = "white", high = "blue")+
    scale_color_viridis_c(direction = 1, name = "Future Weight") +
    xlab("Contribution from niche position") + 
    ylab("Contributions from niche breadth") +
    theme_linedraw()+
    guides(color = guide_colorbar(title.position = "top",
                                  title.hjust = 0.5))+
    theme(legend.position="bottom", legend.justification = "center",
          legend.text=element_text(size=6),
          legend.title=element_text(size=8))
)

#Extract and save legend
leg <- cowplot::get_legend(gad_bivar)
# ggsave(filename = "out/figs/leg.png", as_ggplot(leg), dpi = 600)
ggsave(filename = "out/figs/leg.png", leg, dpi = 600)


# Make legend-less plot
gad_bivar_nl <- gad_bivar + theme(legend.position = "none")

# save plot
ggsave(filename = "out/figs/gad_bivar.png", gad_bivar_nl, width = fig_w, 
       height = fig_h, units = "in")


#---- Climate Vulnerability ----#

message("Creating future population distribution plot...")

# Create df with sim values for plotting
dens_sims_gad <- ind_sum_gad %>% 
  group_by(individual.local.identifier) %>% 
  group_modify(~data.frame(sims = rnorm(100000, mean = .$mu_i, sd = sqrt(.$var)))) %>% 
  full_join(tot2, by = c("individual.local.identifier" = "ID"))

#TODO: this should probably be in the calc_niches or calc vulnerabilities script...
message("FUTURE POPULATION MEAN AND CIs:")
(mix_mean_gad_fut <- tot2 %>%
    mutate(w_mu = fut_weight*mu_i) %>%
    summarise(w_sum = sum(w_mu, na.rm = T),
              pop_mu_fut = w_sum/sum(fut_weight))
)
(gad_mean_ci <- mean_CIs(ind_sum_gad))

message("FUTURE POPULATION VARIANCE AND CIs:")

(mix_var_gad_fut <- estPopVar(var= na.omit(tot2$var),
                              means = na.omit(tot2$mu_i),
                              pop_mean = mix_mean_gad_fut$pop_mu_fut,
                              w = na.omit(tot2$fut_weight)))

(gad_var_ci <- var_CIs(ind_sum_gad))


#- Future population niche density plot

fut_gad_pop_sim <- data.frame(type = c(rep("Future", 10000), 
                                       rep("Current", 10000)),
                              sims = c(rnorm(10000, mean = mix_mean_gad_fut$pop_mu_fut, 
                                             sd = sqrt(mix_var_gad_fut)),
                                       rnorm(10000, mean = mix_mean_gad, 
                                             sd = sqrt(mix_var_gad)))
)


# dev.new(width=fig_w, height=fig_h)

# Create plot
(fut_gad_plot <- ggplot(fut_gad_pop_sim) +
  stat_density(aes(x = sims, color = type), geom = "line", position = "identity")+
  scale_color_manual(values = c("black", "red"), name = "Pop. Niche") +
  theme_linedraw() +
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_rect(linetype="solid", 
                                         color= "black"),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8)) +
  xlab("Temperature"))
#Extract and save legend
leg2 <- cowplot::get_legend(fut_gad_plot)
# ggsave(filename = "out/figs/leg.png", as_ggplot(leg), dpi = 600)
ggsave(filename = "out/figs/leg_fut.png", leg2, dpi = 600)


# Make legend-less plot
fut_gad_nl <-fut_gad_plot + theme(legend.position = "none")

# Save plot
ggsave(filename = "out/figs/gad_fut.png", fut_gad_nl, width = fig_w, 
       height = fig_h, units = "in")

#---- Quick Maps

message("Creating maps...")

world <- ne_countries(returnclass = "sf")

europe <- world[world$continent == "Europe",]

gad <- read.csv("data/gadwall_annotated.csv") %>% 
  left_join(tot2, by = c("individual.local.identifier"="ID")) %>% 
  st_as_sf(coords = c("lng", "lat"), crs = 4326)

gad_cent_sf <- st_as_sfc(st_bbox(gad)) %>% 
  st_centroid() 

gad_cent <- gad_cent_sf %>% 
  st_coordinates()

gad_bg <- get_googlemap(center = gad_cent, zoom = 8, maptype = "terrain", 
                        color = "bw")

gad_inset <- ggplot() +
  geom_sf(data = europe)+
  geom_sf(data = gad_cent_sf, color = "firebrick2", size = 3)+
  coord_sf(crs = 4326,
           xlim = c(-20, 45), 
           ylim = c(30, 73))+
  theme_bw(base_size = 12) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.ticks.length = unit(0, "mm"),
        panel.background = element_rect(fill='white'),
        plot.background = element_rect(color = "transparent", fill = "transparent"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)
        # panel.spacing = unit(c(0, 0, 0, 0), "null")
  )

(gad_main <- ggmap(gad_bg, darken = c(0.8, "white")) +
    geom_sf(data = gad, aes(color = fut_weight), 
            alpha = 1, size = 1, inherit.aes = F)+
    scale_color_viridis_c(direction = 1) +
    theme_bw(base_size = 12) +
    coord_sf(crs = 4326,
             # xlim = c(10.7, 13),
             # ylim = c(48, 49)
    ) +
    annotation_custom(grob = ggplotGrob(gad_inset),
                      xmin = 13.7, xmax = 14.7,
                      ymin = 47.5, ymax = 48.5)+
    theme(legend.position = "none",
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.border=element_blank(),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0)))

# gad_main + inset_element(gad_inset, 0.6, -0.56, 0.99, 1)

ggsave(filename = "out/figs/gadmap_main.png", gad_main, height = fig_h, width = fig_w)
ggsave(filename = "out/figs/gadmap_inset.png", gad_inset, height = 1.25, width = 1.25)


# 
# layout <- "
# AABBBB
# CCCDDD
# ##EE##
# "
# layout <- c(
#   area(t = 2, l = 1, b = 5, r = 4),
#   area(t = 1, l = 3, b = 3, r = 5)
# )
# (out_plot <- 
#     ((gad_main | fut_gad_plot) + plot_layout(widths = c(2,2))) / (gad_bivar | gad_dens) / as_ggplot(leg)+
#     plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D', '' )))+
#     plot_layout(guides = "keep"))
# 
# (out_plot <- 
#     ((((gad_main / gad_bivar) | (fut_gad_plot / gad_dens)) + 
#         plot_layout(widths = c(3,5))) / (as_ggplot(leg)))+
#     plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D', '' )))+
#     plot_layout(guides = "keep",
#                 heights = c(5,5,1)))
# 
# (
#   (((gad_main+inset_element(gad_inset, 0.6, -0.56, 0.99, 1)) | fut_gad_plot)+ 
#      plot_layout(widths = NA, heights = NA))/ 
#     (((gad_bivar|gad_dens)+ 
#         plot_layout(widths = NA, heights = NA))+plot_layout(guides = 'collect')) + plot_layout(heights = c(20,1))
# ) +
#   plot_layout(heights = c(8, 10))
# 
# 
# (gad_main+inset_element(gad_inset, 0.6, -0.56, 0.99, 1)
#   |fut_gad_plot)
