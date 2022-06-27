#---- Quick Maps

# use 'plots' conda env

library(rnaturalearth)
library(ggplot2)
library(raster)
library(glue)
library(tidyverse)
library(sf)
library(ggmap)
library(RColorBrewer)

# Separate file to register static google maps api key not included in repo)
source("src/init/google_reg.r")

world <- ne_countries(returnclass = "sf")

europe <- world[world$continent == "Europe",]

gad <- read.csv("data/gadwall_annotated.csv") %>% 
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

gad_main <- ggmap(gad_bg, darken = c(0.6, "white")) +
    geom_sf(data = gad, aes(color = individual.local.identifier), 
            alpha = 0.1, size = 0.4, inherit.aes = F)+
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
                               l = 0))

# gad_main + inset_element(gad_inset, 0.6, -0.56, 0.99, 1)

ggsave(filename = "out/figs/gadmap_main.png", gad_main, height = 3, width = 3)
ggsave(filename = "out/figs/gadmap_inset.png", gad_inset, height = 1.25, width = 1.25)

library(patchwork)
library(ggpubr)

layout <- "
AABBBB
CCCDDD
##EE##
"
layout <- c(
  area(t = 2, l = 1, b = 5, r = 4),
  area(t = 1, l = 3, b = 3, r = 5)
)
(out_plot <- 
  ((gad_main | fut_gad_plot) + plot_layout(widths = c(2,2))) / (gad_bivar | gad_dens) / as_ggplot(leg)+
    plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D', '' )))+
    plot_layout(guides = "keep"))
  
(out_plot <- 
    ((((gad_main / gad_bivar) | (fut_gad_plot / gad_dens)) + 
       plot_layout(widths = c(3,5))) / (as_ggplot(leg)))+
    plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D', '' )))+
    plot_layout(guides = "keep",
                heights = c(5,5,1)))

(
  (((gad_main+inset_element(gad_inset, 0.6, -0.56, 0.99, 1)) | fut_gad_plot)+ 
    plot_layout(widths = NA, heights = NA))/ 
(((gad_bivar|gad_dens)+ 
    plot_layout(widths = NA, heights = NA))+plot_layout(guides = 'collect')) + plot_layout(heights = c(20,1))
  ) +
  plot_layout(heights = c(8, 10))


(gad_main+inset_element(gad_inset, 0.6, -0.56, 0.99, 1)
|fut_gad_plot)
