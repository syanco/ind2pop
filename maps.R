#---- Quick Maps

library(rnaturalearth)
library(ggplot2)
library(raster)
library(glue)
library(tidyverse)
library(sf)
library(ggmap)
library(RColorBrewer)

# Seperate file to register static google maps api key not inlcuded in repo)
source("src/init/google_reg.r")

world <- ne_countries(returnclass = "sf")

##---- Elephants ----##
africa = world[world$continent == "Africa",]

ele <- read.csv("data/elephants_annotated.csv") %>% 
  st_as_sf(coords = c("lng", "lat"), crs = 4326)

cent_sf <- st_as_sfc(st_bbox(ele)) %>% 
  st_centroid() 

cent <- cent_sf %>% 
  st_coordinates()

ele_bg <- get_googlemap(center = cent, zoom = 8, maptype = "terrain",
                        color = "bw")

(ele_inset <- ggplot() +
    geom_sf(data = africa)+
    geom_sf(data = cent_sf, color = "firebrick2", size = 5)+
    coord_sf(crs = 4326)+
    theme_bw(base_size = 12) +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank()))

(ele_main <- ggmap(ele_bg, darken = c(0.6, "white")) +
    geom_sf(data = ele, aes(color = individual.local.identifier), 
            alpha = 0.1, size = 0.2, inherit.aes = F)+
    theme_bw(base_size = 12) +
    coord_sf(crs = 4326,
             # xlim = c(-20, 45), 
             ylim = c(-18, -20)) +
    theme(legend.position = "none",
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())
)

ggsave(filename = "out/figs/elemap_main.png", ele_main)
ggsave(filename = "out/figs/elemap_inset.png", ele_inset)


##---- Gawall ----##

europe <- world[world$continent == "Europe",]

gad <- read.csv("data/gadwall_annotated.csv") %>% 
  st_as_sf(coords = c("lng", "lat"), crs = 4326)

gad_cent_sf <- st_as_sfc(st_bbox(gad)) %>% 
  st_centroid() 

gad_cent <- gad_cent_sf %>% 
  st_coordinates()

gad_bg <- get_googlemap(center = gad_cent, zoom = 8, maptype = "terrain", 
                        color = "bw")

(gad_inset <- ggplot() +
    geom_sf(data = europe)+
    geom_sf(data = gad_cent_sf, color = "firebrick2", size = 5)+
    coord_sf(crs = 4326,
             xlim = c(-20, 45), 
             ylim = c(30, 73))+
    theme_bw(base_size = 12) +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank()))

(gad_main <- ggmap(gad_bg, darken = c(0.6, "white")) +
    geom_sf(data = gad, aes(color = individual.local.identifier), 
            alpha = 0.1, size = 0.2, inherit.aes = F)+
    theme_bw(base_size = 12) +
    coord_sf(crs = 4326,
             xlim = c(10.7, 13),
             ylim = c(48, 49)
             ) +
    theme(legend.position = "none",
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())
)

ggsave(filename = "out/figs/gadmap_main.png", gad_main, height = 2, width = 3)
ggsave(filename = "out/figs/gadmap_inset.png", gad_inset, height = 1.25, width = 1.25)
