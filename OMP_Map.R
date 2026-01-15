# how to clean your global environment
rm(list = ls())

# load packages
library(googlesheets4)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(sf)                 
library(rnaturalearth)
library(rnaturalearthdata)
install.packages("remotes")
remotes::install_github("ropensci/rnaturalearthhires")
library(rnaturalearthhires)
library(ggspatial)
library(viridis)
library(MetBrewer)
library(hrbrthemes)
library(mapdata)
library(ggrepel)
library(sp) # For converting to decimal degrees
library(maps)
library(mapproj)
library(ggmap)
library(rworldmap)
library(lubridate)


# Emma working directory
setwd("~/Dropbox/_Projects/PEI Oysters/Data/OMP")
# Maddy working directory 
setwd("~/Data/OMP/")


# load dataset you want
mon_2025_co <- read.csv("bay_list_lat_long.csv", header= TRUE)
head(mon_2025_co)


# convert CSV coordinates to spatial points
points_sf <- st_as_sf(
  mon_2025_co,
  coords = c("longitude", "latitude"),
  crs = 4326
)

# load map a Canada and PEI
canada <- ne_states(country = "Canada", returnclass = "sf")
pei <- canada[canada$name == "Prince Edward Island", ]

ggplot(pei) +
  geom_sf(fill = "gray95", color = "gray60") +
  theme_minimal()


# plot PEI map
ggplot(pei) +
  geom_sf(data = pei, fill = "gray95", color = "gray60") +
  geom_sf(data = points_sf, color = "blue", size = 2) +
  theme_minimal() +
  labs(
    title = "OMP Tow Locations across Prince Edward Island",
    subtitle = "CSV latitude and longitude data"
  )














