
library("googlesheets4")
library(tidyverse)
library("ggplot2")
library(patchwork)
library("sf")                 
library("rnaturalearth")
library("rnaturalearthdata")
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

# set your own personal working directory below
# Emma MBP
setwd("~/Dropbox/_Academic/Teaching/UPEI/Data/PEI - Oysters/")
# Leah PC
setwd("C://Users//lrmeister//OneDrive - University of Prince Edward Island//ACC 3100//Oyster Code")
# Camille PC
setwd("T:/Biodiversity Oyster Data")
# Courtney PC
setwd("C://Users//cmcbride13723//Desktop//R Code Oysters//Oysters")
# Romy

# load dataset you want
mon <- read.csv("Oyster Monitoring Results 2013- present.csv", header= TRUE)
latlon_temp <- read.csv("Oyster Growth Temp 2021 - present.csv", header= TRUE)

head(mon)
head(latlon_temp)

latlon_temp %>% select(name) %>% distinct() %>%  arrange(name)


# have a look at monitoring data locations
mon %>% select(area, location) %>% distinct() %>% arrange(location, area)

# clean it up
# use case when to change location spelling
mon_clean <- mon %>%  mutate( location_clean = case_when( location == "Savage Hbr" ~ "Savage Harbour" ,
                                                          location == "St. Peter's Bay" ~ "St. Peters Bay",
                                                          location == "Pinette" ~ "Pinette River",
                                                          TRUE ~ location)) %>%
  # normalise areas of each location, also using case when
  mutate( area_clean = case_when( location == "East River - MacWilliams Seafood" ~ 6,
                                  location == "Foxley - Gibb's Creek" ~ 1,
                                  location == "Foxley - Goff Bridge" ~ 1,
                                  location == "Foxley - Lot 6 Pt." ~ 1, 
                                  location == "Grand River" ~ 4,
                                  location == "Rustico" ~ 8,
                                  # we made decisions, created new areas for this section, and these locations
                                  location == "Dunk River" ~ 9,
                                  location == "Wilmot River" ~ 9,
                                  location_clean == "Pinette River" ~ 7,
                                  location_clean == "St. Peters Bay" ~ 10,
                                  location_clean == "Savage Harbour" ~ 10,
                                  TRUE ~ area 
                                  
  ))

# lets have a look at our work
mon_clean %>% select(area, area_clean, location, location_clean) %>% distinct() %>% arrange(area, location)
mon_clean %>% select(area_clean, location_clean) %>% distinct() %>% arrange( area_clean)

# look at the headers (top 6 rows)
head(mon_clean)

temp_clean <- latlon_temp %>% mutate( location_clean = name) 


head(temp_clean)

mon_prep <- mon_clean %>% full_join(temp_clean)
View(mon_prep)

mon_prep %>% select(location_clean) %>% distinct() %>% arrange(location_clean)






# theme
theme_set(theme_bw())
# map
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)          
head(world)
ggplot(data = world) +
  geom_sf()
# Get the world polygon
world <- map_data("world")
class(world)
head(world)

pei <- world %>% filter(region == "Canada") %>% filter(subregion == "Prince Edward Island")

pei %>% select(subregion) %>% distinct() %>% arrange(subregion)
head(pei)


# coord_equal version
pei_map <- temp_clean %>%
  ggplot() +
  geom_polygon(data = pei, aes(x=long, y = lat, group = group), fill="grey", alpha=0.7) +
  geom_point(aes(x=longitude, y=latitude,
                 shape= location_clean,
                 color=`location_clean`
  ), size=3, #alpha=0.5
  ) +
  scale_color_viridis(discrete =T )+
   coord_equal() +
  theme_void(base_size=18) +
  theme(
    legend.position = 'none',
  ) +
  #  ggplot2::annotate("text", x = -190, y = -44, hjust = 0, size = 4, label = paste("Study Locations"), color = "black", alpha = 0.5) +
  labs(color= "Realm_Biome")+
  scale_x_continuous(expand = c(0.006, 0.006)) #+ guides(col = guide_legend(ncol = 3))

pei_map


