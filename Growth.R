

library(tidyverse)

setwd("~/Dropbox/_Academic/Teaching/UPEI/Data/Oysters/")
growth <- read.csv("Oyster Growth Temp 2021 - present.csv", header= TRUE)

head(growth)


# get locations
locations <- growth %>% select(latitude, longitude, name) %>% distinct()

# clean it up
area_locations <- growth %>% 
  mutate( area = case_when( name == "Bideford River" ~ 3,
                            name == "Foxley River - Gibb's Creek" ~ 1,
                            name == "Foxley River - Lot 6 Point" ~ 1,
                            name == "Foxley River - Portage" ~ 1,
                            name == "Foxley River - Roxbury" ~ 1,
                            name == "Orwell River" ~ 7,
                            name == "Percival River" ~ 7,
                            name == "Savage Harbour" ~ 8,
                            name == "Souris River" ~ 9,
                            name == "Rustico Bay" ~ 8,
                            )) %>%
 separate( reading_datetime, c( "date", "time" ) , sep = " " , remove = F) %>%
  mutate( parsed_date = parse_datetime(date, format= "%m/%d/%Y")) %>% 
  separate(parsed_date, c("year", "month", "day"), sep = "-", remove = F)

head(area_locations)
