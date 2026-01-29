# how to clean your global environment
rm(list = ls())


# testing out git with students!
library(tidyverse)
library(ggplot2)
library(viridis)
library(lubridate)
library(brms)
library(patchwork)
library(tidybayes)



# set your own personal working directory below
# Emma MBP
setwd("~/Dropbox/_Projects/PEI Oysters/Data/OMP")
# Maddy
setwd("~/Data/OMP/")


# load dataset you want
# mon_2013_2024 <- read.csv("Oyster Monitoring Results 2013- present.csv", header= TRUE)
mon_2025 <- read.csv("Oyster Monitoring Results (10)_2025.csv", header= TRUE)

# head(mon_2013_2024)
head(mon_2025)
view(mon_2025)


# lets prep the data for plotting
mon_date_2025 <- mon_2025 %>% 
  # parse the date time
  mutate( parsed_date = parse_datetime(date_collected, format= "%Y-%m-%d")) %>% 
  # separate the column, don't remove the original
  separate(parsed_date, c("year", "month", "day"), sep = "-", remove = F) %>%
  # change year to factor
  mutate(f_year = as.factor(year)) %>%
  # unite month and day back together
  unite("month_day", month:day, remove= F, sep="") %>% mutate( n_month_day = as.numeric(month_day)) %>%
  # create a numeric year and location clean as a factor
  mutate(n_year = as.numeric(year)) %>% #mutate(location_clean = as.factor(location_clean)) %>%
  # mutate( x. = as.Date(parsed_date)) %>% 
  mutate( julian_date = yday(parsed_date))

head(mon_date_2025)
mon_date_2025 %>% select(f_year) %>% distinct()

# have a look at monitoring data locations
# mon_dat <- mon_date_2013_2024 %>% bind_rows(mon_date_2025)

mon_date_2025 %>% select(area, location) %>% distinct() %>% arrange(location, area)

# clean it up
# use case when to change location spelling
mon_clean <- mon_date_2025 %>%
  mutate( location_clean = case_when( location == "Savage Hbr" ~ "Savage Harbour" ,
                                      location == "St. Peter's Bay" ~ "St. Peters Bay",
                                      location == "Pinette" ~ "Pinette River",
                                      location == "River South Lake" ~ "South Lake",
                                      location == "Clark's Bay" ~ "Clarks Bay",
                                      location == "Foxley River - Goff's Bridge" ~ "Foxley River - Goff Bridge",
                                      location == "Foxley - Goff Bridge" ~ "Foxley River - Goff Bridge",
                                      location == "Foxley - Gibb's Creek" ~ "Foxley River - Gibb's Creek",
                                      location == "Foxley - Lot 6 Pt." ~ "Foxley River - Lot 6 Pt.",
                                      location == "Bideford - Green Park" ~ "Bideford River - Green Park",
                                      location == "Bideford - Old Wharf" ~ "Bideford River - Old Wharf",
                                      location == "Bideford - Paugh's Creek" ~ "Bideford River - Paugh's Creek",
                                      location == "Bideford - Station" ~ "Bideford River - Station",
                                      location == "Rustico" ~ "Rustico Bay",
                                      location == "Oyster Bed Bridge" ~ "Rustico Bay",
                                      location == "Kildare River" ~ "Montrose Bridge",
                                      TRUE ~ location)) %>%
  
  # locations are grouped according to what Bay they flow into
  mutate( bay = case_when( location_clean == "Foxley River - Gibb's Creek" ~ "Cascumpec Bay",
                           location_clean == "Foxley River - Goff Bridge" ~ "Cascumpec Bay",
                           location_clean == "Foxley River - Lot 6 Pt." ~ "Cascumpec Bay",
                           #location_clean == "Kildare River" ~ "Cascumpec Bay",
                           location_clean == "Montrose Bridge" ~ "Cascumpec Bay",
                           location_clean == "Mill River" ~ "Cascumpec Bay",
                           location_clean == "Enmore River" ~ "Egmont Bay",
                           location_clean == "Percival River" ~ "Egmont Bay",
                           location_clean == "Bideford River - Green Park" ~ "Malpeque Bay",
                           location_clean == "Bideford River - Old Wharf" ~ "Malpeque Bay",
                           location_clean == "Bideford River - Paugh's Creek" ~ "Malpeque Bay",
                           location_clean == "Bideford River - Station" ~ "Malpeque Bay",
                           location_clean == "Bentick Cove" ~ "Malpeque Bay",
                           location_clean == "Grand River" ~ "Malpeque Bay",
                           location_clean == "Darnley Basin" ~ "Malpeque Bay",
                           location_clean == "Dunk River" ~ "Bedeque Bay",
                           location_clean == "Wilmot River" ~ "Bedeque Bay",
                           location_clean == "Bedeque Bay" ~ "Bedeque Bay",
                           location_clean == "East River - Cranberry Wharf" ~ "Hillsborough Bay",
                           location_clean == "East River - MacWilliams Seafood" ~ "Hillsborough Bay",
                           location_clean == "North River" ~ "Hillsborough Bay",
                           location_clean == "West River" ~ "Hillsborough Bay",
                           location_clean == "Orwell River" ~ "Hillsborough Bay",
                           location_clean == "Pownal Bay" ~ "Hillsborough Bay",
                           location_clean == "Vernon River" ~ "Hillsborough Bay",
                           location_clean == "Rustico Bay" ~ "Rustico Bay",
                           location_clean == "St. Peters Bay" ~ "St. Peters Bay",
                           location_clean == "Pinette River" ~ "Hillsborough Bay",
                           location_clean == "Clarks Bay" ~ "Clarks Bay",
                           location_clean == "Launching Pond" ~ "Cardigan Bay",
                           location_clean == "Little Basin" ~ "Malpeque Bay",
                           location_clean == "New London" ~ "New London Bay",
                           location_clean == "North Lake" ~ "North Lake",
                           location_clean == "South Lake" ~ "Malpeque Bay",
                           location_clean == "Bedeque Bay" ~ "Bedeque Bay",
                           location_clean == "Tracadie Bay" ~ "Tracadie Bay",
                           location_clean == "Souris River" ~ "Colville Bay",
                           location_clean == "Savage Harbour" ~ "Savage Harbour",
                           location_clean == "Dock River" ~ "Savage Harbour",
                           #location_clean == "Oyster Bed Bridge" ~ "Rustico Bay",
                           location_clean == "Brackley Bay" ~ "Brackley Bay",
                           TRUE ~ location_clean)) %>%
  filter(water_temp != 722.00) %>% filter(water_temp != 99.90) %>% filter(water_temp != 247.00) %>%
  filter(water_temp != 0.0) %>% filter(water_temp != 34.0) %>% 
  filter(salinity != 0.0) %>% filter(salinity != 90.0) %>% filter(salinity != 50.0) %>%
  filter(salinity != 70.0) %>% filter(salinity != 95.0) %>% filter(salinity != 40.0) %>%
  filter(salinity != 85.0) 

# any other crazy values?
summary(mon_clean)
mon_clean %>% filter(water_temp == 722.00)
mon_clean %>% filter(water_temp == 99.90)
mon_clean %>% filter(water_temp == 247.00)
mon_clean %>% filter(water_temp == 0.0)
mon_clean %>% filter(water_temp == 34.0)
mon_clean %>% filter(salinity == 0.0)
mon_clean %>% filter(salinity == 90.0)
mon_clean %>% filter(salinity == 50.0)
mon_clean %>% filter(salinity == 70.0)
mon_clean %>% filter(salinity == 95.0)
mon_clean %>% filter(salinity == 40.0)
mon_clean %>% filter(salinity == 85.0)
colnames(mon_clean)


bay_list <- mon_clean %>% select(bay, location_clean) %>% distinct() %>% arrange(bay, location_clean)
write.csv(bay_list, "~/Data/OMP/bay_list.csv", row.names = FALSE)

view(bay_list)

write.csv(mon_clean, "~/Data/OMP/OMP_clean_2025.csv", row.names = FALSE)
write.csv(mon_clean, "~/Dropbox/_Projects/PEI Oysters/Data/OMP/OMP_clean_2025.csv", row.names = FALSE)

# investigate data for quality and strange values
omp_dat <- mon_clean

omp_dat %>% select(bay, location_clean) %>% distinct() %>% arrange(bay, location_clean)  # 38 locations

omp_dat %>% select(bay) %>% distinct() %>% arrange(bay) # 15 bays

omp_dat %>% select(f_year) %>% distinct() %>% arrange(f_year) # 2013-2025
omp_dat %>% select(julian_date) %>% distinct() %>% arrange(julian_date) # julian date 170 (june 18)- 254 (sept 11)
omp_dat %>% select(water_temp) %>% distinct() %>% arrange(water_temp) # 0 - 34 -- units = deg C
omp_dat %>% select(salinity) %>% distinct() %>% arrange(salinity) # 0 - 30.1 -- units = ppt (parts per thousand)

