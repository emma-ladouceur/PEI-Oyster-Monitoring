

# testing out git with students!
library(tidyverse)
library(ggplot2)
library(viridis)
library(lubridate)
library(brms)
library(patchwork)


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

# # lets prep the data for plotting
# mon_date_2013_2024 <- mon_2013_2024 %>% 
#   # parse the date time
#   mutate( parsed_date = parse_datetime(date_collected, format= "%m/%d/%Y")) %>% 
#   # separate the column, don't remove the original
#   separate(parsed_date, c("year", "month", "day"), sep = "-", remove = F) %>%
#   # change year to factor
#   mutate(f_year = as.factor(year)) %>%
#   # unite month and day back together
#   unite("month_day", month:day, remove= F, sep="") %>% mutate( n_month_day = as.numeric(month_day)) %>%
#   # create a numeric year and location clean as a factor
#   mutate(n_year = as.numeric(year)) %>% #mutate(location_clean = as.factor(location_clean)) %>%
#   # mutate( x. = as.Date(parsed_date)) %>% 
#   mutate( julian_date = yday(parsed_date))

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
                                      TRUE ~ location)) %>%
  # my take on numeric area IDs (locations are grouped according to what Bay they flow into)
  mutate( bay = case_when( location_clean == "Foxley River - Gibb's Creek" ~ "Cascumpec",
                           location_clean == "Foxley River - Goff Bridge" ~ "Cascumpec",
                           location_clean == "Foxley River - Lot 6 Pt." ~ "Cascumpec",
                           location_clean == "Kildare River" ~ "Cascumpec",
                           location_clean == "Montrose Bridge" ~ "Cascumpec",
                           location_clean == "Enmore River" ~ "Egmont",
                           location_clean == "Percival River" ~ "Egmont",
                           location_clean == "Bideford River - Green Park" ~ "Malpeque",
                           location_clean == "Bideford River - Old Wharf" ~ "Malpeque",
                           location_clean == "Bideford River - Paugh's Creek" ~ "Malpeque",
                           location_clean == "Bideford River - Station" ~ "Malpeque",
                           location_clean == "Bentick Cove" ~ "Malpeque",
                           location_clean == "Grand River" ~ "Malpeque",
                           location_clean == "Darnley Basin" ~ "Malpeque",
                           location_clean == "Dunk River" ~ "Bedeque",
                           location_clean == "Wilmot River" ~ "Bedeque",
                           location_clean == "Bedeque Bay" ~ "Bedeque",
                           location_clean == "East River - Cranberry Wharf" ~ "Hillsborough",
                           location_clean == "East River - MacWilliams Seafood" ~ "Hillsborough",
                           location_clean == "North River" ~ "Hillsborough",
                           location_clean == "West River" ~ "Hillsborough",
                           location_clean == "Orwell River" ~ "Hillsborough",
                           location_clean == "Pownal Bay" ~ "Hillsborough",
                           location_clean == "Vernon River" ~ "Hillsborough",
                           TRUE ~ location_clean))

# lets have a look at our work
mon_clean %>% select( location_clean, bay, location, area) %>% distinct() %>% arrange( location_clean, area)

# look at the headers (top 6 rows)
head(mon_clean)

# have a look at our timeline
mon_clean %>% select(year) %>% distinct()
mon_clean %>% select(month,day) %>% distinct() %>% arrange(month, day) # june 20 - sept 11
View( mon_clean %>% select(year, month, day) %>% distinct() %>% arrange(year, month, day))
summary(mon_clean)



