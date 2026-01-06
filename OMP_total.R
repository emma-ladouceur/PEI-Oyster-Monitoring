

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
  
  # locations are grouped according to what Bay they flow into
  mutate( bay = case_when( location_clean == "Foxley River - Gibb's Creek" ~ "Cascumpec Bay",
                           location_clean == "Foxley River - Goff Bridge" ~ "Cascumpec Bay",
                           location_clean == "Foxley River - Lot 6 Pt." ~ "Cascumpec Bay",
                           location_clean == "Kildare River" ~ "Cascumpec Bay",
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
                           location_clean == "Oyster Bed Bridge" ~ "Rustico Bay",
                           location_clean == "Brackley Bay" ~ "Brackley Bay",
                           TRUE ~ location_clean)) %>%
  
  mon_clean %>% count(bay)

# zeros are not true zeros, change to NA
mon_clean <- mon_clean %>%  
  mutate(larvae_size = as.character(larvae_size),
         larvae_size = na_if(larvae_size, "0"))

mon_clean <- mon_clean %>%
   mutate(larvae_size = as.factor(larvae_size)) %>%
  #see below line 101, I have got you started here
  mutate( larvae_size_clean = case_when( larvae_size == "80*-150" ~ "80-150",
                                         larvae_size == "100-120/260-280" ~ "100-280",
                                         larvae_size == "90-120, 260-290" ~ "20-290",
                                         larvae_size == "664" ~ "0",
                                         larvae_size == "90-120,270-300" ~ "90-300",
                                         larvae_size == "90-140,280-300" ~ "90-300",
                                         larvae_size == "100-140, 320" ~ "100-320",
                                         larvae_size == "90-130, 250-300" ~ "30-300",
                                         larvae_size == "100-110/280-300" ~ "100-300",
                                         larvae_size == "90-130, 250-300" ~ "90-300",
                                         larvae_size == "90-120, 250-330" ~ "90-330",
                           TRUE ~ larvae_size))

view(mon_clean$larvae_size)
   
# collapsing bays into counties / regions
mon_clean <- mon_clean %>%
mutate( county = case_when( bay %in% c("Cascumpec Bay", "Egmont Bay", "Malpeque Bay", "Bedeque Bay", "New London Bay") ~ "Prince County (West PEI)",
                            bay %in% c("Hillsborough Bay", "Rustico Bay", "Brackley Bay", "Tracadie Bay", "Clarks Bay") ~ "Queens County (Central PEI)",
                            bay %in% c("St. Peters Bay", "Cardigan Bay", "North Lake", "Colville Bay", "Savage Harbour") ~ "Kings County (East PEI)",
                            TRUE ~ "Other"))

# lets have a look at our work
mon_clean %>% select( location_clean, county, bay, location, area) %>% distinct() %>% arrange( location_clean, area)
mon_clean %>% select(larvae_size) %>% distinct()
# need to clean column larval_size and put larval sizes into 3-5 groups
# what will you classes be?
# ask jesse: what are the size groups thta make sense? what is a size group of zero? (or maybe you know)



mon_clean %>% select( location_clean, bay, location, area) %>% distinct() %>% arrange( location_clean, area)

# look at the headers (top 6 rows)
head(mon_clean)

# have a look at our timeline
mon_clean %>% select(year) %>% distinct()
mon_clean %>% select(month,day) %>% distinct() %>% arrange(month, day) # june 20 - sept 11
View( mon_clean %>% select(year, month, day) %>% distinct() %>% arrange(year, month, day))
summary(mon_clean)
head(mon_clean)


larvae_size_mod <- brm( larvae_total ~ water_temp * larvae_size  + (water_temp * larvae_size | bay/location_clean ) + (1 | julian_date/f_year),
                         data = mon_clean , family = lognormal,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99))

save(larvae_size_mod, file = '~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/larvae_size_mod.Rdata')
load("~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/larvae_size_mod.Rdata") 

pp_check(larvae_size_mod)
conditional_effects(larvae_size_mod)

