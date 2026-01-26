

#just water and salinity over tiiiiiiime



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


omp_dat <- read.csv("OMP_clean_2025.csv", header= TRUE)
head(omp_dat)

# take the max larvae above 250 microns for every location and year
mon_prep <- omp_dat %>% group_by(location_clean, f_year) %>% # group by location and year
  #filter( larvae_above_250_microns == max(larvae_above_250_microns))  %>% # take max of larvae for each group specified above
  select(location_clean, area, bay, f_year, n_year, month_day, larvae_above_250_microns, parsed_date, julian_date, water_temp, salinity) %>% # select some columns
  # create numeric month day
  mutate(n_month_day = as.numeric(month_day)) %>%
  mutate(bay = as.factor(bay)) %>% 
  mutate(n_month_day = as.numeric(n_month_day)) %>% filter(larvae_above_250_microns > 0) %>%
  group_by(bay) %>%
  #center continuous variables for modelling
  mutate(water_temp.m = water_temp - mean(water_temp, na.rm = TRUE),
         salinity.m = salinity - mean(salinity, na.rm = TRUE),
         n_year.m = n_year - mean(n_year, na.rm = TRUE)
  ) %>% filter(bay != "Tracadie Bay") 

head(mon_prep)

temp_time <- brm( water_temp ~  n_year.m + (  n_year.m  | bay/location_clean ) ,
                         data = mon_prep, iter = 10000, warmup = 1000, family= student(),
                         #control = list(max_treedepth = 15)
                         control = list(adapt_delta = 0.99, max_treedepth = 15)
)

#Emma's paths
save(temp_time, file = "~/Dropbox/_Projects/PEI Oysters/Model_fits/Enviro/temp_time.Rdata")
load("~/Dropbox/_Projects/PEI Oysters/Model_fits/Enviro/temp_time.Rdata")


conditional_effects(temp_time)

sal_time <- brm( salinity ~  n_year.m + (  n_year.m  | bay/location_clean ) ,
                  data = mon_prep, iter = 10000, warmup = 1000, family= student(),
                  #control = list(max_treedepth = 15)
                  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

#Emma's paths
save(sal_time, file = "~/Dropbox/_Projects/PEI Oysters/Model_fits/Enviro/sal_time.Rdata")
load("~/Dropbox/_Projects/PEI Oysters/Model_fits/Enviro/sal_time.Rdata")


conditional_effects(sal_time)


