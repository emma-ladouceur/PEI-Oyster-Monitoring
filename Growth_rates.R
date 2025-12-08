
library(readxl)
library(tidyverse)


# Emma MBP
setwd("~/Dropbox/_Projects/PEI Oysters/Data/Growth")
# Maddy
setwd("~/Data/OMP/")

growth_rates <- read_excel("Growth_Rates_2013_2025.xlsx", sheet = 3)

head(growth_rates)
colnames(growth_rates)

# get locations, just have a look
gather_gr <- growth_rates %>% select(-c(`...1`, `...43`, `...44`)) %>% 
  slice(-1) %>%
  mutate(across(2:41, as.numeric)) %>%
pivot_longer(
  cols = 2:41,
  names_to = "variable",
  values_to = "value"
)   

head(gather_gr)
gather_gr %>% select(variable) %>% distinct()


av_temp <- gather_gr %>% filter(str_detect(variable, "Temp")) %>%
  mutate(av_temp = value) %>%
  mutate(
    location = case_when(
      variable == "Average Temp...6" ~ "Portage",
      variable == "Average Temp...10" ~ "Lot 6 Pt",
      variable == "Average Temp...14" ~ "Dump Rd",
      variable == "Average Temp...18" ~ "AnotherSite",
      variable == "Average Temp...22" ~ "AnotherSite",
      variable == "Average Temp...26" ~ "AnotherSite",
      variable == "Average Temp...30" ~ "AnotherSite",
      variable == "Average Temp...34" ~ "AnotherSite",
      variable == "Average Temp...38" ~ "AnotherSite",
      variable == "Average Temp...42" ~ "AnotherSite",
      TRUE ~ "Unknown"
    ),
    location = factor(location, levels = c("Portage", "AnotherSite", "Unknown"))
  )

head(av_temp)
av_temp %>% select(variable) %>% distinct()

gr_dat <- gather_gr %>% filter(!str_detect(variable, "Temp")) %>%
  mutate(
    stage = case_when(
      str_detect(variable, "seed") ~ "Seed",
      str_detect(variable, "1yr")  ~ "1-year",
      str_detect(variable, "2yr")  ~ "2-year",
      TRUE ~ NA_character_
    ),
    stage = factor(stage, levels = c("Seed", "1-year", "2-year")),
    
    # remove those substrings from `variable`
    variable = str_remove_all(variable, " seed| 1yr| 2yr")
  ) %>% mutate(location = variable, growth_rate= value) %>% select(-c(variable, value)) 

head(gr_dat)
gr_dat %>% select(location) %>% distinct()


# clean it up
area_locations <- growth %>% 
  # use case when to change
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
