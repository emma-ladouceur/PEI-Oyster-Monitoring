
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
  rename_with(
    ~ case_when(
      .x == "Average Temp...6"  ~ "Portage Temp",
      .x == "Average Temp...10" ~ "Lot 6 Pt Temp",
      .x == "Average Temp...14" ~ "Dump Rd Temp",
      .x == "Average Temp...18" ~ "Gibb's Creek Temp",
      .x == "Average Temp...22" ~ "Bideford Temp",
      .x == "Average Temp...26" ~ "Percival Temp",
      .x == "Average Temp...30" ~ "Savage Temp",
      .x == "Average Temp...34" ~ "Orwell Temp",
      .x == "Average Temp...38" ~ "Souris Temp",
      .x == "Average Temp...42" ~ "Rustico Temp",
      TRUE ~ .x                 # keep original name if no match
    )
  ) %>% pivot_longer(
  cols = 2:41,
  #cols = -matches("temp"),
  names_to = "variable",
  values_to = "value"
)   %>% mutate(
  stage = case_when(
    str_detect(variable, "Seed") ~ "Seed",
    str_detect(variable, "seed") ~ "Seed",
    str_detect(variable, "1yr")  ~ "1-year",
    str_detect(variable, "2yr")  ~ "2-year",
    str_detect(variable, "Temp")  ~ "Temp",
    TRUE ~ NA_character_
  ),
  stage = factor(stage, levels = c("Seed", "1-year", "2-year", "Temp")),
  # remove those substrings from `variable`
  variable = str_remove_all(variable, " Seed| seed| 1yr| 2yr| Temp|")
) %>% mutate(location = variable, growth_rate = value) %>% select(-c(variable, value)) 



head(gather_gr)
View(gather_gr)
gather_gr %>% select(variable) %>% distinct()



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
