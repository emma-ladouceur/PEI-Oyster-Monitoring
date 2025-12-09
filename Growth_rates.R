
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

av_temp <- gather_gr %>% 
  filter( stage == "Temp") %>%
  select(Year, location, stage, growth_rate) %>%
  rename( average_temp = growth_rate) %>%
  select(-c(stage))%>%
  distinct() %>%
  arrange(Year, location)

head(av_temp)


growth <- gather_gr %>%
  filter( stage != "Temp") %>%
  left_join(av_temp, by= c("Year", "location")) %>%
  rename( reading_datetime = Year,
          name = location,
          growth_rate_mm_per_day = growth_rate) %>%
  select( reading_datetime, stage, name, growth_rate_mm_per_day, average_temp)

View(growth)





