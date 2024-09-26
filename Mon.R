
library(tidyverse)
library(ggplot2)
library(viridis)

setwd("~/Dropbox/_Academic/Teaching/UPEI/Data/Oysters/")

mon <- read.csv("Oyster Monitoring Results 2013- present.csv", header= TRUE)

head(mon)

# have a look at monitoring data locations
mon %>% select(area, location) %>% distinct() %>% arrange(location, area)

# clean it up
mon_clean <- mon %>%  mutate( location_clean = case_when( location == "Savage Hbr" ~ "Savage Harbour" ,
                                                          location == "St. Peter's Bay" ~ "St. Peters Bay",
                                                          location == "Pinette" ~ "Pinette River",
                                                          TRUE ~ location)) %>%
# normalise all locations ?
mutate( area_clean = case_when( location == "East River - MacWilliams Seafood" ~ 6,
                                 location == "Foxley - Gibb's Creek" ~ 1,
                                location == "Foxley - Goff Bridge" ~ 1,
                                location == "Foxley - Lot 6 Pt." ~ 1,
                                location == "Grand River" ~ 4,
                                location == "Rustico" ~ 8,
                                # we made decisions
                                location == "Dunk River" ~ 9,
                                location == "Wilmot River" ~ 9,
                                location_clean == "Pinette River" ~ 7,
                                location_clean == "St. Peters Bay" ~ 10,
                                location_clean == "Savage Harbour" ~ 10,
                                TRUE ~ area 

))
# or if more complicated (dependent on two factors?)
# mutate( area_clean = case_when(location == "Savage Harbour" & area == 1 ~ 8, #?
#                                      location == "Rustico" & area == 6 ~ 8, #?
#                                      location == "Rustico" & area == 5 ~ 8, #?
#                                      location == "Grand River" & area == 5 ~ 4, #?
#                                      location == "Foxley - Gibb's Creek" & area == 2 ~ 1,
#                                      location == "Foxley - Gibb's Creek" & area == 3 ~ 1,
#                                      location == "Foxley - Goff Bridge" & area == 3 ~ 1,
#                                      location == "Foxley - Lot 6 Pt." & area == 3 ~ 1,
#                                      location == "Foxley - Lot 6 Pt." & area == 3 ~ 1,
#                                      TRUE ~ location 
#                                      )) %>%


mon_clean %>% select(area, area_clean, location, location_clean) %>% distinct() %>% arrange(area, location)
mon_clean %>% select(area_clean, location_clean) %>% distinct() %>% arrange( location_clean, area_clean)

head(mon_clean)

mon_prep <- mon_clean %>% 
  mutate( parsed_date = parse_datetime(date_collected, format= "%m/%d/%Y")) %>% 
  separate(parsed_date, c("year", "month", "day"), sep = "-", remove = F) %>%
  mutate(f_year = as.factor(year)) %>%
  unite("month_day", month:day, remove= F, sep="") %>% mutate( n_month_day = as.numeric(month_day)) %>%
  mutate(n_year = as.numeric(year)) %>% mutate(location_clean = as.factor(location_clean))


head(mon_prep)

# have a look at our timeline
mon_prep %>% select(year) %>% distinct()
mon_prep %>% select(month,day) %>% distinct() %>% arrange(month, day)
View( mon_prep %>% select(year, month, day) %>% distinct() %>% arrange(year, month, day))
summary(mon_prep)


max_l <- mon_prep %>% group_by(location_clean, f_year) %>% filter( larvae_above_250_microns == max(larvae_above_250_microns))  %>%
  select(location_clean, area_clean, f_year, n_year, month_day, larvae_above_250_microns) 

max_l %>% filter(location_clean == "Dock River")


ggplot(data = max_l #%>% filter(location_clean == "Dock River")
       , aes( y= month_day, x = n_year, 
            group= location_clean, color = as.factor(location_clean) )) + 
  facet_wrap(~area_clean) +
  geom_point( ) +
  geom_line( ) + 
  scale_x_continuous( breaks=c(2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024)) +
  theme_classic()


ggplot(data = max_l #%>% filter(location_clean == "Dock River")
       , aes( y= month_day, x = n_year, group= location_clean, color = as.factor(location_clean)
              )) + 
  facet_wrap(~location_clean) +
  #geom_point( ) +
  geom_smooth(method = lm, se=TRUE) + 
  scale_x_continuous( breaks=c(2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024))+
  theme_classic()



