

# testing out git with students!
library(tidyverse)
library(ggplot2)
library(viridis)

setwd("~/Dropbox/_Academic/Teaching/UPEI/Data/Oysters/")

mon <- read.csv("Oyster Monitoring Results 2013- present.csv", header= TRUE)

head(mon)

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
                                # we made decisions
                                location == "Dunk River" ~ 9,
                                location == "Wilmot River" ~ 9,
                                location_clean == "Pinette River" ~ 7,
                                location_clean == "St. Peters Bay" ~ 10,
                                location_clean == "Savage Harbour" ~ 10,
                                TRUE ~ area 

))

# lets have a look at our work
mon_clean %>% select(area, area_clean, location, location_clean) %>% distinct() %>% arrange(area, location)
mon_clean %>% select(area_clean, location_clean) %>% distinct() %>% arrange( location_clean, area_clean)

# look at the headers (top 6 rows)
head(mon_clean)

# lets prep the data for plotting
mon_prep <- mon_clean %>% 
  # parse the date time
  mutate( parsed_date = parse_datetime(date_collected, format= "%m/%d/%Y")) %>% 
  # separate the column, dont remove the original
  separate(parsed_date, c("year", "month", "day"), sep = "-", remove = F) %>%
  # change year to factor
  mutate(f_year = as.factor(year)) %>%
  # unite month and day back together
  unite("month_day", month:day, remove= F, sep="") %>% mutate( n_month_day = as.numeric(month_day)) %>%
  # create a numeric year and location clean as a factor
  mutate(n_year = as.numeric(year)) %>% mutate(location_clean = as.factor(location_clean))

# have a look at our work
head(mon_prep)

# have a look at our timeline
mon_prep %>% select(year) %>% distinct()
mon_prep %>% select(month,day) %>% distinct() %>% arrange(month, day)
View( mon_prep %>% select(year, month, day) %>% distinct() %>% arrange(year, month, day))
summary(mon_prep)

#
max_l <- mon_prep %>% group_by(location_clean, f_year) %>% # group by location and year
  filter( larvae_above_250_microns == max(larvae_above_250_microns))  %>% # take max of larvae for each group specified above
  select(location_clean, area_clean, f_year, n_year, month_day, larvae_above_250_microns) %>% # select some columns
  # create numeric month day
  mutate(n_month_day = as.numeric(month_day))

# explore our work
head(max_l)
print(max_l %>% ungroup() %>% select(month_day, n_month_day) %>% distinct() %>% arrange(n_month_day), n=50)
max_l %>% filter(location_clean == "Dock River")

# lets use ggplot to make some plots!

ggplot(data = max_l , aes( y= n_month_day, x = n_year,  group= location_clean, color = as.factor(location_clean) )) + 
  facet_wrap(~area_clean) + #facet by each area
  geom_point( ) + # plot raw data points
  geom_line( ) + # use geomline to connect the dots
  # below we manipulate the x and y axis to tell it what to show
  scale_x_continuous( breaks=c(2012, 2014,  2016,  2018,  2020,  2022,  2024))+
  scale_y_continuous( limits = c(700,  825), breaks=c( 710,  725, 810, 825))+
  # use a nice minimalistoverall theme format
  theme_classic()

# this time we use geom_smooth to fit a very simple linear model for every area overall
ggplot(data = max_l  , aes( y= n_month_day, x = n_year, group= area_clean )) + 
  facet_wrap(~area_clean) +
  geom_point(  aes(  group= location_clean, color = as.factor(location_clean) ), alpha =0.5  ) +
  geom_line( aes(  group= location_clean, color = as.factor(location_clean) ), alpha =0.5 ) + 
  geom_smooth(method = lm, se=TRUE, color= "black", alpha =0.5 ) + 
  scale_x_continuous( breaks=c(2012, 2014,  2016,  2018,  2020,  2022,  2024))+
  scale_y_continuous( limits = c(700,  825), breaks=c( 710,  725, 810, 825))+
  theme_classic()

# this time we use grom smooth to fit a line for every area and every location
ggplot(data = max_l , aes( y= n_month_day, x = n_year, group= area_clean)  ) + 
  facet_wrap(~area_clean) +
  geom_point( aes(  group= location_clean, color = as.factor(location_clean) ) , alpha =0.5 ) +
  geom_smooth( method = lm, se=FALSE, aes(  group= location_clean, color = as.factor(location_clean) ), alpha =0.5 ) + 
  geom_smooth(method = lm, se=TRUE, color= "black") + 
  scale_x_continuous( breaks=c(2012, 2014,  2016,  2018,  2020,  2022,  2024))+
  scale_y_continuous( limits = c(700,  825), breaks=c( 710,  725, 810, 825))+
  theme_classic()

