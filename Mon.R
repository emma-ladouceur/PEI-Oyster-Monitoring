

# testing out git with students!
library(tidyverse)
library(ggplot2)
library(viridis)
library(lubridate)
library(brms)
library(patchwork)


# set your own personal working directory below
# Emma MBP
setwd("~/Dropbox/_Academic/Teaching/UPEI/Data/PEI - Oysters/")
# Leah PC
# setwd("C://Users//lrmeister//OneDrive - University of Prince Edward Island//ACC 3100//Oyster Code")
# # Camille PC
# setwd("T:/Biodiversity Oyster Data")
# # Courtney PC
# setwd("C://Users//cmcbride13723//Desktop//R Code Oysters//Oysters")
# Maddy
setwd("~/Desktop/Data/Oysters/")



# load dataset you want
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
                                # we made decisions, created new areas for this section, and these locations
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
  mutate(n_year = as.numeric(year)) %>% mutate(location_clean = as.factor(location_clean)) %>%
 # mutate( x. = as.Date(parsed_date)) %>% 
  mutate( julian_date = yday(parsed_date))

# have a look at our work
head(mon_prep)

# have a look at our timeline
mon_prep %>% select(year) %>% distinct()
mon_prep %>% select(month,day) %>% distinct() %>% arrange(month, day) # june 20 - sept 11
View( mon_prep %>% select(year, month, day) %>% distinct() %>% arrange(year, month, day))
summary(mon_prep)

# take the max larvae above 250 microns for every location and year
max_l <- mon_prep %>% group_by(location_clean, f_year) %>% # group by location and year
  filter( larvae_above_250_microns == max(larvae_above_250_microns))  %>% # take max of larvae for each group specified above
  select(location_clean, area_clean, f_year, n_year, month_day, larvae_above_250_microns, parsed_date, julian_date) %>% # select some columns
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
  scale_y_continuous( limits = c(700,  825), breaks=c( 710,  725,  810, 825))+
  # use a nice minimalist overall theme format
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

# this time we use geom smooth to fit a line for every area and every location
ggplot(data = max_l , aes( y= n_month_day, x = n_year, group= area_clean)  ) + 
  facet_wrap(~area_clean) +
  geom_point( aes(  group= location_clean, color = as.factor(location_clean) ) , alpha =0.5 ) +
  geom_smooth( method = lm, se=FALSE, aes(  group= location_clean, color = as.factor(location_clean) ), alpha =0.5 ) + 
  geom_smooth(method = lm, se=TRUE, color= "black") + 
  scale_x_continuous( breaks=c(2012, 2014,  2016,  2018,  2020,  2022,  2024))+
  scale_y_continuous( limits = c(700,  825), breaks=c( 710,  725, 810, 825))+
  theme_classic()

head(max_l)

# lets use a real model......
library(brms)
summary(max_l)
head(max_l)
max_l_dat <- max_l %>% mutate(area_clean = as.factor(area_clean)) %>% 
  mutate(n_month_day = as.numeric(n_month_day)) %>% filter(larvae_above_250_microns > 0)

head(max_l_dat)
summary(max_l_dat)
print(max_l_dat %>% ungroup() %>% select(n_month_day) %>% distinct() %>% arrange())
max_l_dat %>% filter(larvae_above_250_microns == 0)
print(max_l_dat %>% ungroup() %>% select(julian_date, parsed_date) %>% distinct() %>% arrange(julian_date), n=50)

# we ask, how does timing of max oyster spat vary across time.....(for every area and location)
# oyster_time <- brm( julian_date ~ n_year + (n_year | area_clean/location_clean ), 
#                     data = max_l_dat, iter = 5000, warmup = 1000, control = list(adapt_delta = 0.99))
# 
# 
# save(oyster_time, file = '~/Dropbox/_Academic/Teaching/UPEI/Data/PEI - Oysters/Model_fits/oyster_time.Rdata')
load("~/Dropbox/_Academic/Teaching/UPEI/Data/PEI - Oysters/Model_fits/oyster_time.Rdata") 


summary(oyster_time)

pp_check(oyster_time)
conditional_effects(oyster_time)

oyster_time_fitted <- cbind(oyster_time$data,
                           fitted(oyster_time, re_formula = NA)) %>% 
  as_tibble() %>% 
  # join with plot data for figures
  inner_join(max_l_dat)

head(oyster_time_fitted)

oyster_time_nest <- max_l_dat %>% 
  mutate(area_clean_group = area_clean) %>%
  group_by(area_clean_group, area_clean) %>% 
  summarise(n_year = seq(min(n_year), max(n_year), length.out = 6 )) %>%
  nest(data = c(area_clean, n_year)) %>%
  mutate(predicted = map(data, ~predict(oyster_time, newdata= .x, re_formula = ~( n_year | area_clean) )))

oyster_time.df<- oyster_time_nest %>% unnest(cols = c(data, predicted)) %>% mutate(area_clean = as.numeric(area_clean)) %>% arrange(area_clean)

head(oyster_time.df)

oyster_time_fig <- ggplot() +
  # geom_point(data = max_l_dat,
  #            aes(x = n_year, y = julian_date,
  #                colour = `area_clean`),
  #            size = 1.2, shape=1, position = position_jitter(width = 0.5, height=0.5) ) +
  geom_line(data = oyster_time.df %>% mutate(area_clean =as.factor(`area_clean`)),
            aes(x = n_year, y= (predicted[,1]) ,
                group = `area_clean`,
                colour = `area_clean`),
            size = 0.75) +
  # uncertainy in fixed effect
  geom_ribbon(data = oyster_time_fitted,
              aes(x = n_year, ymin = Q2.5, ymax = Q97.5),
              alpha = 0.3) +
  # fixed effect
  geom_line(data = oyster_time_fitted,
            aes(x = n_year, y = Estimate),
            size = 1.25) +
  scale_y_continuous( limits=c(200,215), breaks = c(204, 208, 212), labels = c("July 23", "July 27", "July 31" )) +
  scale_x_continuous( breaks=c(2012, 2014,  2016,  2018,  2020,  2022,  2024))+
  #scale_color_manual(values = mycolors) +
  scale_color_viridis(discrete = T, option="D")  + 
  #scale_fill_viridis(discrete = T, option="D")  + 
  theme_bw(base_size=18 ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),
                                  legend.direction = "horizontal", legend.position="bottom")  +
   labs(color = "Areas", subtitle = "a)") +
   ylab('Date of max oyster \n larvae > 250 microns') +  xlab("Year of monitoring") 
  # guides(col = guide_legend(ncol = 9))

oyster_time_fig



head(mon_prep)
#next we ask, is max oyster larvae above 250 microns different in different areas?
oyster_max_abund <- brm( larvae_above_250_microns ~ n_year  + (1 +  n_year | area_clean/location_clean ),
                    data = max_l_dat , family = lognormal,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99))

 save(oyster_max_abund, file = '~/Dropbox/_Academic/Teaching/UPEI/Data/PEI - Oysters/Model_fits/oyster_max_abund.Rdata')
load("~/Dropbox/_Academic/Teaching/UPEI/Data/PEI - Oysters/Model_fits/oyster_max_abund.Rdata") 

summary(oyster_area)
summary(oyster_max_abund)

pp_check(oyster_max_abund)
conditional_effects(oyster_area)

oyster_max_abund_fitted <- cbind(oyster_max_abund$data,
                            fitted(oyster_max_abund, re_formula = NA)) %>% 
  as_tibble() %>% 
  # join with plot data for figures
  inner_join(max_l_dat)

head(oyster_max_abund_fitted)

oyster_max_abund_nest <- max_l_dat %>% 
  mutate(area_clean_group = area_clean) %>%
  group_by(area_clean_group, area_clean) %>% 
  summarise(n_year = seq(min(n_year), max(n_year), length.out = 6 )) %>%
  nest(data = c(area_clean, n_year)) %>%
  mutate(predicted = map(data, ~predict(oyster_max_abund, newdata= .x, re_formula = ~(1 + n_year | area_clean) )))

oyster_max_abund.df<- oyster_max_abund_nest %>% unnest(cols = c(data, predicted)) %>% mutate(area_clean = as.numeric(area_clean)) %>% arrange(area_clean)

head(oyster_max_abund.df)

oyster_max_abund_fig <- ggplot() +
  # geom_point(data = max_l_dat,
  #            aes(x = n_year, y = julian_date,
  #                colour = `area_clean`),
  #            size = 1.2, shape=1, position = position_jitter(width = 0.5, height=0.5) ) +
  geom_line(data = oyster_max_abund.df %>% mutate(area_clean =as.factor(`area_clean`)),
            aes(x = n_year, y= (predicted[,1]) ,
                group = `area_clean`,
                colour = `area_clean`),
            size = 0.75) +
  # uncertainy in fixed effect
  geom_ribbon(data = oyster_max_abund_fitted,
              aes(x = n_year, ymin = Q2.5, ymax = Q97.5),
              alpha = 0.3) +
  # fixed effect
  geom_line(data = oyster_max_abund_fitted,
            aes(x = n_year, y = Estimate),
            size = 1.25) +
  # scale_y_continuous( limits=c(200,215), breaks = c(204, 208, 212), labels = c("July 23", "July 27", "July 31" )) +
   scale_x_continuous( breaks=c(2012, 2014,  2016,  2018,  2020,  2022,  2024))+
  #scale_color_manual(values = mycolors) +
  scale_color_viridis(discrete = T, option="D")  + 
  #scale_fill_viridis(discrete = T, option="D")  + 
  theme_bw(base_size=18 ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(colour="black", fill="white"),
                                  legend.direction = "horizontal", legend.position="none")  +
  labs(color = "Areas", subtitle ="b)") +
  ylab('log[ Abundance of max oyster \n larvae > 250 microns ]') +  xlab("Year of monitoring") 
# guides(col = guide_legend(ncol = 9))



g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

area.legend <- g_legend(oyster_time_fig)


( (oyster_time_fig +  theme(legend.position="none")) + (oyster_max_abund_fig) )/(area.legend) + plot_layout(heights = c(10,2)) 



