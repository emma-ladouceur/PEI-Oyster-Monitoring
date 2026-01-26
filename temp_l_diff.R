

library(tidyverse)
library(ggplot2)
library(viridis)
library(lubridate)
library(brms)
library(patchwork)
library(tidybayes)
library(scales)
library(tidybayes)
library(emmeans)


# set your own personal working directory below
# Emma MBP
setwd("~/Dropbox/_Projects/PEI Oysters/Data/OMP")
# Maddy
#setwd("~/Data/OMP/")


omp_dat <- read.csv("OMP_clean_2025.csv", header= TRUE)
head(omp_dat)
summary(omp_dat)



l_total_dat <- omp_dat %>%  # create numeric month day
  #group_by(bay, location_clean, f_year) %>% # group by location and year 
 # filter( larvae_total == max(larvae_total))  %>% # take max of larvae for each group specified above
 filter(larvae_total > 0) %>%      
  mutate(n_month_day = as.numeric(month_day)) %>%
  mutate(bay = as.factor(bay)) %>% 
  mutate(n_month_day = as.numeric(n_month_day)) %>% 
  group_by(bay, location_clean) %>%
  #center continuous variables for modelling
  mutate(water_temp.m = water_temp - mean(water_temp, na.rm = TRUE),
         salinity.m = salinity - mean(salinity, na.rm = TRUE),
         n_year.m = n_year - mean(n_year, na.rm = TRUE),
         julian_date.m = julian_date - mean(julian_date, na.rm = TRUE)
  ) %>% # filter(water_temp != 722.00) %>% filter(water_temp != 99.90) %>% #filter(water_temp != 34.00)
  mutate(log_larvae_total = log(larvae_total)) 



temp_time <- brm(
  water_temp ~ n_year.m +
    s(julian_date, k = 10) +
  (1 + n_year.m | bay/location_clean) ,
  family  = student(),
  data= l_total_dat,
 # autocor = brms::cor_ar(~ julian_date | bay/location_clean:year)
)

save(temp_time, file = "~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/temp_time.Rdata")
load("~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/temp_time.Rdata")

pp_check(temp_time, nsamples = 50)
conditional_effects(temp_time)

lt_time <- brm(
  larvae_total ~ n_year.m  +
    s(julian_date, k = 10) +
  (1 + n_year.m | bay/location_clean) ,
  family  = lognormal(),
  data= l_total_dat,
  # autocor = brms::cor_ar(~ julian_date | bay/location_clean:year)
)

save(lt_time, file = "~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/lt_time.Rdata")
load("~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/lt_time.Rdata")

pp_check(lt_time, nsamples = 50) + xlim(0,1000)
conditional_effects(lt_time)


# ============================================================
# 1) BAY-LEVEL slope draws for TEMP model
# ============================================================

temp_bay_slope_draws <- temp_time %>%
  spread_draws(
    b_n_year.m,
    r_bay[bay, n_year.m]
  ) %>%
  transmute(
    .draw,
    bay,
    slope_temp = b_n_year.m + r_bay
  )

# ============================================================
# 2) BAY-LEVEL slope draws for LARVAE model
# ============================================================

larv_bay_slope_draws <- lt_time %>%
  spread_draws(
    b_n_year.m,
    r_bay[bay, n_year.m]
  ) %>%
  transmute(
    .draw,
    bay,
    slope_larvae = b_n_year.m + r_bay
  )

# ============================================================
# 3) Summaries: median + 50% + 90% CrI for each bay
# ============================================================

temp_bay_slope_summ <- temp_bay_slope_draws %>%
  group_by(bay) %>%
  summarise(
    x_med = median(slope_temp),
    x_l50 = quantile(slope_temp, 0.25),
    x_u50 = quantile(slope_temp, 0.75),
    x_l90 = quantile(slope_temp, 0.05),
    x_u90 = quantile(slope_temp, 0.95),
    .groups = "drop"
  )

larv_bay_slope_summ <- larv_bay_slope_draws %>%
  group_by(bay) %>%
  summarise(
    y_med = median(slope_larvae),
    y_l50 = quantile(slope_larvae, 0.25),
    y_u50 = quantile(slope_larvae, 0.75),
    y_l90 = quantile(slope_larvae, 0.05),
    y_u90 = quantile(slope_larvae, 0.95),
    .groups = "drop"
  )

slope_xy_bay <- inner_join(temp_bay_slope_summ, larv_bay_slope_summ, by = "bay")

# ============================================================
# 4) Plot: x = temp bay slope, y = larvae bay slope
# ============================================================

ggplot(slope_xy_bay, aes(x = x_med, y = y_med, colour = bay)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  
  # 90% CrI (thin)
  geom_errorbarh(aes(xmin = x_l90, xmax = x_u90),
                 height = 0, linewidth = 0.7, alpha = 0.55) +
  geom_linerange(aes(ymin = y_l90, ymax = y_u90),
                 linewidth = 0.7, alpha = 0.55) +
  
  # 50% CrI (thick)
  geom_errorbarh(aes(xmin = x_l50, xmax = x_u50),
                 height = 0, linewidth = 2.1, alpha = 0.9) +
  geom_linerange(aes(ymin = y_l50, ymax = y_u50),
                 linewidth = 2.1, alpha = 0.9) +
  
  geom_point(size = 2.8) +
  scale_colour_viridis_d(option = "viridis") +
  labs(
    x = "Bay-specific temperature trend: d(water_temp)/d(n_year.m)",
    y = "Bay-specific larvae trend: d(log(larvae_total))/d(n_year.m)",
    title = "Bay-level coupling of warming and larvae trends",
    subtitle = "Points = posterior medians; thick = 50% CrI; thin = 90% CrI"
  ) +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")
