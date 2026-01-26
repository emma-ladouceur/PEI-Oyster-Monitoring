library(tidyverse)
library(brms)
library(purrr)

library(tidybayes)
library(dplyr)



# set your own personal working directory below
# Emma MBP
setwd("~/Dropbox/_Projects/PEI Oysters/Data/OMP")
# Maddy
setwd("~/Data/OMP/")


omp_dat <- read.csv("OMP_clean_2025.csv", header= TRUE)
head(omp_dat)

omp_dat %>% select(bay, location_clean) %>% distinct() %>% arrange(bay, location_clean) 

l_total_dat <- omp_dat %>%  # create numeric month day
  group_by(bay, location_clean, f_year) %>% # group by location and year 
  # filter( larvae_total == max(larvae_total))  %>% # take max of larvae for each group specified above
  filter(larvae_total > 0) %>%      
  mutate(n_month_day = as.numeric(month_day)) %>%
  mutate(bay = as.factor(bay)) %>% 
  mutate(n_month_day = as.numeric(n_month_day)) %>% #at least 3 julian dates-
  group_by(location_clean, n_year) %>%
  filter(n_distinct(julian_date) > 3) %>%
  ungroup() %>%
  # 2) now require each location to have at least 3 years remaining
  group_by(bay, location_clean) %>%
  filter(n_distinct(n_year) >= 3) %>%
  ungroup() %>%
  group_by(bay, location_clean) %>%
  #center continuous variables for modelling
  mutate(water_temp.m = water_temp - mean(water_temp, na.rm = TRUE),
         salinity.m = salinity - mean(salinity, na.rm = TRUE),
         n_year.m = n_year - mean(n_year, na.rm = TRUE),
         julian_date.m = julian_date - mean(julian_date, na.rm = TRUE)
  ) %>% ungroup() %>%
  mutate(log_larvae_total = log(larvae_total))  
# (optional) if you want to run faster for many locations
# options(mc.cores = parallel::detectCores())

l_total_dat %>% select(bay, location_clean, n_year, julian_date) %>% distinct() %>% arrange(bay, location_clean, n_year, julian_date) %>% filter(location_clean == "Bentick Cove")

print(l_total_dat %>% select(bay, location_clean) %>% distinct() %>% arrange(bay, location_clean), n=Inf)  # 26 locations
l_total_dat %>% select(bay) %>% distinct() %>% arrange(bay) # 9 bays



# 1) model formula (keep it readable)
temp_time_form <- bf(
  water_temp ~ n_year.m #+ s(julian_date, k = 10)
)

# 2) fit separate model per location_clean
temp_time_by_loc <- l_total_dat %>%
  filter(!is.na(location_clean)) %>%
  group_by(location_clean) %>%
  nest() %>%
  mutate(
    fit = map(
      data,
      ~ brm(
        formula = temp_time_form,
        data    = .x,
        family  = student(),
        chains  = 4, iter = 4000, warmup = 1000,
        control = list(adapt_delta = 0.95)
        # autocor = brms::cor_ar(~ julian_date | bay/location_clean:year)
      )
    )
  ) %>%
  ungroup()

temp_time_by_loc
saveRDS(temp_time_by_loc, file = "~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/temp_time_by_loc.rds")
temp_time_by_loc <- readRDS("~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/temp_time_by_loc.rds")



# 1) model formula (keep it readable)
abund_time_form <- bf(
  larvae_total ~ n_year.m #+ s(julian_date, k = 10)
)

# 2) fit separate model per location_clean
abund_time_by_loc <- l_total_dat %>%
  filter(!is.na(location_clean)) %>%
  group_by(location_clean) %>%
  nest() %>%
  mutate(
    fit = map(
      data,
      ~ brm(
        formula = abund_time_form,
        data    = .x,
        family  = lognormal(),
        chains  = 4, iter = 4000, warmup = 1000,
        control = list(adapt_delta = 0.95)
        # autocor = brms::cor_ar(~ julian_date | bay/location_clean:year)
      )
    )
  ) %>%
  ungroup()

abund_time_by_loc
saveRDS(abund_time_by_loc, file = "~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/tabund_time_by_loc.rds")
abund_time_by_loc <- readRDS("~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/abund_time_by_loc.rds")




temp_slopes <- temp_time_by_loc %>%
  mutate(
    draws = map(
      fit,
      ~ spread_draws(.x, b_n_year.m) %>%
        transmute(
          .draw,
          temp_slope = b_n_year.m
        )
    )
  ) %>%
  select(location_clean, draws) %>%
  unnest(draws)

abund_slopes <- abund_time_by_loc %>%
  mutate(
    draws = map(
      fit,
      ~ spread_draws(.x, b_n_year.m) %>%
        transmute(
          .draw,
          abund_slope = b_n_year.m
        )
    )
  ) %>%
  select(location_clean, draws) %>%
  unnest(draws)

slope_joint <- temp_slopes %>%
  inner_join(
    abund_slopes,
    by = c("location_clean", ".draw")
  )

slope_joint %>%
  summarise(
    cor_median = median(cor(temp_slope, abund_slope)),
    cor_l90    = quantile(cor(temp_slope, abund_slope), 0.05),
    cor_u90    = quantile(cor(temp_slope, abund_slope), 0.95)
  )


slope_meta_mod <- brm(
  abund_slope ~ temp_slope,
  data   = slope_joint,
  family = gaussian(),
  chains = 4, iter = 4000, warmup = 1000,
  control = list(adapt_delta = 0.95)
)

summary(slope_meta_mod)


slope_meta_draws <- slope_meta_mod %>%
  spread_draws(b_Intercept, b_temp_slope) %>%
  summarise(
    intercept_med = median(b_Intercept),
    intercept_l90 = quantile(b_Intercept, 0.05),
    intercept_u90 = quantile(b_Intercept, 0.95),
    
    slope_med = median(b_temp_slope),
    slope_l90 = quantile(b_temp_slope, 0.05),
    slope_u90 = quantile(b_temp_slope, 0.95)
  )

slope_meta_draws


library(dplyr)
library(tidyr)
library(ggplot2)
library(tidybayes)

# summarize per-location posterior for each slope
temp_summ <- temp_slopes %>%
  group_by(location_clean) %>%
  summarise(
    temp_med = median(temp_slope),
    temp_l50 = quantile(temp_slope, 0.25),
    temp_u50 = quantile(temp_slope, 0.75),
    temp_l90 = quantile(temp_slope, 0.05),
    temp_u90 = quantile(temp_slope, 0.95),
    .groups = "drop"
  )

abund_summ <- abund_slopes %>%
  group_by(location_clean) %>%
  summarise(
    abund_med = median(abund_slope),
    abund_l50 = quantile(abund_slope, 0.25),
    abund_u50 = quantile(abund_slope, 0.75),
    abund_l90 = quantile(abund_slope, 0.05),
    abund_u90 = quantile(abund_slope, 0.95),
    .groups = "drop"
  )

loc_summ <- temp_summ %>%
  inner_join(abund_summ, by = "location_clean")


meta_line <- slope_joint %>%
  summarise(
    x_min = quantile(temp_slope, 0.02),
    x_max = quantile(temp_slope, 0.98)
  ) %>%
  tidyr::expand_grid(
    temp_slope = seq(x_min, x_max, length.out = 120)
  ) %>%
  select(temp_slope) %>%
  add_epred_draws(
    slope_meta_mod,
    newdata = .,
    re_formula = NA
  ) %>%
  group_by(temp_slope) %>%
  summarise(
    y_med = median(.epred),
    y_l90 = quantile(.epred, 0.05),
    y_u90 = quantile(.epred, 0.95),
    .groups = "drop"
  )

slope_slope_fig <- ggplot() +
  # meta-regression uncertainty band (90%)
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") + 
  geom_ribbon(
    data = meta_line,
    aes(x = temp_slope, ymin = y_l90, ymax = y_u90),
    alpha = 0.18
  ) +
  # meta-regression median line
  geom_line(
    data = meta_line,
    aes(x = temp_slope, y = y_med),
    linewidth = 1
  ) +
  
  # location uncertainty (90% then 50% so it “reads” nicely)
  geom_linerange(
    data = loc_summ,
    aes(x = temp_med, ymin = abund_l90, ymax = abund_u90),
    alpha = 0.35
  ) +
  geom_linerange(
    data = loc_summ,
    aes(x = temp_med, ymin = abund_l50, ymax = abund_u50),
    linewidth = 2
  ) +
  geom_errorbarh(
    data = loc_summ,
    aes(y = abund_med, xmin = temp_l90, xmax = temp_u90),
    height = 0,
    alpha = 0.35
  ) +
  geom_errorbarh(
    data = loc_summ,
    aes(y = abund_med, xmin = temp_l50, xmax = temp_u50),
    height = 0,
    linewidth = 2
  ) +
  geom_point(
    data = loc_summ,
    aes(x = temp_med, y = abund_med),
    size = 2.7
  ) +
  
  labs(
    x = "Temperature trend (slope vs centered year)",
    y = "Larval abundance trend (slope vs centered year)",
    title = "Do warming trends explain temporal changes in larval abundance?",
    subtitle = "Points = location posterior medians; thick bars = 50% CI; thin bars = 90% CI. Line/band = meta-regression median and 90% CI."
  )

slope_slope_fig

