# how to clean your global environment
rm(list = ls())


library(tidyverse)
library(ggplot2)
library(viridis)
library(lubridate)
library(brms)
library(patchwork)
library(tidybayes)
library(scales)


# set your own personal working directory below
# Emma MBP
setwd("~/Dropbox/_Projects/PEI Oysters/Data/OMP")
# Maddy
setwd("~/Data/OMP/")


omp_dat <- read.csv("OMP_clean_2025.csv", header= TRUE)
head(omp_dat)


l_total_dat <- omp_dat %>%  # create numeric month day
  #group_by(bay, location_clean, f_year) %>% # group by location and year 
  #filter( larvae_above_250_microns == max(larvae_above_250_microns))  %>% # take max of larvae for each group specified above
  filter(larvae_total > 0) %>%      
  filter(larvae_above_250_microns > 0) %>% 
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

View(l_total_dat)
summary(l_total_dat)

# ==================================================
# MODEL
# ==================================================

#inspect high values
#l_total_dat %>% filter(water_temp == 722.00)
# l_total_dat %>% filter(water_temp == 99.90)
# l_total_dat %>% filter(water_temp == 34.00)

# lstotal_temp_time_sal <- brm(
#   larvae_above_250_microns ~ water_temp.m * n_year.m * salinity.m +
#     (1 + n_year.m * water_temp.m  || bay/location_clean) + (1 | julian_date.m),
#   data    = l_total_dat,
#   iter    = 5000,
#   warmup  = 1000,
#   family  = lognormal(),
#   #prior   = prior(normal(0, 0.5), class = "sigma", lb = 0),
#   control = list(adapt_delstA = 0.999, max_treedepth = 20)
# )


lstotal_temp_time_sal <- brm(
  larvae_above_250_microns ~ water_temp.m * n_year.m * salinity.m +
    s(julian_date,  k = 10) +
    (1 + n_year.m | bay/location_clean),
  data    = l_total_dat,
  family  = lognormal(),
  iter    = 5000, warmup = 1000,
  control = list(adapt_delta = 0.999, max_treedepth = 20)
)
# ================================================

# Emma's paths
save(lstotal_temp_time_sal, file = "~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/lstotal_temp_time_sal.Rdata")
load("~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/lstotal_temp_time_sal.Rdata")

# Maddy's path
#save(lstotal_temp_time_sal, file = "~/Data/Model_fits/OMP/lstotal_temp_time_sal.Rdata")
load("~/Data/Model_fits/OMP/lstotal_temp_time_sal.Rdata")


summary(lstotal_temp_time_sal)
pp_check(lstotal_temp_time_sal) + xlim(1,1500) 
conditional_effects(lstotal_temp_time_sal)

# ============================================================
# ALIGNMENT SETUP (shared constants)
# ============================================================

set.seed(123)

lstA_ref <- l_total_dat %>%
  ungroup() %>%
  summarise(
    wt_mean   = mean(water_temp, na.rm = TRUE),
    yr_mean   = mean(n_year,     na.rm = TRUE),
    sal_m_q25 = quantile(salinity.m, 0.25, na.rm = TRUE),
    sal_m_q50 = quantile(salinity.m, 0.50, na.rm = TRUE),
    sal_m_q75 = quantile(salinity.m, 0.75, na.rm = TRUE),
    wt_ref    = median(water_temp, na.rm = TRUE),
    jd_ref    = median(julian_date, na.rm = TRUE)   # <-- NEW (raw julian_date)
  )

# salinity slices on centered scale (matches model input)
lstA_sal_df <- tibble(
  salinity.m = c(lstA_ref$sal_m_q25, lstA_ref$sal_m_q50, lstA_ref$sal_m_q75),
  sal_label  = factor(
    c("Low salinity (25th pct)", "Median salinity", "High salinity (75th pct)"),
    levels = c("Low salinity (25th pct)", "Median salinity", "High salinity (75th pct)")
  )
)

# ============================================================
# lstA_W) TEMP ON X — prediction grid (year-specific temp support)
# ============================================================

lstA_W_grid <- l_total_dat %>%
  distinct(n_year, water_temp) %>%
  group_by(n_year) %>%
  summarise(
    water_temp = seq(
      quantile(water_temp, 0.05, na.rm = TRUE),
      quantile(water_temp, 0.95, na.rm = TRUE),
      length.out = 60
    ),
    .groups = "drop"
  ) %>%
  unnest(water_temp) %>%
  crossing(lstA_sal_df) %>%
  mutate(
    water_temp.m  = water_temp - lstA_ref$wt_mean,
    n_year.m      = n_year     - lstA_ref$yr_mean,
    julian_date   = lstA_ref$jd_ref,  
    year_group    = factor(n_year),
    row_id        = row_number()
  )

# posterior draws (population-level)
lstA_W_ep_mat <- posterior_epred(
  lstotal_temp_time_sal,
  newdata    = lstA_W_grid,
  re_formula = NA
)

# draw pool (memory-safe)
lstA_W_pool_size <- min(2000, nrow(lstA_W_ep_mat))
lstA_W_draw_ids  <- sample(seq_len(nrow(lstA_W_ep_mat)), size = lstA_W_pool_size)

# to long + join grid meta
lstA_W_ep_long <- lstA_W_ep_mat[lstA_W_draw_ids, ] %>%
  as_tibble() %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols = - .draw,
    names_to = "row_id",
    values_to = "epred"
  ) %>%
  mutate(row_id = readr::parse_number(row_id)) %>%
  left_join(
    lstA_W_grid %>% select(row_id, n_year, water_temp, year_group, sal_label),
    by = "row_id"
  )

# ribbons: 50% + 90% + median
lstA_W_summ <- lstA_W_ep_long %>%
  group_by(n_year, water_temp, year_group, sal_label) %>%
  summarise(
    med     = median(epred),
    low50   = quantile(epred, 0.25),
    up50    = quantile(epred, 0.75),
    low90   = quantile(epred, 0.05),
    up90    = quantile(epred, 0.95),
    .groups = "drop"
  )

# log-safe crop (consistent)
lstA_W_y_max_obs  <- max(l_total_dat$larvae_total, na.rm = TRUE)
lstA_W_y_min_plot <- max(1, min(lstA_W_summ$low90, na.rm = TRUE))
lstA_W_y_max_plot <- min(lstA_W_y_max_obs * 1.1, max(lstA_W_summ$up90, na.rm = TRUE))

# plot
lstA_W_fig_panel <- ggplot(
  lstA_W_summ,
  aes(x = water_temp, y = med, group = year_group, colour = year_group)
) +
  geom_ribbon(aes(ymin = low90, ymax = up90), alpha = 0.10, colour = NA) +
  geom_ribbon(aes(ymin = low50, ymax = up50), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ sal_label, nrow = 1) +
  scale_colour_viridis_d(option = "viridis" ) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30, 35)) +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(4, 8, 16, 32, 64, 128, 256, 512, 1024),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(lstA_W_y_min_plot, lstA_W_y_max_plot)) +
  labs(
    x = "Surface water temperature (°C)",
    y = "Total oyster larvae",
    #title = "Temperature–abundance relationships vary through time and with salinity",
    #subtitle = "Lines = posterior medians; ribbons = 50% (inner) and 90% (outer) CrI",
    colour= "Year"
    
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

lstA_W_fig_panel

# ============================================================
# lstA_WMED) Median salinity only — filtered draws pipeline
#   - spaghetti uses sampled filtered draws (for readability)
#   - slopes/intercepts use ALL filtered draws (stable)
# ============================================================

lstA_W_med_df <- tibble(
  salinity.m = lstA_ref$sal_m_q50,
  sal_label  = factor("Median salinity", levels = "Median salinity")
)

lstA_W_grid_med <- l_total_dat %>%
  distinct(n_year, water_temp) %>%
  group_by(n_year) %>%
  summarise(
    water_temp = seq(
      quantile(water_temp, 0.05, na.rm = TRUE),
      quantile(water_temp, 0.95, na.rm = TRUE),
      length.out = 60
    ),
    .groups = "drop"
  ) %>%
  unnest(water_temp) %>%
  crossing(lstA_W_med_df) %>%
  mutate(
    water_temp.m  = water_temp - lstA_ref$wt_mean,
    n_year.m      = n_year     - lstA_ref$yr_mean,
    julian_date   = lstA_ref$jd_ref,          # <-- REPLACE julian_date.m = 0
    year_group    = factor(n_year),
    row_id        = row_number()
  )

lstA_W_ep_mat_med <- posterior_epred(
  lstotal_temp_time_sal,
  newdata    = lstA_W_grid_med,
  re_formula = NA
)

lstA_W_pool_size_med <- min(2000, nrow(lstA_W_ep_mat_med))
lstA_W_draw_ids_med  <- sample(seq_len(nrow(lstA_W_ep_mat_med)), size = lstA_W_pool_size_med)

lstA_W_ep_long_med <- lstA_W_ep_mat_med[lstA_W_draw_ids_med, ] %>%
  as_tibble() %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols = - .draw,
    names_to = "row_id",
    values_to = "epred"
  ) %>%
  mutate(row_id = readr::parse_number(row_id)) %>%
  left_join(
    lstA_W_grid_med %>% select(row_id, n_year, water_temp, year_group, sal_label),
    by = "row_id"
  )

# envelope 10–90 at each x-point (row_id), then keep draws with >=95% inside
lstA_W_bounds80_med <- lstA_W_ep_long_med %>%
  group_by(row_id) %>%
  summarise(
    q10 = quantile(epred, 0.10, na.rm = TRUE),
    q90 = quantile(epred, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

lstA_W_draw_keep_med <- lstA_W_ep_long_med %>%
  left_join(lstA_W_bounds80_med, by = "row_id") %>%
  mutate(in80 = epred >= q10 & epred <= q90) %>%
  group_by(.draw) %>%
  summarise(prop_in80 = mean(in80, na.rm = TRUE), .groups = "drop") %>%
  filter(prop_in80 >= 0.95)

# ALL kept draws (for slopes/intercepts + stable ribbons/medians)
lstA_W_ep_keep_med <- lstA_W_ep_long_med %>%
  semi_join(lstA_W_draw_keep_med, by = ".draw")

# Sample a subset ONLY for spaghetti plotting (not for slopes/intercepts)
lstA_W_n_keep_med <- lstA_W_draw_keep_med %>% summarise(n = n()) %>% pull(n)
lstA_W_n_spag_med <- min(150, lstA_W_n_keep_med)

lstA_W_spag_ids_med <- lstA_W_draw_keep_med %>%
  slice_sample(n = lstA_W_n_spag_med) %>%
  pull(.draw)

lstA_W_ep_spag_med <- lstA_W_ep_keep_med %>%
  filter(.draw %in% lstA_W_spag_ids_med)

# thick lines (medians) computed from ALL kept draws (aligned)
lstA_W_summ_med <- lstA_W_ep_keep_med %>%
  group_by(n_year, water_temp, year_group, sal_label) %>%
  summarise(
    med   = median(epred),
    low90 = quantile(epred, 0.05),
    up90  = quantile(epred, 0.95),
    .groups = "drop"
  )

# log-safe crop (aligned)
lstA_W_y_max_obs_med  <- max(l_total_dat$larvae_total, na.rm = TRUE)
lstA_W_y_min_plot_med <- max(1, min(lstA_W_summ_med$low90, na.rm = TRUE))
lstA_W_y_max_plot_med <- min(lstA_W_y_max_obs_med * 1.1, max(lstA_W_summ_med$up90, na.rm = TRUE))

# spaghetti figure (thin filtered draws + thick medians FROM SAME kept draw set)
lstA_W_fig_spag_med <- ggplot() +
  geom_line(
    data = lstA_W_ep_spag_med,
    aes(x = water_temp, y = epred, group = interaction(.draw, year_group), colour = year_group),
    linewidth = 0.4,
    alpha = 0.15
  ) +
  geom_line(
    data = lstA_W_summ_med,
    aes(x = water_temp, y = med, group = year_group, colour = year_group),
    linewidth = 0.8
  ) +
  facet_wrap(~ sal_label, nrow = 1) +
  scale_colour_viridis_d(option = "viridis") +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30, 35)) +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(4, 8, 16, 32, 64, 128, 256, 512, 1024),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(lstA_W_y_min_plot_med, lstA_W_y_max_plot_med)) +
  labs(
    x = "Surface water temperature (°C)",
    y = "Total oyster larvae",
    title = "Temperature–abundance relationships at median salinity",
    subtitle = "Thin lines = filtered posterior draws; thick lines = posterior medians (all from the same kept-draw set)"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

lstA_W_fig_spag_med


# ============================================================
# INTERCEPTS + SLOPES (ALIGNED: computed from ALL kept draws)
# ============================================================

# intercept-like: predicted larvae at reference temperature (median observed water temp)
lstA_temp_ref <- lstA_ref$wt_ref

lstA_intercept_draws <- lstA_W_ep_keep_med %>%
  group_by(.draw, n_year, year_group) %>%
  mutate(dist_to_ref = abs(water_temp - lstA_temp_ref)) %>%
  slice_min(order_by = dist_to_ref, n = 1, with_ties = FALSE) %>%
  ungroup()

lstA_intercept_summ <- lstA_intercept_draws %>%
  group_by(n_year, year_group) %>%
  summarise(
    estimate = median(epred),
    lower50  = quantile(epred, 0.25),
    upper50  = quantile(epred, 0.75),
    lower90  = quantile(epred, 0.05),
    upper90  = quantile(epred, 0.95),
    .groups = "drop"
  )

lstA_fig_intercepts <- ggplot(lstA_intercept_summ, aes(x = n_year, y = estimate)) +
  geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.6) +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(4, 8, 16, 32, 64, 128, 256, 512, 1024),
    labels = scales::label_number()
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  labs(
    x = "Monitoring year",
    y = "Predicted total larvae at reference temperature",
    title = "Intercept-like differences through time (median salinity)",
    subtitle = paste0(
      "Reference temperature = median observed water temp (",
      round(lstA_temp_ref, 1), " °C). Thick = 50% CrI; thin = 90% CrI."
    )
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank())

lstA_fig_intercepts


# slope: within-year temperature sensitivity on log10 scale (stable for lognormal)
lstA_slope_draws <- lstA_W_ep_keep_med %>%
  mutate(log10_epred = log10(epred)) %>%
  group_by(.draw, n_year, year_group) %>%
  summarise(
    slope_log10 = cov(water_temp, log10_epred) / var(water_temp),
    .groups = "drop"
  )

lstA_slope_summ <- lstA_slope_draws %>%
  group_by(n_year, year_group) %>%
  summarise(
    estimate = median(slope_log10),
    lower50  = quantile(slope_log10, 0.25),
    upper50  = quantile(slope_log10, 0.75),
    lower90  = quantile(slope_log10, 0.05),
    upper90  = quantile(slope_log10, 0.95),
    .groups = "drop"
  )

lstA_fig_slopes <- ggplot(lstA_slope_summ, aes(x = n_year, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.6) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  labs(
    x = "Monitoring year",
    y = "Slope: d(log10(larvae))/d(°C)",
    title = "Temperature sensitivity through time (median salinity)",
    subtitle = "Computed from the same kept-draw set used for the spaghetti plot. Thick = 50% CrI; thin = 90% CrI."
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank())

lstA_fig_slopes

# Optional combined view (patchwork)
lstA_W_fig_spag_med + lstA_fig_intercepts + lstA_fig_slopes

# ============================================================
# lstA_T) YEAR ON X × TEMP LEVELS @ MEDIAN SALINITY (aligned)
# ============================================================

lstA_T_temp_df <- l_total_dat %>%
  ungroup() %>%
  summarise(
    q25 = quantile(water_temp, 0.25, na.rm = TRUE),
    q50 = quantile(water_temp, 0.50, na.rm = TRUE),
    q75 = quantile(water_temp, 0.75, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = everything(), names_to = "q", values_to = "water_temp") %>%
  mutate(
    temp_label = factor(
      c("Cool temp (25th pct)", "Median temp", "Warm temp (75th pct)"),
      levels = c("Cool temp (25th pct)", "Median temp", "Warm temp (75th pct)")
    )
  ) %>%
  select(water_temp, temp_label)


lstA_T_year_seq <- seq(
  floor(min(l_total_dat$n_year, na.rm = TRUE)),
  ceiling(max(l_total_dat$n_year, na.rm = TRUE)),
  by = 1
)

lstA_T_grid <- crossing(
  n_year = lstA_T_year_seq,
  lstA_T_temp_df,
  tibble(salinity.m = lstA_ref$sal_m_q50, sal_label = factor("Median salinity"))
) %>%
  mutate(
    n_year.m      = n_year     - lstA_ref$yr_mean,
    water_temp.m  = water_temp - lstA_ref$wt_mean,
    julian_date   = lstA_ref$jd_ref,          # <-- REPLACE julian_date.m = 0
    row_id        = row_number()
  )

lstA_T_ep_mat <- posterior_epred(
  lstotal_temp_time_sal,
  newdata    = lstA_T_grid,
  re_formula = NA
)

lstA_T_pool_size <- min(2000, nrow(lstA_T_ep_mat))
lstA_T_draw_ids  <- sample(seq_len(nrow(lstA_T_ep_mat)), size = lstA_T_pool_size)

lstA_T_ep_long <- lstA_T_ep_mat[lstA_T_draw_ids, ] %>%
  as_tibble() %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(cols = - .draw, names_to = "row_id", values_to = "epred") %>%
  mutate(row_id = readr::parse_number(row_id)) %>%
  left_join(
    lstA_T_grid %>% select(row_id, n_year, temp_label, sal_label),
    by = "row_id"
  )

# unfiltered summary ribbons (50% + 90%)
lstA_T_summ <- lstA_T_ep_long %>%
  group_by(n_year, temp_label, sal_label) %>%
  summarise(
    med    = median(epred),
    low50  = quantile(epred, 0.25),
    up50   = quantile(epred, 0.75),
    low90  = quantile(epred, 0.05),
    up90   = quantile(epred, 0.95),
    .groups = "drop"
  )

# crop on log scale
lstA_T_y_max_obs  <- max(l_total_dat$larvae_total, na.rm = TRUE)
lstA_T_y_min_plot <- max(1, min(lstA_T_summ$low90, na.rm = TRUE))
lstA_T_y_max_plot <- min(lstA_T_y_max_obs * 1.1, max(lstA_T_summ$up90, na.rm = TRUE))

lstA_T_fig_time_temp <- ggplot(
  lstA_T_summ,
  aes(x = n_year, y = med, colour = temp_label, group = temp_label)
) +
  geom_ribbon(aes(ymin = low90, ymax = up90, fill = temp_label), alpha = 0.10, colour = NA) +
  geom_ribbon(aes(ymin = low50, ymax = up50, fill = temp_label), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.1) +
  # facet_wrap(~ sal_label, nrow = 1) +
  # scale_colour_viridis_d(option = "viridis") +
  #scale_fill_viridis_d(option = "viridis") +
  # --- Dark2 palette with explicit mapping ---
  scale_colour_manual(values = c( "Cool temp (25th pct)" = "#1B9E77",  # green (Dark2)
                                  "Median temp"          = "#7570B3",  # purple (Dark2)
                                  "Warm temp (75th pct)" = "#D95F02"   # orange (Dark2)
  ) ) +
  scale_fill_manual(values = c( "Cool temp (25th pct)" = "#1B9E77",  # green (Dark2)
                                "Median temp"          = "#7570B3",  # purple (Dark2)
                                "Warm temp (75th pct)" = "#D95F02"   # orange (Dark2)
  ) ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(4, 8, 16, 32, 64, 128, 256, 512, 1024),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(lstA_T_y_min_plot, lstA_T_y_max_plot)) +
  labs(
    x = "Monitoring year",
    y = "Total oyster larvae",
    subtitle= "a)", colour= "Temp", fill= "Temp"
    # title = "Temporal trends vary with temperature (median salinity)",
    # subtitle = "Lines = posterior medians; ribbons = 50% (inner) and 90% (outer) CrI"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

lstA_T_fig_time_temp


# ============================================================
# lstA_T) INTERCEPTS + SLOPES for COOL / MEDIAN / WARM temps
# - uses the SAME lstA_T_ep_long object from the year-on-x workflow
# - tidyverse only + no custom functions + no if/else
# - uncertainty: 50% + 90% CrI everywhere
# ============================================================

# ------------------------------------------------------------
# 0) Use your existing temp levels table (must have: water_temp, temp_label)
#     and your year sequence grid already used for lstA_T_ep_long
# ------------------------------------------------------------
# Assumes these already exist from the aligned lstA_T workflow:
#   lstA_T_temp_df, lstA_T_summ, lstA_T_ep_long, lstA_ref, lstA_T_year_seq

# ------------------------------------------------------------
# 1) INTERCEPTS: predicted larvae at a reference year (mean year)
#    - within each draw × temp_label, take epred at year closest to ref
# ------------------------------------------------------------
lstA_T_ref_year <- l_total_dat %>%
  ungroup() %>%
  summarise(ref_year = mean(n_year, na.rm = TRUE)) %>%
  pull(ref_year)

lstA_T_intercept_draws <- lstA_T_ep_long %>%
  group_by(.draw, temp_label) %>%
  mutate(dist_to_ref_year = abs(n_year - lstA_T_ref_year)) %>%
  slice_min(order_by = dist_to_ref_year, n = 1, with_ties = FALSE) %>%
  ungroup()

lstA_T_intercept_summ <- lstA_T_intercept_draws %>%
  group_by(temp_label) %>%
  summarise(
    estimate = median(epred),
    lower50  = quantile(epred, 0.25),
    upper50  = quantile(epred, 0.75),
    lower90  = quantile(epred, 0.05),
    upper90  = quantile(epred, 0.95),
    .groups  = "drop"
  )

lstA_T_fig_intercepts <- ggplot(
  lstA_T_intercept_summ,
  aes(x = temp_label, y = estimate, colour = temp_label)
) +
  geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.9, alpha = 0.6) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.2) +
  geom_point(size = 3) +
  #scale_colour_viridis_d(option = "viridis") +
  # --- Dark2 palette with explicit mapping ---
  scale_colour_manual(values = c( "Cool temp (25th pct)" = "#1B9E77",  # green (Dark2)
                                  "Median temp"          = "#7570B3",  # purple (Dark2)
                                  "Warm temp (75th pct)" = "#D95F02"   # orange (Dark2)
  ) ) +
  #coord_flip() +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(4, 8, 16, 32, 64, 128, 256, 512, 1024),
    labels = scales::label_number()
  ) +
  labs(
    x = NULL,
    y = "Predicted larvae at reference year",
    #title = "Baseline larvae abundance by temperature (median salinity)",
    subtitle = "b)"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x  = element_blank(),
    legend.position = "none"
  )

lstA_T_fig_intercepts


# ------------------------------------------------------------
# 2) SLOPES: change over time within each temp level
#    - compute slope within each draw × temp_label using cov/var
#    - (optional) do on log10 scale for stability/interpretability (recommended)
# ------------------------------------------------------------
lstA_T_slope_draws <- lstA_T_ep_long %>%
  mutate(log10_epred = log10(epred)) %>%
  group_by(.draw, temp_label) %>%
  summarise(
    slope_log10 = cov(n_year, log10_epred) / var(n_year),
    .groups = "drop"
  )

lstA_T_slope_summ <- lstA_T_slope_draws %>%
  group_by(temp_label) %>%
  summarise(
    estimate = median(slope_log10),
    lower50  = quantile(slope_log10, 0.25),
    upper50  = quantile(slope_log10, 0.75),
    lower90  = quantile(slope_log10, 0.05),
    upper90  = quantile(slope_log10, 0.95),
    .groups  = "drop"
  )

view(lstA_T_slope_summ)

lstA_T_fig_slopes <- ggplot(
  lstA_T_slope_summ,
  aes(x = temp_label, y = estimate, colour = temp_label)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.9, alpha = 0.6) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.2) +
  geom_point(size = 3) +
  #scale_colour_viridis_d(option = "viridis") +
  # --- Dark2 palette with explicit mapping ---
  scale_colour_manual(values = c( "Cool temp (25th pct)" = "#1B9E77",  # green (Dark2)
                                  "Median temp"          = "#7570B3",  # purple (Dark2)
                                  "Warm temp (75th pct)" = "#D95F02"   # orange (Dark2)
  ) ) +
  #coord_flip() +
  labs(
    x = NULL,
    y = "Slope:(log(larvae))/(year)",
    #title = "Temporal trends by temperature (median salinity)",
    subtitle = "c)"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x  = element_blank(),
    legend.position = "none"
  )

lstA_T_fig_slopes


# Optional: combine with your year-on-x ribbon plot (patchwork)
lstA_T_fig_time_temp + lstA_T_fig_intercepts + lstA_T_fig_slopes



g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

abund_temp_legend <- g_legend(lstA_T_fig_time_temp)

bottom_row_abund <- ( (lstA_T_fig_time_temp + theme(legend.position="none") ) + lstA_T_fig_intercepts + lstA_T_fig_slopes)/abund_temp_legend + plot_layout(heights = c(10, 1))

top_row_abund <- ( lstA_W_fig_spag_med + lstA_fig_intercepts + lstA_fig_slopes)


top_row_abund / bottom_row_abund
