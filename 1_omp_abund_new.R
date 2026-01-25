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

View(l_total_dat)
summary(l_total_dat)

# ==================================================
# MODEL
# ==================================================

#inspect high values
#l_total_dat %>% filter(water_temp == 722.00)
# l_total_dat %>% filter(water_temp == 99.90)
# l_total_dat %>% filter(water_temp == 34.00)

# ltotal_temp_time_sal <- brm(
#   larvae_total ~ water_temp.m * n_year.m * salinity.m +
#     (1 + n_year.m * water_temp.m  || bay/location_clean) + (1 | julian_date.m),
#   data    = l_total_dat,                
#   iter    = 5000,
#   warmup  = 1000,
#   family  = lognormal(),
#   #prior   = prior(normal(0, 0.5), class = "sigma", lb = 0),
#   control = list(adapt_delta = 0.999, max_treedepth = 20)
# )

# ================================================

# Emma's paths
save(ltotal_temp_time_sal, file = "~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/ltotal_temp_time_sal.Rdata")
load("~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/ltotal_temp_time_sal.Rdata")

# Maddy's path
#save(ltotal_temp_time_sal, file = "~/Data/Model_fits/OMP/ltotal_temp_time_sal.Rdata")
load("~/Data/Model_fits/OMP/ltotal_temp_time_sal.Rdata")


summary(ltotal_temp_time_sal)
pp_check(ltotal_temp_time_sal) + xlim(1,1500) 
conditional_effects(ltotal_temp_time_sal)




# ============================================================
# ALIGNMENT SETUP (shared constants)
# ============================================================

set.seed(123)

ltA_ref <- l_total_dat %>%
  ungroup() %>%
  summarise(
    wt_mean   = mean(water_temp, na.rm = TRUE),
    yr_mean   = mean(n_year,     na.rm = TRUE),
    sal_m_q25 = quantile(salinity.m, 0.25, na.rm = TRUE),
    sal_m_q50 = quantile(salinity.m, 0.50, na.rm = TRUE),
    sal_m_q75 = quantile(salinity.m, 0.75, na.rm = TRUE),
    wt_ref    = median(water_temp, na.rm = TRUE)
  )

# salinity slices on centered scale (matches model input)
ltA_sal_df <- tibble(
  salinity.m = c(ltA_ref$sal_m_q25, ltA_ref$sal_m_q50, ltA_ref$sal_m_q75),
  sal_label  = factor(
    c("Low salinity (25th pct)", "Median salinity", "High salinity (75th pct)"),
    levels = c("Low salinity (25th pct)", "Median salinity", "High salinity (75th pct)")
  )
)

# ============================================================
# ltA_W) TEMP ON X — prediction grid (year-specific temp support)
# ============================================================

ltA_W_grid <- l_total_dat %>%
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
  crossing(ltA_sal_df) %>%
  mutate(
    water_temp.m  = water_temp - ltA_ref$wt_mean,
    n_year.m      = n_year     - ltA_ref$yr_mean,
    julian_date.m = 0,
    year_group    = factor(n_year),
    row_id        = row_number()
  )

# posterior draws (population-level)
ltA_W_ep_mat <- posterior_epred(
  ltotal_temp_time_sal,
  newdata    = ltA_W_grid,
  re_formula = NA
)

# draw pool (memory-safe)
ltA_W_pool_size <- min(2000, nrow(ltA_W_ep_mat))
ltA_W_draw_ids  <- sample(seq_len(nrow(ltA_W_ep_mat)), size = ltA_W_pool_size)

# to long + join grid meta
ltA_W_ep_long <- ltA_W_ep_mat[ltA_W_draw_ids, ] %>%
  as_tibble() %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols = - .draw,
    names_to = "row_id",
    values_to = "epred"
  ) %>%
  mutate(row_id = readr::parse_number(row_id)) %>%
  left_join(
    ltA_W_grid %>% select(row_id, n_year, water_temp, year_group, sal_label),
    by = "row_id"
  )

# ribbons: 50% + 90% + median
ltA_W_summ <- ltA_W_ep_long %>%
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
ltA_W_y_max_obs  <- max(l_total_dat$larvae_total, na.rm = TRUE)
ltA_W_y_min_plot <- max(1, min(ltA_W_summ$low90, na.rm = TRUE))
ltA_W_y_max_plot <- min(ltA_W_y_max_obs * 1.1, max(ltA_W_summ$up90, na.rm = TRUE))

# plot
ltA_W_fig_panel <- ggplot(
  ltA_W_summ,
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
  coord_cartesian(ylim = c(ltA_W_y_min_plot, ltA_W_y_max_plot)) +
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

ltA_W_fig_panel

# ============================================================
# ltA_WMED) Median salinity only — filtered draws pipeline
#   - spaghetti uses sampled filtered draws (for readability)
#   - slopes/intercepts use ALL filtered draws (stable)
# ============================================================

ltA_W_med_df <- tibble(
  salinity.m = ltA_ref$sal_m_q50,
  sal_label  = factor("Median salinity", levels = "Median salinity")
)

ltA_W_grid_med <- l_total_dat %>%
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
  crossing(ltA_W_med_df) %>%
  mutate(
    water_temp.m  = water_temp - ltA_ref$wt_mean,
    n_year.m      = n_year     - ltA_ref$yr_mean,
    julian_date.m = 0,
    year_group    = factor(n_year),
    row_id        = row_number()
  )

ltA_W_ep_mat_med <- posterior_epred(
  ltotal_temp_time_sal,
  newdata    = ltA_W_grid_med,
  re_formula = NA
)

ltA_W_pool_size_med <- min(2000, nrow(ltA_W_ep_mat_med))
ltA_W_draw_ids_med  <- sample(seq_len(nrow(ltA_W_ep_mat_med)), size = ltA_W_pool_size_med)

ltA_W_ep_long_med <- ltA_W_ep_mat_med[ltA_W_draw_ids_med, ] %>%
  as_tibble() %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols = - .draw,
    names_to = "row_id",
    values_to = "epred"
  ) %>%
  mutate(row_id = readr::parse_number(row_id)) %>%
  left_join(
    ltA_W_grid_med %>% select(row_id, n_year, water_temp, year_group, sal_label),
    by = "row_id"
  )

# envelope 10–90 at each x-point (row_id), then keep draws with >=95% inside
ltA_W_bounds80_med <- ltA_W_ep_long_med %>%
  group_by(row_id) %>%
  summarise(
    q10 = quantile(epred, 0.10, na.rm = TRUE),
    q90 = quantile(epred, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

ltA_W_draw_keep_med <- ltA_W_ep_long_med %>%
  left_join(ltA_W_bounds80_med, by = "row_id") %>%
  mutate(in80 = epred >= q10 & epred <= q90) %>%
  group_by(.draw) %>%
  summarise(prop_in80 = mean(in80, na.rm = TRUE), .groups = "drop") %>%
  filter(prop_in80 >= 0.95)

# ALL kept draws (for slopes/intercepts + stable ribbons/medians)
ltA_W_ep_keep_med <- ltA_W_ep_long_med %>%
  semi_join(ltA_W_draw_keep_med, by = ".draw")

# Sample a subset ONLY for spaghetti plotting (not for slopes/intercepts)
ltA_W_n_keep_med <- ltA_W_draw_keep_med %>% summarise(n = n()) %>% pull(n)
ltA_W_n_spag_med <- min(150, ltA_W_n_keep_med)

ltA_W_spag_ids_med <- ltA_W_draw_keep_med %>%
  slice_sample(n = ltA_W_n_spag_med) %>%
  pull(.draw)

ltA_W_ep_spag_med <- ltA_W_ep_keep_med %>%
  filter(.draw %in% ltA_W_spag_ids_med)

# thick lines (medians) computed from ALL kept draws (aligned)
ltA_W_summ_med <- ltA_W_ep_keep_med %>%
  group_by(n_year, water_temp, year_group, sal_label) %>%
  summarise(
    med   = median(epred),
    low90 = quantile(epred, 0.05),
    up90  = quantile(epred, 0.95),
    .groups = "drop"
  )

# log-safe crop (aligned)
ltA_W_y_max_obs_med  <- max(l_total_dat$larvae_total, na.rm = TRUE)
ltA_W_y_min_plot_med <- max(1, min(ltA_W_summ_med$low90, na.rm = TRUE))
ltA_W_y_max_plot_med <- min(ltA_W_y_max_obs_med * 1.1, max(ltA_W_summ_med$up90, na.rm = TRUE))

# spaghetti figure (thin filtered draws + thick medians FROM SAME kept draw set)
ltA_W_fig_spag_med <- ggplot() +
  geom_line(
    data = ltA_W_ep_spag_med,
    aes(x = water_temp, y = epred, group = interaction(.draw, year_group), colour = year_group),
    linewidth = 0.4,
    alpha = 0.15
  ) +
  geom_line(
    data = ltA_W_summ_med,
    aes(x = water_temp, y = med, group = year_group, colour = year_group),
    linewidth = 0.8
  ) +
  facet_wrap(~ sal_label, nrow = 1) +
  scale_colour_viridis_d(option = "viridis", name = "Monitoring year") +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30, 35)) +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(4, 8, 16, 32, 64, 128, 256, 512, 1024),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(ltA_W_y_min_plot_med, ltA_W_y_max_plot_med)) +
  labs(
    x = "Surface water temperature (°C)",
    y = "Total oyster larvae",
    # title = "Temperature–abundance relationships at median salinity",
    subtitle = "a)"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ltA_W_fig_spag_med


# ============================================================
# INTERCEPTS + SLOPES (ALIGNED: computed from ALL kept draws)
# ============================================================

# intercept-like: predicted larvae at reference temperature (median observed water temp)
ltA_temp_ref <- ltA_ref$wt_ref

ltA_intercept_draws <- ltA_W_ep_keep_med %>%
  group_by(.draw, n_year, year_group) %>%
  mutate(dist_to_ref = abs(water_temp - ltA_temp_ref)) %>%
  slice_min(order_by = dist_to_ref, n = 1, with_ties = FALSE) %>%
  ungroup()

ltA_intercept_summ <- ltA_intercept_draws %>%
  group_by(n_year, year_group) %>%
  summarise(
    estimate = median(epred),
    lower50  = quantile(epred, 0.25),
    upper50  = quantile(epred, 0.75),
    lower90  = quantile(epred, 0.05),
    upper90  = quantile(epred, 0.95),
    .groups = "drop"
  )

ltA_fig_intercepts <- ggplot(ltA_intercept_summ, aes(x = n_year, y = estimate, colour = year_group)) +
  geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.6) +
  scale_colour_viridis_d(option = "viridis", guide = "none") +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(4, 8, 16, 32, 64, 128, 256, 512, 1024),
    labels = scales::label_number()
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  labs(
    x = "Monitoring year",
    y = "Predicted total larvae at reference temperature",
    # title = "Intercept-like differences through time (median salinity)",
    subtitle = "b)"
    ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank())

ltA_fig_intercepts


# slope: within-year temperature sensitivity on log10 scale (stable for lognormal)
ltA_slope_draws <- ltA_W_ep_keep_med %>%
  mutate(log10_epred = log10(epred)) %>%
  group_by(.draw, n_year, year_group) %>%
  summarise(
    slope_log10 = cov(water_temp, log10_epred) / var(water_temp),
    .groups = "drop"
  )

ltA_slope_summ <- ltA_slope_draws %>%
  group_by(n_year, year_group) %>%
  summarise(
    estimate = median(slope_log10),
    lower50  = quantile(slope_log10, 0.25),
    upper50  = quantile(slope_log10, 0.75),
    lower90  = quantile(slope_log10, 0.05),
    upper90  = quantile(slope_log10, 0.95),
    .groups = "drop"
  )

ltA_fig_slopes <- ggplot(ltA_slope_summ, aes(x = n_year, y = estimate, colour = year_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.6) +
  scale_colour_viridis_d(option = "viridis", guide = "none") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  labs(
    x = "Monitoring year",
    y = "Slope: d(log10(larvae))/d(°C)",
    # title = "Temperature sensitivity through time (median salinity)",
    subtitle = "c)"
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank())

ltA_fig_slopes

# Optional combined view (patchwork)
ltA_W_fig_spag_med + ltA_fig_intercepts + ltA_fig_slopes

ltA_W_fig_spag_med +
    ltA_fig_intercepts +
    ltA_fig_slopes +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")


# ============================================================
# ltA_T) YEAR ON X × TEMP LEVELS @ MEDIAN SALINITY (aligned)
# ============================================================

ltA_T_temp_df <- l_total_dat %>%
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


ltA_T_year_seq <- seq(
  floor(min(l_total_dat$n_year, na.rm = TRUE)),
  ceiling(max(l_total_dat$n_year, na.rm = TRUE)),
  by = 1
)

ltA_T_grid <- crossing(
  n_year = ltA_T_year_seq,
  ltA_T_temp_df,
  tibble(salinity.m = ltA_ref$sal_m_q50, sal_label = factor("Median salinity"))
) %>%
  mutate(
    n_year.m      = n_year     - ltA_ref$yr_mean,
    water_temp.m  = water_temp - ltA_ref$wt_mean,
    julian_date.m = 0,
    row_id        = row_number()
  )

ltA_T_ep_mat <- posterior_epred(
  ltotal_temp_time_sal,
  newdata    = ltA_T_grid,
  re_formula = NA
)

ltA_T_pool_size <- min(2000, nrow(ltA_T_ep_mat))
ltA_T_draw_ids  <- sample(seq_len(nrow(ltA_T_ep_mat)), size = ltA_T_pool_size)

ltA_T_ep_long <- ltA_T_ep_mat[ltA_T_draw_ids, ] %>%
  as_tibble() %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(cols = - .draw, names_to = "row_id", values_to = "epred") %>%
  mutate(row_id = readr::parse_number(row_id)) %>%
  left_join(
    ltA_T_grid %>% select(row_id, n_year, temp_label, sal_label),
    by = "row_id"
  )

# unfiltered summary ribbons (50% + 90%)
ltA_T_summ <- ltA_T_ep_long %>%
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
ltA_T_y_max_obs  <- max(l_total_dat$larvae_total, na.rm = TRUE)
ltA_T_y_min_plot <- max(1, min(ltA_T_summ$low90, na.rm = TRUE))
ltA_T_y_max_plot <- min(ltA_T_y_max_obs * 1.1, max(ltA_T_summ$up90, na.rm = TRUE))

ltA_T_fig_time_temp <- ggplot(
  ltA_T_summ,
  aes(x = n_year, y = med, colour = temp_label, group = temp_label)
) +
  geom_ribbon(aes(ymin = low90, ymax = up90, fill = temp_label), alpha = 0.10, colour = NA) +
  geom_ribbon(aes(ymin = low50, ymax = up50, fill = temp_label), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.1) +
 # facet_wrap(~ sal_label, nrow = 1) +
 # scale_colour_viridis_d(option = "viridis") +
 #scale_fill_viridis_d(option = "viridis") +
  # --- Dark2 palette with explicit mapping ---
  scale_colour_manual(name = "Temp",
    values = c( "Cool temp (25th pct)" = "#1B9E77",  # green (Dark2)
                                  "Median temp"          = "#7570B3",  # purple (Dark2)
                                  "Warm temp (75th pct)" = "#D95F02"   # orange (Dark2)
  ) ) +
  scale_fill_manual(values = c( "Cool temp (25th pct)" = "#1B9E77",  # green (Dark2)
                                  "Median temp"          = "#7570B3",  # purple (Dark2)
                                  "Warm temp (75th pct)" = "#D95F02"   # orange (Dark2)
  ) ) +
  guides(fill = "none") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(4, 8, 16, 32, 64, 128, 256, 512, 1024),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(ltA_T_y_min_plot, ltA_T_y_max_plot)) +
  labs(
    x = "Monitoring year",
    y = "Total oyster larvae",
    subtitle= "d)", colour= "Temp", fill= "Temp"
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

ltA_T_fig_time_temp


# ============================================================
# ltA_T) INTERCEPTS + SLOPES for COOL / MEDIAN / WARM temps
# - uses the SAME ltA_T_ep_long object from the year-on-x workflow
# - tidyverse only + no custom functions + no if/else
# - uncertainty: 50% + 90% CrI everywhere
# ============================================================

# ------------------------------------------------------------
# 0) Use your existing temp levels table (must have: water_temp, temp_label)
#     and your year sequence grid already used for ltA_T_ep_long
# ------------------------------------------------------------
# Assumes these already exist from the aligned ltA_T workflow:
#   ltA_T_temp_df, ltA_T_summ, ltA_T_ep_long, ltA_ref, ltA_T_year_seq

# ------------------------------------------------------------
# 1) INTERCEPTS: predicted larvae at a reference year (mean year)
#    - within each draw × temp_label, take epred at year closest to ref
# ------------------------------------------------------------
ltA_T_ref_year <- l_total_dat %>%
  ungroup() %>%
  summarise(ref_year = mean(n_year, na.rm = TRUE)) %>%
  pull(ref_year)

ltA_T_intercept_draws <- ltA_T_ep_long %>%
  group_by(.draw, temp_label) %>%
  mutate(dist_to_ref_year = abs(n_year - ltA_T_ref_year)) %>%
  slice_min(order_by = dist_to_ref_year, n = 1, with_ties = FALSE) %>%
  ungroup()

ltA_T_intercept_summ <- ltA_T_intercept_draws %>%
  group_by(temp_label) %>%
  summarise(
    estimate = median(epred),
    lower50  = quantile(epred, 0.25),
    upper50  = quantile(epred, 0.75),
    lower90  = quantile(epred, 0.05),
    upper90  = quantile(epred, 0.95),
    .groups  = "drop"
  )

ltA_T_fig_intercepts <- ggplot(
  ltA_T_intercept_summ,
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
    subtitle = "e)"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x  = element_blank(),
    legend.position = "none"
  )

ltA_T_fig_intercepts


# ------------------------------------------------------------
# 2) SLOPES: change over time within each temp level
#    - compute slope within each draw × temp_label using cov/var
#    - (optional) do on log10 scale for stability/interpretability (recommended)
# ------------------------------------------------------------
ltA_T_slope_draws <- ltA_T_ep_long %>%
  mutate(log10_epred = log10(epred)) %>%
  group_by(.draw, temp_label) %>%
  summarise(
    slope_log10 = cov(n_year, log10_epred) / var(n_year),
    .groups = "drop"
  )

ltA_T_slope_summ <- ltA_T_slope_draws %>%
  group_by(temp_label) %>%
  summarise(
    estimate = median(slope_log10),
    lower50  = quantile(slope_log10, 0.25),
    upper50  = quantile(slope_log10, 0.75),
    lower90  = quantile(slope_log10, 0.05),
    upper90  = quantile(slope_log10, 0.95),
    .groups  = "drop"
  )

view(ltA_T_slope_summ)

ltA_T_fig_slopes <- ggplot(
  ltA_T_slope_summ,
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
    subtitle = "f)"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x  = element_blank(),
    legend.position = "none"
  )

ltA_T_fig_slopes


# Setting the legend from graph a) to be centered below the 3 graphs patchworked together
ltA_W_fig_spag_med + ltA_fig_intercepts + ltA_fig_slopes

final_top_row <- (ltA_W_fig_spag_med +
                    ltA_fig_intercepts +
                    ltA_fig_slopes) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom") 

final_top_row

# Setting the legend from graph d) to ve centered below the 3 graphs patchworked together
ltA_T_fig_time_temp + ltA_T_fig_intercepts + ltA_T_fig_slopes
 
 ltA_T_fig_intercepts <- ltA_T_fig_intercepts +
   guides(colour = "none", fill = "none") +
   theme(legend.position = "none")
 
 ltA_T_fig_slopes <- ltA_T_fig_slopes +
   guides(colour = "none", fill = "none") +
   theme(legend.position = "none")
 
 final_bottom_row <- (ltA_T_fig_time_temp +
   ltA_T_fig_intercepts +
   ltA_T_fig_slopes) +
   plot_layout(ncol = 3, guides = "collect") &
   theme(
     legend.position = "bottom",
     legend.justification = "center")
 
 final_bottom_row


 # Final combined figures
 final_top_row / final_bottom_row
 
 

