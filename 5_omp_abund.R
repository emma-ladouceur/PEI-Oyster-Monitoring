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
# UNIQUE-NAME WORKFLOW: TOTAL LARVAE FIGURES (NO NAME COLLISIONS)
# Prefix convention:
# ============================================================
# ============================================================
# ltW) WATER TEMP ON X-AXIS — 3 SALINITY PANELS (median + 90% CrI)
# FIX: join by row_id 
# ============================================================


# -----------------------------
# ltW_1) Salinity slices (centered) using quantiles (25/50/75)
# -----------------------------
ltW_sal_q <- quantile(l_total_dat$salinity.m, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

ltW_sal_df <- tibble(
  salinity.m = as.numeric(ltW_sal_q),
  sal_label  = factor(
    c("Low salinity (25th pct)", "Median salinity", "High salinity (75th pct)"),
    levels = c("Low salinity (25th pct)", "Median salinity", "High salinity (75th pct)")
  )
)

# -----------------------------
# ltW_2) Prediction grid (year-specific temp support) + row_id FIRST
# -----------------------------
ltW_pred_grid <- l_total_dat %>%
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
  crossing(ltW_sal_df) %>%
  mutate(
    # center exactly as in model
    water_temp.m  = water_temp - mean(l_total_dat$water_temp, na.rm = TRUE),
    n_year.m      = n_year     - mean(l_total_dat$n_year,     na.rm = TRUE),
    
    # group/colour by year
    year_group    = factor(n_year),
    
    # required by (1 | julian_date.m)
    julian_date.m = 0,
    
    # IMPORTANT: stable id for joining predictions back to rows
    row_id        = row_number()
  )

# -----------------------------
# ltW_3) Posterior expected values (population-level)
# -----------------------------
ltW_ep_mat <- posterior_epred(
  ltotal_temp_time_sal,
  newdata    = ltW_pred_grid,
  re_formula = NA
)

# -----------------------------
# ltW_4) Draw pool BEFORE pivoting long (memory-safe)
# -----------------------------
set.seed(123)
ltW_pool_size <- min(2000, nrow(ltW_ep_mat))
ltW_draw_ids  <- sample(seq_len(nrow(ltW_ep_mat)), size = ltW_pool_size)

# -----------------------------
# ltW_5) Draws × rows -> long, then JOIN grid metadata by row_id
# -----------------------------
ltW_ep_long <- ltW_ep_mat[ltW_draw_ids, ] %>%
  as_tibble() %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols      = - .draw,
    names_to  = "row_id",
    values_to = "epred"
  ) %>%
  mutate(row_id = readr::parse_number(row_id)) %>%
  left_join(
    ltW_pred_grid %>% select(row_id, n_year, water_temp, year_group, sal_label),
    by = "row_id"
  )

# -----------------------------
# ltW_6) Summarise: median + 90% CrI at each point
# -----------------------------
ltW_summ <- ltW_ep_long %>%
  group_by(n_year, water_temp, year_group, sal_label) %>%
  summarise(
    estimate = median(epred),
    lower90  = quantile(epred, 0.05),
    upper90  = quantile(epred, 0.95),
    .groups  = "drop"
  )

# -----------------------------
# ltW_7) Log-safe crop
# -----------------------------
ltW_y_max_obs  <- max(l_total_dat$larvae_total, na.rm = TRUE)
ltW_y_min_plot <- max(1, min(ltW_summ$lower90, na.rm = TRUE))
ltW_y_max_plot <- min(ltW_y_max_obs * 1.1, max(ltW_summ$upper90, na.rm = TRUE))

# -----------------------------
# ltW_8) Figure
# -----------------------------
ltW_fig_panel <- ggplot(
  ltW_summ,
  aes(x = water_temp, y = estimate, group = year_group, colour = year_group)
) +
  geom_ribbon(aes(ymin = lower90, ymax = upper90), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ sal_label, nrow = 1) +
  scale_colour_viridis_d(option = "viridis") +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30, 35)) +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(24, 48, 96, 192),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(ltW_y_min_plot, ltW_y_max_plot)) +
  labs(
    x = "Surface water temperature (°C)",
    y = "Total oyster larvae",
    title = "Temperature–abundance relationships vary through time and with salinity",
    subtitle = "Population-level medians; shaded bands show 90% credible intervals"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ltW_fig_panel

# ============================================================
# ltW) MEDIAN-SALINITY ONLY — SPAGHETTI FIGURE (median + draws)
# - thin lines: posterior draws (filtered to central 80% region)
# - thick lines: posterior medians by year
# - tidyverse only workflow (no bind_cols after pivot_longer)
# ============================================================


# -----------------------------
# ltW_S1) Median salinity only (centered scale used in model)
# -----------------------------
ltW_sal_med <- quantile(l_total_dat$salinity.m, probs = 0.50, na.rm = TRUE)

ltW_sal_df_med <- tibble(
  salinity.m = as.numeric(ltW_sal_med),
  sal_label  = factor("Median salinity", levels = "Median salinity")
)

# -----------------------------
# ltW_S2) Prediction grid (year-specific temp support) + row_id FIRST
# -----------------------------
ltW_pred_grid_med <- l_total_dat %>%
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
  crossing(ltW_sal_df_med) %>%
  mutate(
    water_temp.m  = water_temp - mean(l_total_dat$water_temp, na.rm = TRUE),
    n_year.m      = n_year     - mean(l_total_dat$n_year,     na.rm = TRUE),
    year_group    = factor(n_year),
    julian_date.m = 0,
    row_id        = row_number()
  )

# -----------------------------
# ltW_S3) Posterior expected values (population-level)
# -----------------------------
ltW_ep_mat_med <- posterior_epred(
  ltotal_temp_time_sal,
  newdata    = ltW_pred_grid_med,
  re_formula = NA
)

# -----------------------------
# ltW_S4) Sample a draw pool BEFORE pivoting long (memory-safe)
# -----------------------------
set.seed(123)
ltW_pool_size_med <- min(2000, nrow(ltW_ep_mat_med))
ltW_draw_ids_med  <- sample(seq_len(nrow(ltW_ep_mat_med)), size = ltW_pool_size_med)

# -----------------------------
# ltW_S5) Draws × rows -> long, then JOIN grid metadata by row_id
# -----------------------------
ltW_ep_long_med <- ltW_ep_mat_med[ltW_draw_ids_med, ] %>%
  as_tibble() %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols      = - .draw,
    names_to  = "row_id",
    values_to = "epred"
  ) %>%
  mutate(row_id = readr::parse_number(row_id)) %>%
  left_join(
    ltW_pred_grid_med %>% select(row_id, n_year, water_temp, year_group, sal_label),
    by = "row_id"
  )

# -----------------------------
# ltW_S6) Thick lines: median + 90% CrI (for reference + cropping)
# -----------------------------
ltW_summ_med <- ltW_ep_long_med %>%
  group_by(row_id, n_year, water_temp, year_group, sal_label) %>%
  summarise(
    estimate = median(epred),
    lower90  = quantile(epred, 0.05),
    upper90  = quantile(epred, 0.95),
    .groups  = "drop"
  )

# -----------------------------
# ltW_S7) Log-safe crop (keeps plot from exploding)
# -----------------------------
ltW_y_max_obs_med  <- max(l_total_dat$larvae_total, na.rm = TRUE)
ltW_y_min_plot_med <- max(1, min(ltW_summ_med$lower90, na.rm = TRUE))
ltW_y_max_plot_med <- min(ltW_y_max_obs_med * 1.1, max(ltW_summ_med$upper90, na.rm = TRUE))

# -----------------------------
# ltW_S8) Central 80% envelope at each x-point (row_id)
# -----------------------------
ltW_bounds80_med <- ltW_ep_long_med %>%
  group_by(row_id) %>%
  summarise(
    q10 = quantile(epred, 0.10),
    q90 = quantile(epred, 0.90),
    .groups = "drop"
  )

# -----------------------------
# ltW_S9) Keep draws that mostly fall inside the envelope (>=95%)
# -----------------------------
ltW_draw_keep_med <- ltW_ep_long_med %>%
  left_join(ltW_bounds80_med, by = "row_id") %>%
  mutate(in80 = epred >= q10 & epred <= q90) %>%
  group_by(.draw) %>%
  summarise(prop_in80 = mean(in80), .groups = "drop") %>%
  filter(prop_in80 >= 0.95)

# -----------------------------
# ltW_S10) Sample spaghetti draws (up to 150) WITHOUT using n() inside slice_sample()
# -----------------------------
ltW_n_keep_med <- ltW_draw_keep_med %>% summarise(n = n()) %>% pull(n)
ltW_n_spag_med <- min(150, ltW_n_keep_med)

ltW_spag_ids_med <- ltW_draw_keep_med %>%
  slice_sample(n = ltW_n_spag_med) %>%
  pull(.draw)

ltW_ep_spag_med <- ltW_ep_long_med %>%
  filter(.draw %in% ltW_spag_ids_med)

# -----------------------------
# ltW_S11) Spaghetti figure (median salinity only)
# -----------------------------
ltW_fig_spag_med <- ggplot() +
  geom_line(
    data = ltW_ep_spag_med,
    aes(
      x = water_temp,
      y = epred,
      group = interaction(.draw, year_group),
      colour = year_group
    ),
    linewidth = 0.4,
    alpha = 0.15
  ) +
  geom_line(
    data = ltW_summ_med,
    aes(
      x = water_temp,
      y = estimate,
      group = year_group,
      colour = year_group
    ),
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
  coord_cartesian(ylim = c(ltW_y_min_plot_med, ltW_y_max_plot_med)) +
  labs(
    x = "Surface water temperature (°C)",
    y = "Total oyster larvae",
    title = "Temperature–abundance relationships at median salinity",
    subtitle = "Thin lines = filtered posterior draws; thick lines = posterior medians"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ltW_fig_spag_med


# ============================================================
# ltW) INTERCEPT + SLOPE FIGURES (median salinity)
# Uses ltW_ep_long_med (draws in long form) from the spaghetti workflow
# ============================================================


# ------------------------------------------------------------
# 0) Choose a reference temperature for "intercept"
#    - Use the median observed temperature in the original dataset
#    - This is interpretable and always supported by data
# ------------------------------------------------------------
ltW_temp_ref <- l_total_dat %>%
  summarise(temp_ref = median(water_temp, na.rm = TRUE)) %>%
  pull(temp_ref)

# ------------------------------------------------------------
# 1) INTERCEPT FIGURE DATA
#    For each year × draw:
#      - find the prediction point whose water_temp is closest to temp_ref
#      - take that epred as the "intercept-like" value for that year
#
# Notes:
# - We do this WITHIN draw → uncertainty is preserved
# - slice_min() avoids any if/else logic
# ------------------------------------------------------------
ltW_intercept_draws <- ltW_ep_long_med %>%
  group_by(.draw, n_year, year_group) %>%
  mutate(dist_to_ref = abs(water_temp - ltW_temp_ref)) %>%
  slice_min(order_by = dist_to_ref, n = 1, with_ties = FALSE) %>%
  ungroup()

# Summarise across draws for each year
ltW_intercept_summ <- ltW_intercept_draws %>%
  group_by(n_year, year_group) %>%
  summarise(
    estimate = median(epred),
    lower50  = quantile(epred, 0.25),
    upper50  = quantile(epred, 0.75),
    lower90  = quantile(epred, 0.05),
    upper90  = quantile(epred, 0.95),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 2) INTERCEPT FIGURE
# ------------------------------------------------------------
ltW_fig_intercepts <- ggplot(
  ltW_intercept_summ,
  aes(x = n_year, y = estimate)
) +
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
      round(ltW_temp_ref, 1), " °C). Points = posterior medians; thick bars = 50% CrI; thin bars = 90% CrI."
    )
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank())

ltW_fig_intercepts


# ------------------------------------------------------------
# 3) SLOPE FIGURE DATA (temperature sensitivity by year)
#    Goal: within each year, estimate how larvae changes with temperature.
#
# Because response is lognormal (positive, skewed), slopes on a log scale
# are more stable and interpretable:
#   slope = d(log10(larvae))/d(temp)
# Interpretation:
#   +0.10 per °C ≈ 10^0.10 = 1.26× per °C (~26% increase per °C)
#
# Compute slope WITHIN each draw × year using cov/var (no lm, no apply)
# ------------------------------------------------------------
ltW_slope_draws <- ltW_ep_long_med %>%
  mutate(log10_epred = log10(epred)) %>%
  group_by(.draw, n_year, year_group) %>%
  summarise(
    slope_log10 = cov(water_temp, log10_epred) / var(water_temp),
    .groups = "drop"
  )

# Summarise slopes across draws for each year
ltW_slope_summ <- ltW_slope_draws %>%
  group_by(n_year, year_group) %>%
  summarise(
    estimate = median(slope_log10),
    lower50  = quantile(slope_log10, 0.25),
    upper50  = quantile(slope_log10, 0.75),
    lower90  = quantile(slope_log10, 0.05),
    upper90  = quantile(slope_log10, 0.95),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 4) SLOPE FIGURE
# ------------------------------------------------------------
ltW_fig_slopes <- ggplot(
  ltW_slope_summ,
  aes(x = n_year, y = estimate)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.6) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  labs(
    x = "Monitoring year",
    y = "Slope: d(log10(larvae))/d(°C)",
    title = "Temperature sensitivity through time (median salinity)",
    subtitle = "Within-year slopes. Points = posterior medians; thick bars = 50% CrI; thin bars = 90% CrI."
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank())

ltW_fig_slopes

ltW_fig_spag_med + ltW_fig_intercepts + ltW_fig_slopes

# ------------------------------------------------------------
# 5) Optional: view side-by-side if patchwork is loaded
# ------------------------------------------------------------
# ltW_fig_intercepts + ltW_fig_slopes


# ============================================================
# ltT) YEAR ON X-AXIS × TEMP LEVELS (cool/med/warm) @ MEDIAN SALINITY
# Includes:
#   - non-spaghetti (median + 90% CrI)
#   - spaghetti (filtered draws)
#   - slope figure matched to spaghetti
# ============================================================

# -----------------------------
# ltT_1) Median salinity only
# -----------------------------
ltT_sal_med <- quantile(l_total_dat$salinity.m, probs = 0.5, na.rm = TRUE)

ltT_sal_df <- tibble(
  salinity.m = as.numeric(ltT_sal_med),
  sal_label  = factor("Median salinity")
)

# -----------------------------
# ltT_2) Temperature levels (cool/med/warm) on ORIGINAL scale
# -----------------------------
ltT_temp_q <- quantile(l_total_dat$water_temp, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

ltT_temp_df <- tibble(
  water_temp = as.numeric(ltT_temp_q),
  temp_label = factor(
    c("Cool temp (25th pct)", "Median temp", "Warm temp (75th pct)"),
    levels = c("Cool temp (25th pct)", "Median temp", "Warm temp (75th pct)")
  )
)

# -----------------------------
# ltT_3) Year sequence (original scale)
# -----------------------------
ltT_year_seq <- seq(
  floor(min(l_total_dat$n_year, na.rm = TRUE)),
  ceiling(max(l_total_dat$n_year, na.rm = TRUE)),
  by = 1
)

# -----------------------------
# ltT_4) Prediction grid FROM REAL group levels
# (keeps bay/location_clean levels from the dataset)
# -----------------------------
ltT_pred_grid <- l_total_dat %>%
  distinct(bay, location_clean) %>%
  crossing(
    n_year = ltT_year_seq,
    ltT_temp_df,
    ltT_sal_df
  ) %>%
  mutate(
    # center exactly as in model
    n_year.m     = n_year     - mean(l_total_dat$n_year,     na.rm = TRUE),
    water_temp.m = water_temp - mean(l_total_dat$water_temp, na.rm = TRUE),
    
    # required by (1 | julian_date.m)
    julian_date.m = 0,
    
    # stable row id for joins
    row_id = row_number()
  )

# -----------------------------
# ltT_5) Posterior expected values (population-level)
# -----------------------------
ltT_ep_mat <- posterior_epred(
  ltotal_temp_time_sal,
  newdata    = ltT_pred_grid,
  re_formula = NA
)

# -----------------------------
# ltT_6) Draw pool BEFORE going long
# -----------------------------
set.seed(123)
ltT_pool_size <- min(2000, nrow(ltT_ep_mat))
ltT_draw_ids  <- sample(seq_len(nrow(ltT_ep_mat)), size = ltT_pool_size)

ltT_ep_long <- ltT_ep_mat[ltT_draw_ids, ] %>%
  as_tibble() %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols = - .draw,
    names_to = "row_id",
    values_to = "epred"
  ) %>%
  mutate(row_id = readr::parse_number(row_id)) %>%
  left_join(
    ltT_pred_grid %>% select(row_id, n_year, temp_label, sal_label),
    by = "row_id"
  )

# -----------------------------
# ltT_7) Non-spaghetti summary (median + 90% CrI across draws)
# -----------------------------
ltT_summ_time <- ltT_ep_long %>%
  group_by(n_year, temp_label, sal_label) %>%
  summarise(
    estimate = median(epred),
    lower90  = quantile(epred, 0.05),
    upper90  = quantile(epred, 0.95),
    .groups  = "drop"
  )

# -----------------------------
# ltT_8) Log-safe crop for time-temp figs
# -----------------------------
ltT_y_max_obs  <- max(l_total_dat$larvae_total, na.rm = TRUE)
ltT_y_min_plot <- max(1, min(ltT_summ_time$lower90, na.rm = TRUE))
ltT_y_max_plot <- min(ltT_y_max_obs * 1.1, max(ltT_summ_time$upper90, na.rm = TRUE))

# -----------------------------
# ltT_9) Non-spaghetti figure (year on x, temp levels as lines)
# -----------------------------
ltT_fig_time_temp <- ggplot(
  ltT_summ_time,
  aes(x = n_year, y = estimate, colour = temp_label, group = temp_label)
) +
  geom_ribbon(
    aes(ymin = lower90, ymax = upper90, fill = temp_label),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(linewidth = 1.0) +
  facet_wrap(~ sal_label, nrow = 1) +
  scale_colour_brewer(name = "Surface Water\nTemperature (°C)", type = "qual", palette = "Dark2") +
  scale_fill_brewer(name = "Surface Water\nTemperature (°C)", type = "qual", palette = "Dark2") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(24, 48, 96, 192),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(ltT_y_min_plot, ltT_y_max_plot)) +
  labs(
    x = "Monitoring year",
    y = "Total oyster larvae",
    title = "Temporal trends vary with temperature (median salinity)",
    subtitle = "Population-level medians at cool/median/warm temperatures; ribbons show 90% credible intervals"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ltT_fig_time_temp


# ============================================================
# ltT SPAGHETTI VERSION (filtered draws)
# ============================================================

# -----------------------------
# ltT_10) Envelope (10–90%) at each (year × temp level)
# -----------------------------
ltT_bounds80 <- ltT_ep_long %>%
  group_by(n_year, temp_label) %>%
  summarise(
    q10 = quantile(epred, 0.10),
    q90 = quantile(epred, 0.90),
    .groups = "drop"
  )

# -----------------------------
# ltT_11) Keep "well-behaved" draws (>=95% inside envelope)
# -----------------------------
ltT_draw_keep <- ltT_ep_long %>%
  left_join(ltT_bounds80, by = c("n_year", "temp_label")) %>%
  mutate(in80 = epred >= q10 & epred <= q90) %>%
  group_by(.draw) %>%
  summarise(prop_in80 = mean(in80), .groups = "drop") %>%
  filter(prop_in80 >= 0.95)

# -----------------------------
# ltT_12) Sample spaghetti draw IDs (up to 150) WITHOUT n() inside slice_sample()
# -----------------------------
ltT_n_keep <- ltT_draw_keep %>% summarise(n = n()) %>% pull(n)
ltT_n_spag <- min(150, ltT_n_keep)

ltT_spag_ids <- ltT_draw_keep %>%
  slice_sample(n = ltT_n_spag) %>%
  pull(.draw)

ltT_ep_spag <- ltT_ep_long %>%
  filter(.draw %in% ltT_spag_ids)

# -----------------------------
# ltT_13) Spaghetti figure (thin draws + thick medians)
# -----------------------------
ltT_fig_time_temp_spag <- ggplot() +
  geom_line(
    data = ltT_ep_spag,
    aes(
      x = n_year,
      y = epred,
      group = interaction(.draw, temp_label),
      colour = temp_label
    ),
    linewidth = 0.4,
    alpha = 0.15
  ) +
  geom_line(
    data = ltT_summ_time,
    aes(x = n_year, y = estimate, group = temp_label, colour = temp_label),
    linewidth = 1.1
  ) +
  facet_wrap(~ sal_label, nrow = 1) +
  scale_colour_brewer(name = "Surface Water\nTemperature (°C)", type = "qual", palette = "Dark2") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(4, 8, 16, 32, 64, 128, 256, 512, 1024),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(ltT_y_min_plot, ltT_y_max_plot)) +
  labs(
    x = "Monitoring year",
    y = "Total oyster larvae",
    title = "Temporal trends vary with temperature (median salinity)",
    subtitle = "(a) Thin lines = filtered posterior draws; thick lines = posterior medians"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ltT_fig_time_temp_spag


# ============================================================
# ltT SLOPE FIGURE (matched to spaghetti)
# ============================================================

# -----------------------------
# ltT_14) Slopes within each draw × temp level
# slope = cov(x,y)/var(x) (simple linear slope, draw-by-draw)
# -----------------------------
ltT_slope_draws <- ltT_ep_spag %>%
  group_by(.draw, temp_label) %>%
  summarise(
    slope = cov(n_year, epred) / var(n_year),
    .groups = "drop"
  )

# -----------------------------
# ltT_15) Summarise slopes across draws (mean + 50%/90% CrI)
# -----------------------------
ltT_slope_summ <- ltT_slope_draws %>%
  group_by(temp_label) %>%
  summarise(
    estimate = mean(slope),
    lower50  = quantile(slope, 0.25),
    upper50  = quantile(slope, 0.75),
    lower90  = quantile(slope, 0.05),
    upper90  = quantile(slope, 0.95),
    .groups  = "drop"
  )

# -----------------------------
# ltT_16) Slope figure
# -----------------------------
ltT_fig_time_temp_slopes <- ggplot(
  ltT_slope_summ,
  aes(x = temp_label, y = estimate, colour = temp_label)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.9, alpha = 0.6) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.2) +
  geom_point(size = 3) +
  scale_colour_brewer(name = "Surface Water\nTemperature (°C)", type = "qual", palette = "Dark2") +
  coord_flip() +
  labs(
    x = NULL,
    y = "Slope (change in total larvae per year)",
    title = "Temporal trends strengthen with temperature (median salinity)",
    subtitle = "(b) Points = posterior mean slopes; thick bars = 50% CrI; thin bars = 90% CrI"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

ltT_fig_time_temp_slopes




# Optional: show spaghetti + slopes together (requires patchwork loaded)
 ltT_fig_time_temp_spag + ltT_fig_time_temp_slopes


 
 
 
 # ============================================================
 # ltotal: Ribbon fig + slopes + intercepts (year × temp; salinity held at mean)
 # For: ltotal_temp_time_sal (data = l_total_dat)
 # - no custom functions
 # - tidyverse + tidybayes
 # - ribbons = posterior uncertainty (80% + 50%)
 # - family is lognormal(): epred is on the RESPONSE scale (positive)
 # ============================================================
 
 set.seed(123)
 
 # ------------------------------------------------------------
 # lt_01) References for centering + temp quantiles (UNCENTERED vars)
 # ------------------------------------------------------------
 lt_alt_ref <- l_total_dat %>%
   ungroup() %>%
   summarise(
     n_year_mean     = mean(n_year, na.rm = TRUE),
     water_temp_mean = mean(water_temp, na.rm = TRUE),
     sal_mean        = mean(salinity, na.rm = TRUE),
     temp_cool       = quantile(water_temp, 0.10, na.rm = TRUE),
     temp_med        = quantile(water_temp, 0.50, na.rm = TRUE),
     temp_warm       = quantile(water_temp, 0.90, na.rm = TRUE),
     julian_m_mean   = mean(julian_date.m, na.rm = TRUE) # for random effect (1 | julian_date.m)
   )
 
 # ------------------------------------------------------------
 # lt_02) Temperature level lookup table
 # ------------------------------------------------------------
 lt_alt_temp_levels <- tibble(
   temp_level = factor(c("Cool", "Median", "Warm"),
                       levels = c("Cool", "Median", "Warm")),
   water_temp = c(lt_alt_ref$temp_cool,
                  lt_alt_ref$temp_med,
                  lt_alt_ref$temp_warm)
 )
 
 # ------------------------------------------------------------
 # lt_03) Prediction grid (year sequence × temp levels),
 #        salinity fixed at mean, julian_date.m fixed at mean
 # ------------------------------------------------------------
 lt_alt_newdat_year_temp <- crossing(
   n_year = seq(min(l_total_dat$n_year, na.rm = TRUE),
                max(l_total_dat$n_year, na.rm = TRUE),
                length.out = 200),
   temp_level = lt_alt_temp_levels$temp_level
 ) %>%
   left_join(lt_alt_temp_levels, by = "temp_level") %>%
   mutate(
     salinity       = lt_alt_ref$sal_mean,
     julian_date.m  = lt_alt_ref$julian_m_mean,     # holds the julian RE at "typical" value
     salinity.m     = salinity   - lt_alt_ref$sal_mean,
     n_year.m       = n_year     - lt_alt_ref$n_year_mean,
     water_temp.m   = water_temp - lt_alt_ref$water_temp_mean
   )
 
 # ------------------------------------------------------------
 # lt_04) Posterior expected draws (population-level)
 #   - IMPORTANT: we drop ALL random effects here (re_formula = NA)
 #     so (1 | julian_date.m) is not included in the prediction.
 #   - If you *want* the julian_date random effect included, set re_formula = NULL
 #     AND ensure newdata supplies valid julian_date.m levels.
 # ------------------------------------------------------------
 lt_alt_draws <- ltotal_temp_time_sal %>%
   add_epred_draws(newdata = lt_alt_newdat_year_temp, re_formula = NA) %>%
   ungroup() %>%
   rename(lt_alt_epred = .epred)
 
 # ------------------------------------------------------------
 # lt_05) Summaries for ribbons: 80% outer + 50% inner + median line
 # ------------------------------------------------------------
 lt_alt_summ80 <- lt_alt_draws %>%
   group_by(n_year, temp_level) %>%
   median_qi(lt_alt_epred, .width = 0.80) %>%
   ungroup() %>%
   rename(
     lt_alt_med   = lt_alt_epred,
     lt_alt_low80 = .lower,
     lt_alt_up80  = .upper
   )
 
 lt_alt_summ50 <- lt_alt_draws %>%
   group_by(n_year, temp_level) %>%
   median_qi(lt_alt_epred, .width = 0.50) %>%
   ungroup() %>%
   rename(
     lt_alt_low50 = .lower,
     lt_alt_up50  = .upper
   )
 
 lt_alt_ribbon_dat <- lt_alt_summ80 %>%
   left_join(
     lt_alt_summ50 %>%
       select(n_year, temp_level, lt_alt_low50, lt_alt_up50),
     by = c("n_year", "temp_level")
   )
 
 # ------------------------------------------------------------
 # lt_06) Ribbon figure
 # ------------------------------------------------------------
 lt_alt_fig_year3temp_ribbons <- ggplot(lt_alt_ribbon_dat) +
   geom_ribbon(
     aes(x = n_year, ymin = lt_alt_low80, ymax = lt_alt_up80, fill = temp_level),
     alpha = 0.18, colour = NA
   ) +
   geom_ribbon(
     aes(x = n_year, ymin = lt_alt_low50, ymax = lt_alt_up50, fill = temp_level),
     alpha = 0.32, colour = NA
   ) +
   geom_line(
     aes(x = n_year, y = lt_alt_med, colour = temp_level),
     linewidth = 1.0
   ) +
   scale_colour_viridis_d(option = "viridis") +
   scale_fill_viridis_d(option = "viridis") +
   coord_cartesian(
     ylim = quantile(
       lt_alt_draws$lt_alt_epred,
       probs = c(0.01, 0.99),
       na.rm = TRUE
     )
   ) +
   labs(
     x = "Year",
     y = "Predicted total larvae (response scale)",
     title = "Temporal trends in larvae abundance by temperature (salinity held at mean)",
     subtitle = "Lines = posterior medians; ribbons = 50% (inner) and 80% (outer) CrI"
   ) +
   theme_bw(base_size = 11) +
   theme(panel.grid.minor = element_blank())
 
 lt_alt_fig_year3temp_ribbons
 
 
 # ============================================================
 # lt: SLOPE FIGURE (matched to ribbons)
 # - draw-by-draw slopes within temp_level
 # - summarises across draws (mean + 50%/90% CrI)
 # ============================================================
 
 lt_slope_draws <- lt_alt_draws %>%
   group_by(.draw, temp_level) %>%
   summarise(
     slope = cov(n_year, lt_alt_epred) / var(n_year),
     .groups = "drop"
   )
 
 lt_slope_summ <- lt_slope_draws %>%
   group_by(temp_level) %>%
   summarise(
     estimate = mean(slope),
     lower50  = quantile(slope, 0.25),
     upper50  = quantile(slope, 0.75),
     lower90  = quantile(slope, 0.05),
     upper90  = quantile(slope, 0.95),
     .groups  = "drop"
   )
 
 lt_fig_year3temp_slopes <- ggplot(
   lt_slope_summ,
   aes(x = temp_level, y = estimate, colour = temp_level)
 ) +
   geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
   geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.9, alpha = 0.6) +
   geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.2) +
   geom_point(size = 3) +
   scale_colour_viridis_d(option = "viridis") +
   coord_flip() +
   labs(
     x = NULL,
     y = "Slope (change in predicted larvae per year)",
     title = "Temporal trends by temperature (salinity held at mean)",
     subtitle = "(b) Points = posterior mean slopes; thick bars = 50% CrI; thin bars = 90% CrI"
   ) +
   theme_bw(base_size = 18) +
   theme(
     panel.grid.minor = element_blank(),
     legend.position = "none"
   )
 
 lt_fig_year3temp_slopes
 
 
 # ============================================================
 # lt: INTERCEPT FIGURE (matched to ribbons + slopes style)
 # - intercept = predicted value at reference year (mean year)
 # - summarises across draws (mean + 50%/90% CrI)
 # ============================================================
 
 lt_ref_year <- lt_alt_ref$n_year_mean
 
 lt_intercept_draws <- lt_alt_draws %>%
   mutate(.dist = abs(n_year - lt_ref_year)) %>%
   group_by(.draw, temp_level) %>%
   slice_min(.dist, n = 1, with_ties = FALSE) %>%
   ungroup() %>%
   transmute(
     .draw,
     temp_level,
     intercept = lt_alt_epred
   )
 
 lt_intercept_summ <- lt_intercept_draws %>%
   group_by(temp_level) %>%
   summarise(
     estimate = mean(intercept),
     lower50  = quantile(intercept, 0.25),
     upper50  = quantile(intercept, 0.75),
     lower90  = quantile(intercept, 0.05),
     upper90  = quantile(intercept, 0.95),
     .groups  = "drop"
   )
 
 lt_fig_year3temp_intercepts <- ggplot(
   lt_intercept_summ,
   aes(x = temp_level, y = estimate, colour = temp_level)
 ) +
   geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.9, alpha = 0.6) +
   geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.2) +
   geom_point(size = 3) +
   scale_colour_viridis_d(option = "viridis") +
   coord_flip() +
   labs(
     x = NULL,
     y = "Intercept (predicted larvae at reference year)",
     title = "Baseline larvae abundance by temperature (salinity held at mean)",
     subtitle = paste0(
       "(c) Intercepts evaluated at Year = ",
       round(lt_ref_year, 1),
       "; points = posterior mean; thick bars = 50% CrI; thin bars = 90% CrI"
     )
   ) +
   theme_bw(base_size = 18) +
   theme(
     panel.grid.minor = element_blank(),
     legend.position = "none"
   )
 
 lt_fig_year3temp_intercepts
 
 
 # Optional: show ribbons + intercepts + slopes together (requires patchwork loaded)
 # lt_alt_fig_year3temp_ribbons + lt_fig_year3temp_intercepts + lt_fig_year3temp_slopes
 
 
 # ============================================================
 # ltotal: FILTERED-DRAW "SPAGHETTI" VERSION IN RIBBON FORM
 # (same logic as your ltT filtered draws, but rendered as ribbons)
 # - no custom functions
 # - tidyverse + tidybayes
 # - lognormal model: epred is response scale (positive)
 # ============================================================
 
 # -----------------------------
 # ltF_00) Start from the draws you already created
 #   lt_alt_draws has: n_year, temp_level, .draw, lt_alt_epred
 # -----------------------------
 ltF_ep_long <- lt_alt_draws %>%
   transmute(
     .draw,
     n_year,
     temp_label = temp_level,
     epred = lt_alt_epred
   )
 
 # -----------------------------
 # ltF_10) Envelope (10–90%) at each (year × temp)
 # -----------------------------
 ltF_bounds80 <- ltF_ep_long %>%
   group_by(n_year, temp_label) %>%
   summarise(
     q10 = quantile(epred, 0.10, na.rm = TRUE),
     q90 = quantile(epred, 0.90, na.rm = TRUE),
     .groups = "drop"
   )
 
 # -----------------------------
 # ltF_11) Keep "well-behaved" draws (>=95% inside envelope)
 # -----------------------------
 ltF_draw_keep <- ltF_ep_long %>%
   left_join(ltF_bounds80, by = c("n_year", "temp_label")) %>%
   mutate(in80 = epred >= q10 & epred <= q90) %>%
   group_by(.draw) %>%
   summarise(prop_in80 = mean(in80, na.rm = TRUE), .groups = "drop") %>%
   filter(prop_in80 >= 0.95)
 
 # -----------------------------
 # ltF_12) Filter the epred draws to only "kept" draws
 # -----------------------------
 ltF_ep_keep <- ltF_ep_long %>%
   semi_join(ltF_draw_keep, by = ".draw")
 
 # -----------------------------
 # ltF_13) Summarise kept draws into ribbons + median line
 #   - outer ribbon = 10–90%
 #   - inner ribbon = 25–75%
 #   - median line  = 50%
 # -----------------------------
 ltF_rib_1090 <- ltF_ep_keep %>%
   group_by(n_year, temp_label) %>%
   summarise(
     med   = quantile(epred, 0.50, na.rm = TRUE),
     low10 = quantile(epred, 0.10, na.rm = TRUE),
     up90  = quantile(epred, 0.90, na.rm = TRUE),
     .groups = "drop"
   )
 
 ltF_rib_2575 <- ltF_ep_keep %>%
   group_by(n_year, temp_label) %>%
   summarise(
     low25 = quantile(epred, 0.25, na.rm = TRUE),
     up75  = quantile(epred, 0.75, na.rm = TRUE),
     .groups = "drop"
   )
 
 ltF_ribbon_dat <- ltF_rib_1090 %>%
   left_join(ltF_rib_2575, by = c("n_year", "temp_label"))
 
 # -----------------------------
 # ltF_14) Optional plotting y-limits (robust)
 #   - use kept draws so the plot focuses on the "well-behaved" range
 # -----------------------------
 ltF_ylim <- quantile(
   ltF_ep_keep$epred,
   probs = c(0.01, 0.99),
   na.rm = TRUE
 )
 
 # -----------------------------
 # ltF_15) Ribbon plot (filtered-draw version)
 #   - If you want log scale like your ltT spaghetti figure, keep trans="log10"
 # -----------------------------
 ltF_fig_time_temp_ribbon_filtered <- ggplot(ltF_ribbon_dat) +
   
   # outer ribbon (10–90%)
   geom_ribbon(
     aes(x = n_year, ymin = low10, ymax = up90, fill = temp_label),
     alpha = 0.18,
     colour = NA
   ) +
   
   # inner ribbon (25–75%)
   geom_ribbon(
     aes(x = n_year, ymin = low25, ymax = up75, fill = temp_label),
     alpha = 0.32,
     colour = NA
   ) +
   
   # median line
   geom_line(
     aes(x = n_year, y = med, colour = temp_label),
     linewidth = 1.1
   ) +
   
   scale_colour_viridis_d(option = "viridis") +
   scale_fill_viridis_d(option = "viridis") +
   scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
   
   # match your ltT style: log axis (since larvae totals are usually heavy-tailed)
   scale_y_continuous(
     trans  = "log10",
     breaks = c(4, 8, 16, 32, 64, 128, 256, 512, 1024),
     labels = scales::label_number()
   ) +
   
   coord_cartesian(ylim = ltF_ylim) +
   
   labs(
     x = "Monitoring year",
     y = "Total oyster larvae",
     title = "Temporal trends vary with temperature (salinity held at mean)",
     subtitle = "(a) Ribbons/line based on filtered posterior draws (>=95% inside the 10–90% envelope)"
   ) +
   
   theme_bw(base_size = 18) +
   theme(
     panel.grid.minor = element_blank(),
     legend.position = "bottom"
   )
 
 ltT_fig_time_temp_spag + ltF_fig_time_temp_ribbon_filtered
 
 
 