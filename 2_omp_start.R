


# testing out git with students!
library(tidyverse)
library(ggplot2)
library(viridis)
library(lubridate)
library(brms)
library(patchwork)
library(tidybayes)



# set your own personal working directory below
# Emma MBP
setwd("~/Dropbox/_Projects/PEI Oysters/Data/OMP")
# Maddy
setwd("~/Data/OMP/")


omp_dat <- read.csv("OMP_clean_2025.csv", header= TRUE)
head(omp_dat)

# ============================================================
# MODEL 2 (FIRST larvae event): 
# ============================================================

# ------------------------------------------------------------
# A) Build "first event" dataset (renamed)
# ------------------------------------------------------------

m2_first_l <- omp_dat %>%
  ungroup() %>%
  filter(larvae_total > 0) %>%                             # only events above zero
  group_by(bay, location_clean, f_year) %>%                # within year × location
  slice_min(order_by = julian_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(
    location_clean, area, bay, f_year, n_year, month_day,
    larvae_total, parsed_date, julian_date, water_temp, salinity
  ) %>%
  mutate(n_month_day = as.numeric(month_day))

m2_first_l_dat <- m2_first_l %>%
  mutate(
    bay = as.factor(bay),
    n_month_day = as.numeric(n_month_day)
  ) %>%
  group_by(bay) %>%
  mutate(
    water_temp.m = water_temp - mean(water_temp, na.rm = TRUE), #centering variables, to standarduze your data (no bias)
    salinity.m   = salinity   - mean(salinity,   na.rm = TRUE),
    n_year.m     = n_year     - mean(n_year,     na.rm = TRUE)
  )
# %>% filter(bay != "Tracadie Bay")

head(m2_first_l_dat)

# ------------------------------------------------------------
# B) Fit model (KEEP THIS NAME)
# ------------------------------------------------------------
# 
# start_temp_time_sal <- brm(
#   julian_date ~ water_temp.m * n_year.m * salinity.m +
#     (1 + water_temp.m * n_year.m || bay/location_clean),
#   data    = m2_first_l_dat,
#   iter    = 5000,
#   warmup  = 1000,
#   family  = gaussian(),
#   control = list(adapt_delta = 0.999, max_treedepth = 20)
# )

# Emma's paths
#save(start_temp_time_sal, file = "~/Data/Model_fits/OMP/start_temp_time_sal.Rdata")
load("~/Data/Model_fits/OMP/start_temp_time_sal.Rdata")
#load("~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/start_temp_time_sal.Rdata")

summary(start_temp_time_sal)
pp_check(start_temp_time_sal)
conditional_effects(start_temp_time_sal)
# ============================================================
# FIGURE WORKFLOW FOR start_temp_time_sal 
# ============================================================


# ============================================================
# 0) DEFINE SALINITY SLICES (CENTERED SCALE)
# ------------------------------------------------------------
# We define three salinity conditions for visualization:
#  - Low:   −1 SD
#  - Mean:   0
#  - High:  +1 SD
#
# These are representative slices for conditional visualization,
# NOT additional observed salinity values.
# ============================================================

m2_first_sal_sd <- sd(m2_first_l_dat$salinity.m, na.rm = TRUE)

m2_first_sal_df <- tibble(
  salinity.m = c(-m2_first_sal_sd, 0, m2_first_sal_sd),
  sal_label  = factor(
    c("Low salinity (−1 SD)", "Mean salinity", "High salinity (+1 SD)"),
    levels = c("Low salinity (−1 SD)", "Mean salinity", "High salinity (+1 SD)")
  )
)

# ============================================================
# 1) IDENTIFY YEAR LEVELS FROM THE FIT DATA
# ------------------------------------------------------------
# We pull year levels from m2_first_l_dat so colouring is stable
# and matches the years actually present in the fit dataset.
# ============================================================

m2_first_year_levels <- m2_first_l_dat %>%
  distinct(n_year) %>%
  arrange(n_year) %>%
  pull(n_year)

# ============================================================
# 2) BUILD PREDICTION GRID *ONLY FROM OBSERVED DATA*
# ------------------------------------------------------------
# Key restriction:
#   We use ONLY the observed (n_year, water_temp) combinations
#   that occurred in m2_first_l_dat (the fit dataset).
#
# Then we cross with the 3 salinity slices (panels) and re-create
# the centered predictors used in the model.
# ============================================================

m2_first_pred_grid <- m2_first_l_dat %>%
  distinct(n_year, water_temp) %>%      # <-- observed combinations only
  crossing(m2_first_sal_df) %>%
  mutate(
    # Center predictors exactly as used in the model
    water_temp.m = water_temp - mean(m2_first_l_dat$water_temp, na.rm = TRUE),
    n_year.m     = n_year     - mean(m2_first_l_dat$n_year,     na.rm = TRUE),
    
    # SAFE factor creation (prevents duplicated level errors)
    year_group   = factor(n_year, levels = sort(unique(n_year))),
    
    # Row ID needed to map predictions back to grid rows after pivoting
    row_id       = row_number()
  )

# ============================================================
# 3) POSTERIOR EXPECTED VALUES (POPULATION-LEVEL)
# ------------------------------------------------------------
# posterior_epred():
#   - returns expected julian_date on response scale
#   - includes residual variance
#   - re_formula = NA gives population-level predictions (no group RE)
#
# Output: draws × rows matrix
# ============================================================

m2_first_ep <- posterior_epred(
  start_temp_time_sal,
  newdata    = m2_first_pred_grid,
  re_formula = NA
)

# ============================================================
# 4) CONVERT POSTERIOR MATRIX → LONG FORMAT
# ------------------------------------------------------------
# This lets us compute posterior summaries using tidyverse verbs
# without any apply-family functions.
# ============================================================

m2_first_ep_long <- as_tibble(m2_first_ep) %>%
  mutate(.draw = row_number()) %>%       # draw index
  pivot_longer(
    cols      = - .draw,
    names_to  = "row_id",
    values_to = "epred"
  ) %>%
  mutate(row_id = readr::parse_number(row_id)) %>%
  left_join(
    m2_first_pred_grid %>%
      select(row_id, water_temp, n_year, year_group, sal_label),
    by = "row_id"
  )

# ============================================================
# 5) SUMMARISE POSTERIOR:
#    MEAN + 90% CREDIBLE INTERVAL
# ------------------------------------------------------------
# Summarise ACROSS DRAWS for each prediction point.
# ============================================================

m2_first_summ <- m2_first_ep_long %>%
  group_by(row_id, water_temp, n_year, year_group, sal_label) %>%
  summarise(
    estimate = mean(epred),
    lower90  = quantile(epred, 0.05),
    upper90  = quantile(epred, 0.95),
    .groups  = "drop"
  )

# ============================================================
# 6) FIGURE 1:
#    Temperature–phenology relationships (3 salinity panels)
# ------------------------------------------------------------
# Each coloured line = one year
# Ribbon = 90% CrI
# X-values = observed temperatures only
# ============================================================

m2_first_fig1_panel <- ggplot(
  m2_first_summ,
  aes(
    x = water_temp,
    y = estimate,
    group = year_group,
    colour = year_group
  )
) +
  geom_ribbon(
    aes(ymin = lower90, ymax = upper90),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ sal_label, nrow = 1) +
  scale_colour_viridis_d(option = "viridis", guide = "none") +
  labs(
    x = "Surface water temperature (°C)",
    y = "Julian date",
    title = "Temperature–phenology relationships vary through time and with salinity",
    subtitle = "FIRST EVENT: population-level predictions at observed temperatures; shaded bands show 90% credible intervals"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

m2_first_fig1_panel

# ============================================================
# 7) MEAN-SALINITY PANEL ONLY (FOR MAIN TEXT FIGURE)
# ------------------------------------------------------------
# Extract the mean-salinity slice for a cleaner single-panel plot.
# ============================================================

m2_first_summ_mid <- m2_first_summ %>%
  filter(sal_label == "Mean salinity")

m2_first_fig1_mid_start <- ggplot(
  m2_first_summ_mid,
  aes(
    x = water_temp,
    y = estimate,
    group = year_group,
    colour = year_group
  )
) +
  geom_ribbon(
    aes(ymin = lower90, ymax = upper90),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(linewidth = 0.7) +
  scale_colour_viridis_d(option = "viridis", guide = "none") +
  labs(
    x = "Surface water temperature (°C)",
    y = "Julian date (first event)",
    subtitle = "a) Date of first oyster larvae > 250 μm"
  ) +
  # NOTE:
  # If you want the same calendar labels as the max-event figure,
  # keep/adjust the limits/breaks below to match the observed range for this event.
  coord_cartesian() +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank())

m2_first_fig1_mid_start


# ------------------------------------------------------------
# 8) Spaghetti (mean salinity; central 80% region)
# ------------------------------------------------------------
# Same spaghetti logic as max-event figure:
#   - use mean-salinity only
#   - sample a pool of draws
#   - define row-wise 10–90% band
#   - keep draws mostly within that band
#   - plot ~150 “typical” draws
# ------------------------------------------------------------

set.seed(123)

m2_first_ep_mid_long <- m2_first_ep_long %>%
  filter(sal_label == "Mean salinity")

m2_first_ci_pool_size <- min(2000, n_distinct(m2_first_ep_mid_long$.draw))

m2_first_ci_pool_ids <- m2_first_ep_mid_long %>%
  distinct(.draw) %>%
  slice_sample(n = m2_first_ci_pool_size) %>%
  pull(.draw)

m2_first_ep_ci_long <- m2_first_ep_mid_long %>%
  filter(.draw %in% m2_first_ci_pool_ids)

m2_first_bounds80 <- m2_first_ep_ci_long %>%
  group_by(row_id) %>%
  summarise(
    q10 = quantile(epred, 0.10),
    q90 = quantile(epred, 0.90),
    .groups = "drop"
  )

m2_first_draw_keep <- m2_first_ep_ci_long %>%
  left_join(m2_first_bounds80, by = "row_id") %>%
  mutate(in80 = epred >= q10 & epred <= q90) %>%
  group_by(.draw) %>%
  summarise(prop_in80 = mean(in80), .groups = "drop") %>%
  filter(prop_in80 >= 0.95)

m2_first_n_spaghetti <- 150

m2_first_n_keep <- min(m2_first_n_spaghetti, nrow(m2_first_draw_keep))

m2_first_spaghetti_draw_ids <- m2_first_draw_keep %>%
  slice_sample(n = m2_first_n_keep) %>%
  pull(.draw)

m2_first_ep_spag_long <- m2_first_ep_ci_long %>%
  filter(.draw %in% m2_first_spaghetti_draw_ids)

m2_first_fig1_spag <- ggplot() +
  geom_line(
    data = m2_first_ep_spag_long,
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
    data = m2_first_summ_mid,
    aes(x = water_temp, y = estimate, group = year_group, colour = year_group),
    linewidth = 0.7
  ) +
  scale_colour_viridis_d(option = "viridis", guide = "none") +
  scale_y_continuous(
    limits = c(160, 230),
    breaks = c(166, 176, 186, 196, 206, 216, 226),
    labels = c("June 15","June 25","July 5","July 15",
               "July 25","Aug 4","Aug 14")
  ) +
  coord_cartesian() +
  labs(
    x = "Surface water temperature (°C)",
    y = "Julian date (first event)",
    title = "Temperature–phenology relationships at mean salinity",
    subtitle = "FIRST EVENT: thick lines = posterior means by year; spaghetti = posterior draws filtered to central 80% region"
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank())

m2_first_fig1_spag


# ------------------------------------------------------------
# 9) Intercepts-by-year (temperature-averaged)
#    3 salinity panels + mean-salinity-only
# ------------------------------------------------------------
# EXACTLY the same intercept logic:
#   - average across observed temperatures within each year
#   - do temperature averaging WITHIN each draw
#   - then summarise across draws for intervals
# ------------------------------------------------------------

m2_first_year_sal_temp_grid <- m2_first_l_dat %>%
  distinct(n_year, water_temp) %>%
  crossing(m2_first_sal_df) %>%
  mutate(
    n_year.m     = n_year     - mean(m2_first_l_dat$n_year,     na.rm = TRUE),
    water_temp.m = water_temp - mean(m2_first_l_dat$water_temp, na.rm = TRUE),
    row_id       = row_number()
  )

m2_first_ep2 <- posterior_epred(
  start_temp_time_sal,
  newdata    = m2_first_year_sal_temp_grid,
  re_formula = NA
)

m2_first_ep2_long <- as_tibble(m2_first_ep2) %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols      = - .draw,
    names_to  = "row_id",
    values_to = "epred"
  ) %>%
  mutate(row_id = readr::parse_number(row_id)) %>%
  left_join(
    m2_first_year_sal_temp_grid %>% select(row_id, n_year, sal_label, water_temp),
    by = "row_id"
  )

m2_first_draw_means <- m2_first_ep2_long %>%
  group_by(.draw, n_year, sal_label) %>%
  summarise(epred_mean = mean(epred), .groups = "drop")

m2_first_int_summ <- m2_first_draw_means %>%
  group_by(n_year, sal_label) %>%
  summarise(
    estimate = mean(epred_mean),
    lower50  = quantile(epred_mean, 0.25),
    upper50  = quantile(epred_mean, 0.75),
    lower90  = quantile(epred_mean, 0.05),
    upper90  = quantile(epred_mean, 0.95),
    .groups  = "drop"
  )

m2_first_fig_intercepts_3panel <- ggplot(
  m2_first_int_summ,
  aes(x = n_year, y = estimate)
) +
  geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.4) +
  facet_wrap(~ sal_label, nrow = 1) +
  scale_x_continuous(breaks = sort(unique(m2_first_l_dat$n_year))) +
  labs(
    x = "Monitoring year",
    y = "Julian date (temperature-averaged, population-level)",
    title = "Year-specific phenology across salinity regimes",
    subtitle = "FIRST EVENT: points = posterior means; thick bars = 50% CrI; thin bars = 90% CrI\nAveraged over observed temperatures within each year"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

m2_first_fig_intercepts_3panel

m2_first_fig_intercepts_mean <- m2_first_int_summ %>%
  filter(sal_label == "Mean salinity") %>%
  mutate(year_group = factor(n_year, levels = sort(unique(m2_first_l_dat$n_year)))) %>%
  ggplot(aes(x = n_year, y = estimate, colour = year_group)) +
  geom_linerange(aes(ymin = lower90, ymax = upper90, colour = year_group), linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50, colour = year_group), linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.4) +
  scale_colour_viridis_d(option = "viridis", guide = "none") +
  scale_x_continuous(breaks = sort(unique(m2_first_l_dat$n_year))) +
  labs(
    x = "Monitoring year",
    y = "Julian date (first event)",
    subtitle = "b) Date of first oyster larvae > 250 μm (Intercept; temp-averaged)"
  ) +
  scale_y_continuous(
    limits = c(160, 230),
    breaks = c(166, 176, 186, 196, 206, 216, 226),
    labels = c("June 15","June 25","July 5","July 15",
               "July 25","Aug 4","Aug 14")
  ) +
  coord_cartesian() +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank())


m2_first_fig1_spag + m2_first_fig_intercepts_mean



# ------------------------------------------------------------
# 10) Slopes-by-year (temperature effect; temp-averaged logic analog)
# ------------------------------------------------------------
# GOAL:
#   For each monitoring year, estimate the slope of julian_date vs water_temp
#   at each salinity slice (days per °C), using ONLY observed temps in that year.
#
# STRATEGY:
#   - Use m2_first_ep2_long (draws × observed year-temp rows × salinity)
#   - Within each (.draw, n_year, sal_label):
#       fit epred ~ water_temp  (using that year’s observed temps)
#       extract slope
#   - Then summarise across draws: mean + 50% and 90% CrI
#
# NOTE:
#   Slopes are on the response scale: "days change in phenology per 1°C".
# ------------------------------------------------------------

# wouldn't run 
# 10a) Draw-wise slopes per year and salinity
m2_first_draw_slopes <- m2_first_ep2_long %>%
  select(.draw, n_year, sal_label, water_temp, epred) %>%
  group_by(.draw, n_year, sal_label) %>%
  summarise(
    n_pts    = n_distinct(water_temp),
    wt_var   = var(water_temp, na.rm = TRUE),
    wt_mean  = mean(water_temp, na.rm = TRUE),
    yy_mean  = mean(epred,     na.rm = TRUE),
    wt_yy_cov = mean((water_temp - wt_mean) * (epred - yy_mean), na.rm = TRUE),
    slope = wt_yy_cov / wt_var,
    .groups = "drop"
  ) %>%
  filter(n_pts >= 2, wt_var > 0, is.finite(slope)) %>%
  select(.draw, n_year, sal_label, slope)


# 10b) Summarise slopes across draws: mean + 50% and 90% CrI
m2_first_slope_summ <- m2_first_draw_slopes %>%
  group_by(n_year, sal_label) %>%
  summarise(
    estimate = mean(slope),
    lower50  = quantile(slope, 0.25),
    upper50  = quantile(slope, 0.75),
    lower90  = quantile(slope, 0.05),
    upper90  = quantile(slope, 0.95),
    .groups  = "drop"
  )

# 10c) Slope figure: 3 salinity panels (like intercepts)
m2_first_fig_slopes_3panel <- ggplot(
  m2_first_slope_summ,
  aes(x = n_year, y = estimate)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, alpha = 0.5) +
  geom_linerange(aes(ymin = lower90, ymax = upper90),
                 linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50),
                 linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.4) +
  facet_wrap(~ sal_label, nrow = 1) +
  scale_x_continuous(breaks = sort(unique(m2_first_l_dat$n_year))) +
  labs(
    x = "Monitoring year",
    y = "Temperature slope (days per °C)",
    title = "Year-specific temperature sensitivity across salinity regimes",
    subtitle = "FIRST EVENT: points = posterior means; thick bars = 50% CrI; thin bars = 90% CrI\nSlope estimated within each draw using observed temperatures within each year"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

m2_first_fig_slopes_3panel


# 10d) Mean-salinity-only slope figure (styled like your mean intercept panel)
m2_first_fig_slopes_mean <- m2_first_slope_summ %>%
  filter(sal_label == "Mean salinity") %>%
  mutate(year_group = factor(n_year, levels = sort(unique(m2_first_l_dat$n_year)))) %>%
  ggplot(aes(x = n_year, y = estimate, colour = year_group)) +
  geom_hline(yintercept = 0, linewidth = 0.4, alpha = 0.5) +
  geom_linerange(aes(ymin = lower90, ymax = upper90, colour = year_group),
                 linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50, colour = year_group),
                 linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.4) +
  scale_colour_viridis_d(option = "viridis", guide = "none") +
  scale_x_continuous(breaks = sort(unique(m2_first_l_dat$n_year))) +
  labs(
    x = "Monitoring year",
    y = "Temperature slope (days per °C)",
    subtitle = "c) Temperature sensitivity of first larvae date (Slope; temp effect within-year)"
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank())

# Example combo:
m2_first_fig1_spag + m2_first_fig_intercepts_mean + m2_first_fig_slopes_mean


