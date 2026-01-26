# how to clean your global environment
rm(list = ls())


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
    water_temp.m = water_temp - mean(water_temp, na.rm = TRUE), #centering variables, to standardize your data (no bias)
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
#     (1 + n_year.m | bay/location_clean),
#   data    = m2_first_l_dat,
#   iter    = 5000,
#   warmup  = 1000,
#   family  = gaussian(),
#   control = list(adapt_delta = 0.999, max_treedepth = 20)
# )

# Emma's paths
save(start_temp_time_sal, file = "~/Dropbox/_Projects/PEI Oysters/Model_Fits/OMP/start_temp_time_sal.Rdata")
load("~/Dropbox/_Projects/PEI Oysters/Model_Fits/OMP/start_temp_time_sal.Rdata")
#load("~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/start_temp_time_sal.Rdata")

# Maddy's path
#save(start_temp_time_sal, file = "~/Data/Model_fits/OMP/start_temp_time_sal.Rdata")
load("~/Data/Model_fits/OMP/start_temp_time_sal.Rdata")


summary(start_temp_time_sal)
pp_check(start_temp_time_sal)

color_scheme_set("darkgray")
fig_pp_start <- pp_check(start_temp_time_sal) +   
  xlab( "Julian date") + ylab("Density") + 
  ggtitle('First date of oyster larvae (Julian date)')+
  labs(subtitle= "a)")+
  theme_classic() +  theme( plot.title=element_text(size=18, hjust=0.5), legend.position= "none")# predicted vs. observed values
fig_pp_start

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
  scale_colour_viridis_d(option = "viridis", name = "Monitoring year") +
  labs(
    x = "Surface water temperature (°C)",
    y = "Julian date",
    # title = "Temperature–phenology relationships vary through time and with salinity",
    # subtitle = "FIRST EVENT: population-level predictions at observed temperatures; shaded bands show 90% credible intervals"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom")

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
  scale_colour_viridis_d(option = "viridis", name = "Monitoring year") +
  facet_wrap(~ sal_label, nrow = 1) +
  labs(
    x = "Surface water temperature (°C)",
    y = "Julian date",
    subtitle = "a)" # Date of first detection of oyster larvae (> 250 μm)
  ) +
  # NOTE:
  # If you want the same calendar labels as the max-event figure,
  # keep/adjust the limits/breaks below to match the observed range for this event.
  coord_cartesian() +
  theme_bw(base_size = 18) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

m2_first_fig1_mid_start

# what does this figure tell us?
# under average salinity coniditions warmer surface temps are associated with an earlier appearance of oyster larvae above 250um
# the strenght (slope) and timing of this relationships vary from year to year


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
    aes(x = water_temp, y = estimate, group = factor(n_year), colour = factor(n_year)),
    linewidth = 0.7
  ) +
  scale_colour_viridis_d(option = "viridis", name = "Monitoring year") +
  scale_y_continuous(
    limits = c(160, 230),
    breaks = c(166, 176, 186, 196, 206, 216, 226),
    labels = c("June 15","June 25","July 5","July 15",
               "July 25","Aug 4","Aug 14")
  ) +
  coord_cartesian() +
  labs(
    x = "Surface water temperature (°C)",
    y = "Julian date",
    # title = "Temperature effects on first oyster larvae appearance above 250um at mean salinity",
    subtitle = "a)" # "FIRST EVENT: thick lines = posterior means by year; spaghetti = posterior draws filtered to central 80% region"
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

m2_first_fig1_spag

# what is this figure telling us?
# same at part 7)
# under average salinity coniditions warmer surface temps are associated with an earlier appearance of oyster larvae above 250um

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
    estimate = round(mean(epred_mean), 2),
    lower50  = round(quantile(epred_mean, 0.25), 2),
    upper50  = round(quantile(epred_mean, 0.75), 2),
    lower90  = round(quantile(epred_mean, 0.05), 2),
    upper90  = round(quantile(epred_mean, 0.95), 2),
    .groups  = "drop"
  )

view(m2_first_int_summ)
write.csv(m2_first_int_summ, "~/Data/Results/m2_first_int_summ.csv", row.names = FALSE)

m2_first_fig_intercepts_3panel <- ggplot(
  m2_first_int_summ,
  aes(x = n_year, y = estimate, colour = factor(n_year), group = factor(n_year))
) +
  geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.4) +
  facet_wrap(~ sal_label, nrow = 1) +
  scale_colour_viridis_d(option = "viridis", name = "Monitoring year") +
  scale_x_continuous(breaks = sort(unique(m2_first_l_dat$n_year))[seq(1, length(unique(m2_first_l_dat$n_year)), by = 2)]) +
  labs(
    x = "Monitoring year",
    y = "Julian date",
    # title = "Year-specific phenology across salinity regimes",
    # subtitle = "FIRST EVENT: points = posterior means; thick bars = 50% CrI; thin bars = 90% CrI\nAveraged over observed temperatures within each year"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  theme(legend.position = "bottom")

m2_first_fig_intercepts_3panel

m2_first_fig_intercepts_mean <- m2_first_int_summ %>%
  filter(sal_label == "Mean salinity") %>%
  mutate(year_group = factor(n_year, levels = sort(unique(m2_first_l_dat$n_year)))) %>%
  ggplot(aes(x = n_year, y = estimate, colour = year_group)) +
  geom_linerange(aes(ymin = lower90, ymax = upper90, colour = year_group), linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50, colour = year_group), linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.4) +
  scale_colour_viridis_d(option = "viridis", guide = "none") +
  scale_x_continuous(breaks = sort(unique(m2_first_l_dat$n_year))[seq(1, length(unique(m2_first_l_dat$n_year)), by = 2)]) +
  labs(
    x = "Monitoring year",
    y = "Julian date",
    subtitle = "b)" # Date of first oyster larvae > 250 μm (Intercept; temp-averaged)
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

# yes it runs now
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
  aes(x = n_year, y = estimate, colour = factor(n_year), group = factor(n_year))
) +
  geom_hline(yintercept = 0, linewidth = 0.4, alpha = 0.5) +
  geom_linerange(aes(ymin = lower90, ymax = upper90),
                 linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50),
                 linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.4) +
  facet_wrap(~ sal_label, nrow = 1) +
  scale_colour_viridis_d(option = "viridis", name = "Monitoring year") +
  scale_x_continuous(breaks = sort(unique(m2_first_l_dat$n_year))[
    sort(unique(m2_first_l_dat$n_year)) %% 2 == 1 &
      sort(unique(m2_first_l_dat$n_year)) >= 2013]) +
  labs(
    x = "Monitoring year",
    y = "Temperature slope" # (days per °C)",
    # title = "Year-specific temperature sensitivity across salinity regimes",
    # subtitle = "FIRST EVENT: points = posterior means; thick bars = 50% CrI; thin bars = 90% CrI\nSlope estimated within each draw using observed temperatures within each year"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
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
  scale_colour_viridis_d(option = "viridis", guide = "Monitoring year") +
  scale_x_continuous(breaks = sort(unique(m2_first_l_dat$n_year))[seq(1, length(unique(m2_first_l_dat$n_year)), by = 2)]) +
  labs(
    x = "Monitoring year",
    y = "Temperature slope" # (days per °C)",
    subtitle = "c)" # Temperature sensitivity of first larvae date (Slope; temp effect within-year)"
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank())

# Example combo:

# Kill the legends for graph b) and c)
m2_first_fig1_spag <- m2_first_fig1_spag + theme(legend.position = "none")
m2_first_fig_intercepts_mean <- m2_first_fig_intercepts_mean + theme(legend.position = "none")


(m2_first_fig1_spag +
    m2_first_fig_intercepts_mean +
    m2_first_fig_slopes_mean) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  )




# ============================================================
# ALT spaghetti fig: year × temperature (salinity held at mean)
# For: start_temp_time_sal (data = m2_first_l_dat)
# - no custom functions
# - tidyverse + tidybayes
# - namespaced objects to avoid collisions
# ============================================================

set.seed(123)

# --- references for centering + temp quantiles (UNCENTERED vars) ---
m2_start_alt_ref <- m2_first_l_dat %>%
  ungroup() %>%
  summarise(
    n_year_mean     = mean(n_year, na.rm = TRUE),
    water_temp_mean = mean(water_temp, na.rm = TRUE),
    sal_mean        = mean(salinity, na.rm = TRUE),
    temp_cool       = quantile(water_temp, 0.10, na.rm = TRUE),
    temp_med        = quantile(water_temp, 0.50, na.rm = TRUE),
    temp_warm       = quantile(water_temp, 0.90, na.rm = TRUE)
  )

# --- temperature level lookup table ---
m2_start_alt_temp_levels <- tibble(
  temp_level = factor(c("Cool temp (25th pct)", "Median temp", "Warm temp (75th pct)"),
                      levels = c("Cool temp (25th pct)", "Median temp", "Warm temp (75th pct)")),
  water_temp = c(m2_start_alt_ref$temp_cool,
                 m2_start_alt_ref$temp_med,
                 m2_start_alt_ref$temp_warm)
)

# --- prediction grid (year sequence × temp levels), salinity fixed at mean ---
m2_start_alt_newdat_year_temp <- crossing(
  n_year = seq(min(m2_first_l_dat$n_year, na.rm = TRUE),
               max(m2_first_l_dat$n_year, na.rm = TRUE),
               length.out = 200),
  temp_level = m2_start_alt_temp_levels$temp_level
) %>%
  left_join(m2_start_alt_temp_levels, by = "temp_level") %>%
  mutate(
    salinity     = m2_start_alt_ref$sal_mean,
    salinity.m   = salinity   - m2_start_alt_ref$sal_mean,
    n_year.m     = n_year     - m2_start_alt_ref$n_year_mean,
    water_temp.m = water_temp - m2_start_alt_ref$water_temp_mean
  )

# --- posterior expected draws (population-level) ---
m2_start_alt_draws <- start_temp_time_sal %>%
  add_epred_draws(newdata = m2_start_alt_newdat_year_temp, re_formula = NA) %>%
  ungroup() %>%
  rename(m2_start_alt_epred = .epred)

# --- summarise across draws (median + 80% interval) ---
m2_start_alt_summ <- m2_start_alt_draws %>%
  group_by(n_year, temp_level) %>%
  median_qi(m2_start_alt_epred, .width = 0.80) %>%
  ungroup()

# --- thin spaghetti: sample N draws total (not per temp; consistent across temps) ---
m2_start_alt_n_keep <- 60

m2_start_alt_keep_draws <- m2_start_alt_draws %>%
  distinct(.draw) %>%
  slice_sample(n = min(m2_start_alt_n_keep, nrow(.)))

m2_start_alt_draws_thin <- m2_start_alt_draws %>%
  semi_join(m2_start_alt_keep_draws, by = ".draw")

# --- optional: set y-limits to data-informed range (instead of hard-coding) ---
m2_start_alt_ylim <- m2_first_l_dat %>%
  summarise(
    ymin = quantile(julian_date, 0.01, na.rm = TRUE),
    ymax = quantile(julian_date, 0.99, na.rm = TRUE)
  )


# --- plot ---
m2_start_alt_fig_year3temp_spag_thin <- ggplot() +
  geom_line(
    data = m2_start_alt_draws_thin,
    aes(
      x = n_year, y = m2_start_alt_epred,
      group = interaction(.draw, temp_level),
      colour = temp_level
    ),
    linewidth = 0.35,
    alpha = 0.10
  ) +
  geom_line(
    data = m2_start_alt_summ,
    aes(
      x = n_year, y = m2_start_alt_epred,
      group = temp_level,
      colour = temp_level
    ),
    linewidth = 1.0
  ) +
  scale_colour_manual(
    name = "Temp",
    values = c(
      "Cool temp (25th pct)" = "#1B9E77",  # green (Dark2)
      "Median temp"          = "#7570B3",  # purple (Dark2)
      "Warm temp (75th pct)" = "#D95F02"   # orange (Dark2)
    )
  ) +
  #scale_colour_viridis_d(option = "viridis", name = "Temp") +
  scale_x_continuous(
    breaks = seq(
      from = 2013,
      to   = max(m2_start_alt_summ$n_year, na.rm = TRUE),
      by   = 2
    ),
    labels = scales::label_number(accuracy = 1)
  ) +
  coord_cartesian(
    ylim = quantile(
      m2_first_l_dat$julian_date,
      probs = c(0.01, 0.99),
      na.rm = TRUE
    )
  )+
  # if you want the same "date labels" trick, keep this block and adjust breaks/labels if needed:
 # scale_y_continuous(
  #  breaks = c(184, 188, 192, 196, 200, 204, 208, 212, 215, 218, 222),
   # labels = c("Jul 3","Jul 7","Jul 11","Jul 15","Jul 19","Jul 23",
    #           "Jul 27","Jul 31","Aug 2","Aug 5","Aug 9")
  #) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 6),
    labels = scales::label_number(accuracy = 1)
  ) +
  labs(
    x = "Year",
    y = "Julian date",
    # title = "Time–phenology relationships at cool/median/warm temperatures",
    subtitle = "a)" #paste0(
      #"Thick lines = posterior medians; spaghetti = ",
      #m2_start_alt_n_keep,
      #" posterior draws (salinity held at mean)"
    #)
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

m2_start_alt_fig_year3temp_spag_thin





# ============================================================
# ALT ribbon fig: year × temperature (salinity held at mean)
# For: start_temp_time_sal (data = m2_first_l_dat)
# - no custom functions
# - tidyverse + tidybayes
# - ribbons = posterior uncertainty (80% + 50%)
# ============================================================

set.seed(123)

# --- references for centering + temp quantiles (UNCENTERED vars) ---
m2_start_alt_ref <- m2_first_l_dat %>%
  ungroup() %>%
  summarise(
    n_year_mean     = mean(n_year, na.rm = TRUE),
    water_temp_mean = mean(water_temp, na.rm = TRUE),
    sal_mean        = mean(salinity, na.rm = TRUE),
    temp_cool       = quantile(water_temp, 0.10, na.rm = TRUE),
    temp_med        = quantile(water_temp, 0.50, na.rm = TRUE),
    temp_warm       = quantile(water_temp, 0.90, na.rm = TRUE)
  )

# --- temperature level lookup table ---
m2_start_alt_temp_levels <- tibble(
  temp_level = factor(c("Cool", "Median", "Warm"),
                      levels = c("Cool", "Median", "Warm")),
  water_temp = c(m2_start_alt_ref$temp_cool,
                 m2_start_alt_ref$temp_med,
                 m2_start_alt_ref$temp_warm)
)

# --- prediction grid (year sequence × temp levels), salinity fixed at mean ---
m2_start_alt_newdat_year_temp <- crossing(
  n_year = seq(min(m2_first_l_dat$n_year, na.rm = TRUE),
               max(m2_first_l_dat$n_year, na.rm = TRUE),
               length.out = 200),
  temp_level = m2_start_alt_temp_levels$temp_level
) %>%
  left_join(m2_start_alt_temp_levels, by = "temp_level") %>%
  mutate(
    salinity     = m2_start_alt_ref$sal_mean,
    salinity.m   = salinity   - m2_start_alt_ref$sal_mean,
    n_year.m     = n_year     - m2_start_alt_ref$n_year_mean,
    water_temp.m = water_temp - m2_start_alt_ref$water_temp_mean
  )

# --- posterior expected draws (population-level) ---
m2_start_alt_draws <- start_temp_time_sal %>%
  add_epred_draws(newdata = m2_start_alt_newdat_year_temp, re_formula = NA) %>%
  ungroup() %>%
  rename(m2_start_alt_epred = .epred)

# --- summarise across draws:
#     80% interval for wide ribbon + 50% interval for inner ribbon + median line
m2_start_alt_summ80 <- m2_start_alt_draws %>%
  group_by(n_year, temp_level) %>%
  median_qi(m2_start_alt_epred, .width = 0.80) %>%
  ungroup() %>%
  rename(
    m2_start_alt_med   = m2_start_alt_epred,
    m2_start_alt_low80 = .lower,
    m2_start_alt_up80  = .upper
  )

m2_start_alt_summ50 <- m2_start_alt_draws %>%
  group_by(n_year, temp_level) %>%
  median_qi(m2_start_alt_epred, .width = 0.50) %>%
  ungroup() %>%
  rename(
    m2_start_alt_med50 = m2_start_alt_epred,
    m2_start_alt_low50 = .lower,
    m2_start_alt_up50  = .upper
  )

# --- join 80 + 50 for plotting convenience ---
m2_start_alt_ribbon_dat <- m2_start_alt_summ80 %>%
  left_join(
    m2_start_alt_summ50 %>%
      select(n_year, temp_level, m2_start_alt_low50, m2_start_alt_up50),
    by = c("n_year", "temp_level")
  )

# --- plot (ribbons + median line) ---
m2_start_alt_fig_year3temp_ribbons <- ggplot(m2_start_alt_ribbon_dat) +
  
  # outer ribbon (80%)
  geom_ribbon(
    aes(
      x = n_year,
      ymin = m2_start_alt_low80,
      ymax = m2_start_alt_up80,
      fill = temp_level
    ),
    alpha = 0.18,
    colour = NA
  ) +
  
  # inner ribbon (50%)
  geom_ribbon(
    aes(
      x = n_year,
      ymin = m2_start_alt_low50,
      ymax = m2_start_alt_up50,
      fill = temp_level
    ),
    alpha = 0.32,
    colour = NA
  ) +
  
  # median line
  geom_line(
    aes(
      x = n_year,
      y = m2_start_alt_med,
      colour = temp_level
    ),
    linewidth = 1.0
  ) +
  scale_x_continuous(
    breaks = seq(
      from = 2013,
      to   = max(m2_start_alt_ribbon_dat$n_year, na.rm = TRUE),
      by   = 2
    ),
    labels = scales::label_number(accuracy = 1)
  ) +
  scale_colour_manual(
    name = "Temp",
    values = c(
      "Cool"   = "#1B9E77",  # green (Dark2)
      "Median" = "#7570B3",  # purple (Dark2)
      "Warm"   = "#D95F02"   # orange (Dark2)
    )
  ) +
  scale_fill_manual(
    name = "Temp",
    values = c(
      "Cool"   = "#1B9E77",
      "Median" = "#7570B3",
      "Warm"   = "#D95F02"
    )
  ) +

  #scale_colour_viridis_d(option = "viridis") +
  #scale_fill_viridis_d(option = "viridis") +
  
  coord_cartesian(
    ylim = quantile(
      m2_first_l_dat$julian_date,
      probs = c(0.01, 0.99),
      na.rm = TRUE
    )
  ) +
  
  # keep / edit this if your julian-date-to-calendar mapping differs
  #scale_y_continuous(
   # breaks = c(184, 188, 192, 196, 200, 204, 208, 212, 215, 218, 222),
    #labels = c("Jul 3","Jul 7","Jul 11","Jul 15","Jul 19","Jul 23",
     #          "Jul 27","Jul 31","Aug 2","Aug 5","Aug 9")
  #) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 6),
    labels = scales::label_number(accuracy = 1)
  ) +
  labs(
    x = "Year",
    y = "Julian date",
    # title = "Time–phenology relationships at cool/median/warm temperatures",
    subtitle = "a)", #"Lines = posterior medians; ribbons = 50% (inner) and 80% (outer) credible intervals (salinity held at mean)"
  colour = "Temp",
  fill = "Temp"
    ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

m2_start_alt_fig_year3temp_ribbons



# ============================================================
# m2 START: SLOPE FIGURE (matched to year × temp ribbons figure)
# - uses draw-by-draw slopes within temp_level
# - summarises across draws (mean + 50%/90% CrI)
# ============================================================

# -----------------------------
# m2_start_14) Slopes within each draw × temp level
# slope = cov(x,y)/var(x)  (simple linear slope, draw-by-draw)
# -----------------------------
m2_start_slope_draws <- m2_start_alt_draws %>%
  group_by(.draw, temp_level) %>%
  summarise(
    slope = cov(n_year, m2_start_alt_epred) / var(n_year),
    .groups = "drop"
  )

# -----------------------------
# m2_start_15) Summarise slopes across draws (mean + 50%/90% CrI)
# -----------------------------
m2_start_slope_summ <- m2_start_slope_draws %>%
  group_by(temp_level) %>%
  summarise(
    estimate = mean(slope),
    lower50  = quantile(slope, 0.25),
    upper50  = quantile(slope, 0.75),
    lower90  = quantile(slope, 0.05),
    upper90  = quantile(slope, 0.95),
    .groups  = "drop"
  )

# -----------------------------
# m2_start_16) Slope figure
# -----------------------------
m2_start_fig_year3temp_slopes <- ggplot(
  m2_start_slope_summ,
  aes(y = temp_level, x = estimate, colour = temp_level)
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(xmin = lower90, xmax = upper90), linewidth = 0.9, alpha = 0.6) +
  geom_linerange(aes(xmin = lower50, xmax = upper50), linewidth = 2.2) +
  geom_point(size = 3) +
  scale_colour_manual(
    name = "Temp",
    values = c(
      "Cool"   = "#1B9E77",  # green (Dark2)
      "Median" = "#7570B3",  # purple (Dark2)
      "Warm"   = "#D95F02"   # orange (Dark2)
    )
  ) +

  #scale_colour_viridis_d(option = "viridis") +
  coord_flip() +
  labs(
    y = NULL,
    x = "Slope", # (change in predicted Julian date per year)",
    #title = "Temporal trends by temperature (salinity held at mean)",
    subtitle = "c)" # Points = posterior mean slopes; thick bars = 50% CrI; thin bars = 90% CrI"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

m2_start_fig_year3temp_slopes

# Optional: show ribbons + slopes together (requires patchwork loaded)
# m2_start_alt_fig_year3temp_ribbons + m2_start_fig_year3temp_slopes


# ============================================================
# m2 START: INTERCEPT FIGURE (matched to ribbons + slopes style)
# - draw-by-draw intercepts within temp_level
# - intercept = predicted value at a chosen reference year
# - summarises across draws (mean + 50%/90% CrI)
# ============================================================

# -----------------------------
# m2_start_17) Choose the reference year for the intercept
#   - using the mean year keeps it "central" and comparable
# -----------------------------
m2_start_ref_year <- m2_start_alt_ref$n_year_mean

# -----------------------------
# m2_start_18) Intercepts within each draw × temp level
#   - estimate epred at the reference year
#   - we use the existing draw grid; just filter to closest year
# -----------------------------
m2_start_intercept_draws <- m2_start_alt_draws %>%
  mutate(.dist = abs(n_year - m2_start_ref_year)) %>%
  group_by(.draw, temp_level) %>%
  slice_min(.dist, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    .draw,
    temp_level,
    intercept = m2_start_alt_epred
  )

# -----------------------------
# m2_start_19) Summarise intercepts across draws (mean + 50%/90% CrI)
# -----------------------------
m2_start_intercept_summ <- m2_start_intercept_draws %>%
  group_by(temp_level) %>%
  summarise(
    estimate = mean(intercept),
    lower50  = quantile(intercept, 0.25),
    upper50  = quantile(intercept, 0.75),
    lower90  = quantile(intercept, 0.05),
    upper90  = quantile(intercept, 0.95),
    .groups  = "drop"
  )

# -----------------------------
# m2_start_20) Intercept figure
# -----------------------------
m2_start_fig_year3temp_intercepts <- ggplot(
  m2_start_intercept_summ,
  aes(y = temp_level, x = estimate, colour = temp_level)
) +
  geom_linerange(aes(xmin = lower90, xmax = upper90), linewidth = 0.9, alpha = 0.6) +
  geom_linerange(aes(xmin = lower50, xmax = upper50), linewidth = 2.2) +
  geom_point(size = 3) +
  #scale_colour_viridis_d(option = "viridis") +
  scale_colour_manual(
    values = c(
      "Cool"   = "#1B9E77",  # green (Dark2)
      "Median" = "#7570B3",  # purple (Dark2)
      "Warm"   = "#D95F02"   # orange (Dark2)
    )
  ) +
  coord_flip() +
  labs(
    x = "Intercept (predicted Julian date at reference year)",
    y = NULL,
    #title = "Baseline phenology by temperature (salinity held at mean)",
    subtitle = "b)"
    ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

m2_start_fig_year3temp_intercepts

# Optional: show ribbons + slopes + intercepts together (requires patchwork loaded)
# m2_start_alt_fig_year3temp_ribbons + m2_start_fig_year3temp_slopes + m2_start_fig_year3temp_intercepts

#remove legend
m2_start_fig_year3temp_intercepts <- m2_start_fig_year3temp_intercepts +
guides(colour = "none", fill = "none") +
  theme(legend.position = "none")

#remove legend
m2_start_fig_year3temp_slopes <- m2_start_fig_year3temp_slopes +
guides(colour = "none", fill = "none") +
  theme(legend.position = "none")

(m2_start_alt_fig_year3temp_ribbons +
  m2_start_fig_year3temp_intercepts +
  m2_start_fig_year3temp_slopes) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center")

m2_start_alt_fig_year3temp_spag_thin+ m2_start_alt_fig_year3temp_ribbons
