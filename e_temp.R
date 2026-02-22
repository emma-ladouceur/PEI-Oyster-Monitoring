
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
nrow(omp_dat)


e_dat <- omp_dat %>%
  select(bay, location_clean, julian_date, water_temp, salinity, n_year) %>%
  filter(!is.na(water_temp), !is.na(salinity), !is.na(n_year), !is.na(julian_date)) %>%
  group_by(n_year) %>%
  mutate(
    # within-year season index (0 = first sampled day in that year)
    j_season = julian_date - min(julian_date, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    # center + scale for Stan friendliness
    n_year.z    = as.numeric(scale(n_year)),        # mean 0, sd 1
    j_season.z  = as.numeric(scale(j_season)),      # mean 0, sd 1 (or use 0–1 scaling below)
    water_temp.z = as.numeric(scale(water_temp)),
    salinity.z   = as.numeric(scale(salinity))
  )

head(e_dat)

pri <- c(
  prior(exponential(1), class = sd),      # group-level SDs
  prior(exponential(1), class = sds),     # smooth SD
  prior(student_t(3, 0, 5), class = sigma),
  prior(normal(0, 1), class = b)          # fixed effects (n_year.m)
)


e_t_mod <- brm(
  water_temp ~
    n_year.z +
    s(j_season.z, k = 8) +
    (1 | bay/location_clean) +
    (0 + n_year.z | bay/location_clean),
  data = e_dat,
  family = gaussian(),
  prior = pri,
  iter = 5000,
  warmup = 1000,
  control = list(adapt_delta = 0.999, max_treedepth = 20)
)


conditional_effects(e_t_mod)
pp_check(e_t_mod)


# Emma's path
save(e_t_mod, file = "~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/e_t_mod.Rdata")
load("~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/e_t_mod.Rdata")

# Maddy's path
save(e_t_mod, file = "~/Data/Model_fits/OMP/e_t_mod.Rdata")
load("~/Data/Model_fits/OMP/e_t_mod.Rdata")



# ============================================================
# PUBLICATION-READY FIGURES for TEMPERATURE MODEL (e_t_mod)
# - Panel a) Seasonal pattern (population band + bay/location lines)
# - Panel b) Long-term trend through years at mid-season
# - Panel c) Bay-level slopes (°C per year) at mid-season
#   (computed correctly: location slopes within draw -> bay avg within draw)
# ============================================================


# ------------------------------------------------------------
# 0) Exact data used to fit the model (avoids new-level issues)
# ------------------------------------------------------------
t_fit_df <- model.frame(e_t_mod)

# ------------------------------------------------------------
# 1) Fixed bay order + colour palette
# ------------------------------------------------------------
t_bay_levels <- t_fit_df %>%
  distinct(bay) %>%
  arrange(bay) %>%
  pull(bay)

t_bay_colors <- setNames(
  viridisLite::viridis(length(t_bay_levels), option = "plasma"),
  t_bay_levels
)

# ------------------------------------------------------------
# 2) Helpers to convert z-scales back to natural units
#    (must match how you created n_year.z and j_season.z)
# ------------------------------------------------------------
year_mu <- mean(e_dat$n_year, na.rm = TRUE)
year_sd <- sd(e_dat$n_year, na.rm = TRUE)

j_mu    <- mean(e_dat$j_season, na.rm = TRUE)
j_sd    <- sd(e_dat$j_season, na.rm = TRUE)

z_to_year <- function(z) z * year_sd + year_mu
z_to_j    <- function(z) z * j_sd + j_mu

# ============================================================
# PANEL a) Seasonal curves (year held at mean: n_year.z = 0)
# ============================================================

# grid: real bay/location groups only
t_season_grid_loc <- t_fit_df %>%
  distinct(bay, location_clean) %>%
  mutate(bay = factor(bay, levels = t_bay_levels)) %>%
  crossing(
    j_season.z = seq(min(t_fit_df$j_season.z), max(t_fit_df$j_season.z), length.out = 160),
    n_year.z   = 0
  ) %>%
  mutate(j_season = z_to_j(j_season.z))

# fitted incl. bay/location effects
t_season_fit_loc <- fitted(
  e_t_mod,
  newdata    = t_season_grid_loc,
  re_formula = NULL,
  probs      = c(0.10, 0.90)
) %>%
  as_tibble() %>%
  bind_cols(t_season_grid_loc)

# population-only seasonal curve
t_season_grid_pop <- tibble(
  j_season.z = seq(min(t_fit_df$j_season.z), max(t_fit_df$j_season.z), length.out = 200),
  n_year.z   = 0
) %>%
  mutate(j_season = z_to_j(j_season.z))

t_season_fit_pop <- fitted(
  e_t_mod,
  newdata    = t_season_grid_pop,
  re_formula = NA,
  probs      = c(0.10, 0.90)
) %>%
  as_tibble() %>%
  bind_cols(t_season_grid_pop)

t_fig_season <- ggplot() +
  geom_ribbon(data = t_season_fit_pop, aes(x = j_season, ymin = Q10, ymax = Q90), alpha = 0.30) +
  geom_line(  data = t_season_fit_pop, aes(x = j_season, y = Estimate), colour = "black", linewidth = 1.2) +
  geom_line(
    data = t_season_fit_loc,
    aes(x = j_season, y = Estimate, group = interaction(bay, location_clean), colour = bay),
    linewidth = 0.6, alpha = 0.80
  ) +
  scale_colour_manual(values = t_bay_colors) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  labs(
    x = "Day of season (Julian day relative to first sampled day in year)",
    y = "Surface water temperature (°C)",
    colour = "Bay",
    subtitle = "a) Seasonal pattern (year held at mean)"
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom")

t_fig_season


# ============================================================
# PANEL b) Long-term trend through years at MID-SEASON
# ============================================================

j_ref <- median(t_fit_df$j_season.z, na.rm = TRUE)

# grid over years for each bay/location at mid-season
t_year_grid_loc <- t_fit_df %>%
  distinct(bay, location_clean) %>%
  crossing(
    n_year.z   = seq(min(t_fit_df$n_year.z), max(t_fit_df$n_year.z), length.out = 160),
    j_season.z = j_ref
  ) %>%
  mutate(
    bay   = factor(bay, levels = t_bay_levels),
    n_year = z_to_year(n_year.z)
  )

t_year_fit_loc <- fitted(
  e_t_mod,
  newdata    = t_year_grid_loc,
  re_formula = NULL,
  probs      = c(0.10, 0.90)
) %>%
  as_tibble() %>%
  bind_cols(t_year_grid_loc)

# population trend at mid-season
t_year_grid_pop <- tibble(
  n_year.z   = seq(min(t_fit_df$n_year.z), max(t_fit_df$n_year.z), length.out = 200),
  j_season.z = j_ref
) %>%
  mutate(n_year = z_to_year(n_year.z))

t_year_fit_pop <- fitted(
  e_t_mod,
  newdata    = t_year_grid_pop,
  re_formula = NA,
  probs      = c(0.10, 0.90)
) %>%
  as_tibble() %>%
  bind_cols(t_year_grid_pop)

t_fig_trend <- ggplot() +
  geom_ribbon(data = t_year_fit_pop, aes(x = n_year, ymin = Q10, ymax = Q90), alpha = 0.30) +
  geom_line(  data = t_year_fit_pop, aes(x = n_year, y = Estimate), colour = "black", linewidth = 1.2) +
  geom_line(
    data = t_year_fit_loc,
    aes(x = n_year, y = Estimate, group = interaction(bay, location_clean), colour = bay),
    linewidth = 0.6, alpha = 0.80
  ) +
  scale_colour_manual(values = t_bay_colors) +
  scale_x_continuous(breaks = pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  labs(
    x = "Monitoring year",
    y = "Surface water temperature (°C)",
    colour = "Bay",
    subtitle = "b) Long-term trend (evaluated at mid-season)"
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom")

t_fig_trend


# ============================================================
# PANEL c) Bay-level slopes (°C per year) at mid-season
# ============================================================

# endpoint grid per location (observed year range), at mid-season
t_slope_grid_loc <- t_fit_df %>%
  distinct(bay, location_clean, n_year.z) %>%
  group_by(bay, location_clean) %>%
  summarise(
    n_year.z_min = min(n_year.z, na.rm = TRUE),
    n_year.z_max = max(n_year.z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(n_year.z_min, n_year.z_max),
    names_to = "endpoint",
    values_to = "n_year.z"
  ) %>%
  mutate(
    j_season.z = j_ref,
    bay        = factor(bay, levels = t_bay_levels),
    n_year     = z_to_year(n_year.z)
  )

# posterior expected draws at endpoints (INCLUDING bay/location effects)
t_slope_draws_loc <- e_t_mod %>%
  add_epred_draws(
    newdata    = t_slope_grid_loc,
    re_formula = NULL
  ) %>%
  rename(epred = .epred)

# location-level slopes WITHIN each draw (Δy/Δyear)
t_draw_slopes_loc <- t_slope_draws_loc %>%
  group_by(.draw, bay, location_clean) %>%
  arrange(n_year) %>%
  summarise(
    slope = (last(epred) - first(epred)) / (last(n_year) - first(n_year)),
    .groups = "drop"
  ) %>%
  filter(is.finite(slope))

# bay-level slopes: average locations WITHIN each draw
t_draw_slopes_bay <- t_draw_slopes_loc %>%
  group_by(.draw, bay) %>%
  summarise(slope = mean(slope), .groups = "drop")

# summarise across draws
t_slope_summ_bay <- t_draw_slopes_bay %>%
  group_by(bay) %>%
  summarise(
    estimate = mean(slope),
    lower50  = quantile(slope, 0.25),
    upper50  = quantile(slope, 0.75),
    lower90  = quantile(slope, 0.05),
    upper90  = quantile(slope, 0.95),
    .groups  = "drop"
  ) %>%
  mutate(bay = factor(bay, levels = rev(t_bay_levels)))

# population slope summary (fixed effect only)
# NOTE: b_n_year.z is per 1 SD(year) -> convert to per 1 year: / year_sd
t_slope_summ_pop <- e_t_mod %>%
  spread_draws(b_n_year.z) %>%
  mutate(slope_per_year = b_n_year.z / year_sd) %>%
  summarise(
    estimate = mean(slope_per_year),
    lower50  = quantile(slope_per_year, 0.25),
    upper50  = quantile(slope_per_year, 0.75),
    lower90  = quantile(slope_per_year, 0.05),
    upper90  = quantile(slope_per_year, 0.95)
  )

t_fig_slopes <- ggplot(t_slope_summ_bay, aes(x = bay, y = estimate, colour = bay)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = t_slope_summ_pop$lower90, ymax = t_slope_summ_pop$upper90, alpha = 0.12) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = t_slope_summ_pop$lower50, ymax = t_slope_summ_pop$upper50, alpha = 0.25) +
  geom_hline(yintercept = t_slope_summ_pop$estimate, linewidth = 1.2, colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.9) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.2) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = t_bay_colors, name = "Bay") +
  scale_y_continuous(breaks = pretty_breaks(n = 6)) +
  labs(
    x = NULL,
    y = "Change in surface water temperature (°C) per year\n(evaluated at mid-season)",
    subtitle = "c)"
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

t_fig_slopes


# ============================================================
# COMBINE PANELS (ABC)
# ============================================================

# optional: hide legends on b and c and collect once
t_fig_trend_noleg  <- t_fig_trend  + guides(colour = "none") + theme(legend.position = "none")
t_fig_slopes_noleg <- t_fig_slopes + guides(colour = "none") + theme(legend.position = "none")

fig_temp_combo_abc <- (t_fig_season + t_fig_trend_noleg + t_fig_slopes_noleg) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  )

fig_temp_combo_abc
