
# testing out git with students!
library(tidyverse)
library(ggplot2)
library(viridis)
library(lubridate)
library(brms)
library(patchwork)
library(tidybayes)
library(tidyverse)
library(brms)
library(tidybayes)
library(viridisLite)
library(patchwork)
library(scales)




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

e_s_mod <- brm(
  salinity ~
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

conditional_effects(e_s_mod)
pp_check(e_s_mod)

# Emma's path
save(e_s_mod, file = "~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/e_s_mod.Rdata")
load("~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/e_s_mod.Rdata")

# Maddy's path
save(e_s_mod, file = "~/Data/Model_fits/OMP/e_s_mod.Rdata")
load("~/Data/Model_fits/OMP/e_s_mod.Rdata")


# 1) exact data used in model
e_fit_df <- model.frame(e_s_mod)

# 2) bay order + colors
e_bay_levels <- e_fit_df %>%
  distinct(bay) %>%
  arrange(bay) %>%
  pull(bay)

e_bay_colors <- setNames(
  viridisLite::viridis(length(e_bay_levels), option = "plasma"),
  e_bay_levels
)

# 3) helpers to convert z -> natural scales (match YOUR scaling)
year_mu <- mean(e_dat$n_year, na.rm = TRUE)
year_sd <- sd(e_dat$n_year, na.rm = TRUE)

j_mu    <- mean(e_dat$j_season, na.rm = TRUE)
j_sd    <- sd(e_dat$j_season, na.rm = TRUE)

z_to_year <- function(z) z * year_sd + year_mu
z_to_j    <- function(z) z * j_sd + j_mu


# prediction grid: real bay/location only, year fixed at mean
e_season_grid_loc <- e_fit_df %>%
  distinct(bay, location_clean) %>%
  mutate(bay = factor(bay, levels = e_bay_levels)) %>%
  crossing(
    j_season.z = seq(min(e_fit_df$j_season.z), max(e_fit_df$j_season.z), length.out = 160),
    n_year.z   = 0
  ) %>%
  mutate(
    j_season = z_to_j(j_season.z)
  )

# location-level seasonal fits (include RE)
e_season_fit_loc <- fitted(
  e_s_mod,
  newdata    = e_season_grid_loc,
  re_formula = NULL,
  probs      = c(0.10, 0.90)
) %>%
  as_tibble() %>%
  bind_cols(e_season_grid_loc)

# population seasonal fits (no RE)
e_season_grid_pop <- tibble(
  j_season.z = seq(min(e_fit_df$j_season.z), max(e_fit_df$j_season.z), length.out = 200),
  n_year.z   = 0
) %>%
  mutate(j_season = z_to_j(j_season.z))

e_season_fit_pop <- fitted(
  e_s_mod,
  newdata    = e_season_grid_pop,
  re_formula = NA,
  probs      = c(0.10, 0.90)
) %>%
  as_tibble() %>%
  bind_cols(e_season_grid_pop)

# figure
e_fig_season <- ggplot() +
  geom_ribbon(
    data = e_season_fit_pop,
    aes(x = j_season, ymin = Q10, ymax = Q90),
    alpha = 0.30
  ) +
  geom_line(
    data = e_season_fit_pop,
    aes(x = j_season, y = Estimate),
    colour = "black", linewidth = 1.2
  ) +
  geom_line(
    data = e_season_fit_loc,
    aes(
      x = j_season, y = Estimate,
      group  = interaction(bay, location_clean),
      colour = bay
    ),
    linewidth = 0.6, alpha = 0.80
  ) +
  scale_colour_manual(values = e_bay_colors) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6)
  ) +
  labs(
    x = "Day of season (Julian day relative to first sampled day in year)",
    y = "Salinity (PSU)",
    colour = "Bay",
    subtitle = "a) Seasonal pattern (year held at mean)"
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom")

e_fig_season

j_ref <- median(e_fit_df$j_season.z, na.rm = TRUE)

# grid over years for each bay/location at mid-season
e_year_grid_loc <- e_fit_df %>%
  distinct(bay, location_clean) %>%
  crossing(
    n_year.z   = seq(min(e_fit_df$n_year.z), max(e_fit_df$n_year.z), length.out = 160),
    j_season.z = j_ref
  ) %>%
  mutate(
    bay   = factor(bay, levels = e_bay_levels),
    n_year = z_to_year(n_year.z)
  )

e_year_fit_loc <- fitted(
  e_s_mod,
  newdata    = e_year_grid_loc,
  re_formula = NULL,
  probs      = c(0.10, 0.90)
) %>%
  as_tibble() %>%
  bind_cols(e_year_grid_loc)

# population trend at mid-season
e_year_grid_pop <- tibble(
  n_year.z   = seq(min(e_fit_df$n_year.z), max(e_fit_df$n_year.z), length.out = 200),
  j_season.z = j_ref
) %>%
  mutate(n_year = z_to_year(n_year.z))

e_year_fit_pop <- fitted(
  e_s_mod,
  newdata    = e_year_grid_pop,
  re_formula = NA,
  probs      = c(0.10, 0.90)
) %>%
  as_tibble() %>%
  bind_cols(e_year_grid_pop)

e_fig_trend <- ggplot() +
  geom_ribbon(
    data = e_year_fit_pop,
    aes(x = n_year, ymin = Q10, ymax = Q90),
    alpha = 0.30
  ) +
  geom_line(
    data = e_year_fit_pop,
    aes(x = n_year, y = Estimate),
    colour = "black", linewidth = 1.2
  ) +
  geom_line(
    data = e_year_fit_loc,
    aes(
      x = n_year, y = Estimate,
      group  = interaction(bay, location_clean),
      colour = bay
    ),
    linewidth = 0.6, alpha = 0.80
  ) +
  scale_colour_manual(values = e_bay_colors) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 7)
  ) +
  labs(
    x = "Monitoring year",
    y = "Salinity (PSU)",
    colour = "Bay",
    subtitle = "b) Long-term trend (evaluated at mid-season)"
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom")

e_fig_trend


# endpoint grid per location (observed year range), at mid-season j_ref
e_slope_grid_loc <- e_fit_df %>%
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
    bay        = factor(bay, levels = e_bay_levels),
    n_year     = z_to_year(n_year.z)
  )

# draws at endpoints including RE
e_slope_draws_loc <- e_s_mod %>%
  add_epred_draws(
    newdata    = e_slope_grid_loc,
    re_formula = NULL
  ) %>%
  rename(epred = .epred)

# location slopes within draw (Δy/Δyear on natural year scale)
e_draw_slopes_loc <- e_slope_draws_loc %>%
  group_by(.draw, bay, location_clean) %>%
  arrange(n_year) %>%
  summarise(
    slope = (last(epred) - first(epred)) / (last(n_year) - first(n_year)),
    .groups = "drop"
  ) %>%
  filter(is.finite(slope))

# bay slopes: average locations within draw
e_draw_slopes_bay <- e_draw_slopes_loc %>%
  group_by(.draw, bay) %>%
  summarise(slope = mean(slope), .groups = "drop")

# summarise across draws
e_slope_summ_bay <- e_draw_slopes_bay %>%
  group_by(bay) %>%
  summarise(
    estimate = mean(slope),
    lower50  = quantile(slope, 0.25),
    upper50  = quantile(slope, 0.75),
    lower90  = quantile(slope, 0.05),
    upper90  = quantile(slope, 0.95),
    .groups  = "drop"
  ) %>%
  mutate(bay = factor(bay, levels = rev(e_bay_levels)))

# population slope: fixed effect only (convert from z-year to per-year)
# b_n_year.z is per 1 SD(year) change -> convert to per 1 year:
e_slope_summ_pop <- e_s_mod %>%
  spread_draws(b_n_year.z) %>%
  mutate(slope_per_year = b_n_year.z / year_sd) %>%
  summarise(
    estimate = mean(slope_per_year),
    lower50  = quantile(slope_per_year, 0.25),
    upper50  = quantile(slope_per_year, 0.75),
    lower90  = quantile(slope_per_year, 0.05),
    upper90  = quantile(slope_per_year, 0.95)
  )

e_fig_slopes <- ggplot(e_slope_summ_bay, aes(x = bay, y = estimate, colour = bay)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = e_slope_summ_pop$lower90, ymax = e_slope_summ_pop$upper90, alpha = 0.12) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = e_slope_summ_pop$lower50, ymax = e_slope_summ_pop$upper50, alpha = 0.25) +
  geom_hline(yintercept = e_slope_summ_pop$estimate, linewidth = 1.2, colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.9) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.2) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = e_bay_colors, name = "Bay") +
  labs(
    x = NULL,
    y = "Change in salinity (PSU) per year\n(evaluated at mid-season)",
    subtitle = "c)"
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

e_fig_slopes

# ------------------------------------------------------------
# Remove duplicate legends from panels b and c
# ------------------------------------------------------------

e_fig_trend_noleg  <- e_fig_trend  +
  guides(colour = "none") +
  theme(legend.position = "none")

e_fig_slopes_noleg <- e_fig_slopes +
  guides(colour = "none") +
  theme(legend.position = "none")

# ------------------------------------------------------------
# Combine panels
# ------------------------------------------------------------

fig_salinity_combo_abc <- (e_fig_season +
                             e_fig_trend_noleg +
                             e_fig_slopes_noleg) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  )

fig_salinity_combo_abc


