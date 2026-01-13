

# testing out git with students!
library(tidyverse)
library(ggplot2)
library(viridis)
library(lubridate)
library(brms)
library(patchwork)
library(tidybayes)
library(MetBrewer)



# set your own personal working directory below
# Emma MBP
setwd("~/Dropbox/_Projects/PEI Oysters/Data/OMP")
# Maddy
setwd("~/Data/OMP/")


omp_dat <- read.csv("OMP_clean_2025.csv", header= TRUE)
head(omp_dat)

# ============================================================
# MODEL 3 (diff first -> max)
# ============================================================

# ------------------------------------------------------------
# 1) FIRST event (earliest non-zero larvae)
# ------------------------------------------------------------

m3_first <- omp_dat %>%
  ungroup() %>%
  filter(larvae_total > 0) %>%
  group_by(bay, location_clean, f_year) %>%
  slice_min(julian_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(
    area, bay, location_clean, f_year, n_year,
    julian_date, larvae_total, larvae_size, larvae_above_250_microns
  ) %>%
  mutate(
    first_julian_date              = julian_date,
    first_larvae_total             = larvae_total,
    first_larvae_size              = larvae_size,
    first_larvae_above_250_microns = larvae_above_250_microns
  ) %>%
  select(-julian_date, -larvae_total, -larvae_size, -larvae_above_250_microns)

# ------------------------------------------------------------
# 2) MAX event (max larvae >250 µm)  (using omp_dat, per your working version)
# ------------------------------------------------------------

m3_max <- omp_dat %>%
  ungroup() %>%
  group_by(bay, location_clean, f_year) %>%
  filter(larvae_above_250_microns == max(larvae_above_250_microns)) %>%
  filter(larvae_above_250_microns > 0) %>%
  ungroup() %>%
  select(
    area, bay, location_clean, f_year, n_year,
    julian_date, larvae_total, larvae_size, larvae_above_250_microns
  ) %>%
  mutate(
    max_julian_date              = julian_date,
    max_larvae_total             = larvae_total,
    max_larvae_size              = larvae_size,
    max_larvae_above_250_microns = larvae_above_250_microns
  ) %>%
  select(-julian_date, -larvae_total, -larvae_size, -larvae_above_250_microns)

# ------------------------------------------------------------
# 3) JOIN + difference + center year within bay
# ------------------------------------------------------------

m3_dat <- m3_first %>%
  left_join(
    m3_max,
    by = c("area", "bay", "location_clean", "f_year", "n_year")
  ) %>%
  mutate(diff_first_last = max_julian_date - first_julian_date) %>%
  group_by(bay) %>%
  mutate(n_year.m = n_year - mean(n_year, na.rm = TRUE)) %>%
  ungroup()

# ------------------------------------------------------------
# 4) FIT MODEL (KEEP NAME)
# ------------------------------------------------------------

# oyster_first_last <- brm(
#   diff_first_last ~ n_year.m + (n_year.m | bay/location_clean),
#   data    = m3_dat,
#   iter    = 5000,
#   warmup  = 1000,
#   family  = student(), #or gaussian
#   control = list(adapt_delta = 0.99)
# )

save(oyster_first_last, file = "~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/oyster_first_last.Rdata")
load("~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/oyster_first_last.Rdata")

summary(oyster_first_last)
pp_check(oyster_first_last)
conditional_effects(oyster_first_last)
# ============================================================
# POPULATION + BAY TRENDS (NATURAL YEAR SCALE)
# Model is fit on n_year.m, but we PLOT on n_year
# ============================================================


# ------------------------------------------------------------
# 1) Get the exact data used to fit the model
#    (avoids any new-level problems)
# ------------------------------------------------------------
m3_fit_df <- model.frame(oyster_first_last)

# ------------------------------------------------------------
# 2) Fixed bay order + colour palette
# ------------------------------------------------------------
m3_bay_levels <- m3_fit_df %>%
  distinct(bay) %>%
  arrange(bay) %>%
  pull(bay)

m3_bay_colors <- setNames(
  viridis(length(m3_bay_levels), option = "plasma"),
  m3_bay_levels
)

# ------------------------------------------------------------
# 3) Compute the year centering constant ONCE
#    This is exactly what was used when n_year.m was created
# ------------------------------------------------------------
m3_year_center <- mean(m3_dat$n_year, na.rm = TRUE)

# ------------------------------------------------------------
# 4) Build prediction grid for BAY + LOCATION trends
#    • Uses ONLY real data groups
#    • Smooths across the observed n_year.m range per group
# ------------------------------------------------------------
m3_pred_grid_bay <- m3_fit_df %>%
  distinct(bay, location_clean, n_year.m) %>%
  group_by(bay, location_clean) %>%
  summarise(
    n_year.m = seq(
      min(n_year.m, na.rm = TRUE),
      max(n_year.m, na.rm = TRUE),
      length.out = 120
    ),
    .groups = "drop"
  ) %>%
  mutate(
    # convert back to NATURAL year scale for plotting
    n_year = n_year.m + m3_year_center,
    bay    = factor(bay, levels = m3_bay_levels)
  )

# ------------------------------------------------------------
# 5) Fitted values INCLUDING bay/location effects
# ------------------------------------------------------------
m3_bay_fit <- fitted(
  oyster_first_last,
  newdata    = m3_pred_grid_bay,
  re_formula = NULL,        # include bay + location effects
  probs      = c(0.10, 0.90)
) %>%
  as_tibble() %>%
  bind_cols(m3_pred_grid_bay)

# ------------------------------------------------------------
# 6) Population-level trend (no random effects)
# ------------------------------------------------------------
m3_pop_grid <- tibble(
  n_year.m = seq(
    min(m3_fit_df$n_year.m, na.rm = TRUE),
    max(m3_fit_df$n_year.m, na.rm = TRUE),
    length.out = 150
  )
) %>%
  mutate(
    n_year = n_year.m + m3_year_center
  )

m3_pop_fit <- fitted(
  oyster_first_last,
  newdata    = m3_pop_grid,
  re_formula = NA,          # population only
  probs      = c(0.10, 0.90)
) %>%
  as_tibble() %>%
  bind_cols(m3_pop_grid)

# ------------------------------------------------------------
# 7) FINAL FIGURE: NATURAL YEAR ON X-AXIS
# ------------------------------------------------------------
m3_fig_trends <- ggplot() +
  
  # --- population uncertainty band ---
  geom_ribbon(
    data = m3_pop_fit,
    aes(x = n_year, ymin = Q10, ymax = Q90),
    alpha = 0.3
  ) +
  
  # --- population mean ---
  geom_line(
    data = m3_pop_fit,
    aes(x = n_year, y = Estimate),
    colour = "black",
    linewidth = 1.2
  ) +
  
  # --- bay/location trends ---
  geom_line(
    data = m3_bay_fit,
    aes(
      x = n_year,
      y = Estimate,
      group  = interaction(bay, location_clean),
      colour = bay
    ),
    linewidth = 0.7,
    alpha = 0.85
  ) +
  
  scale_colour_manual(values = m3_bay_colors) +
  
  labs(
    x = "Monitoring year",
    y = "Difference in days (first → max)",
    colour = "Bay"
  ) +
  
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

m3_fig_trends

# ============================================================
# m3: TRENDS (bay/location) + SLOPES (BAY-LEVEL ONLY)
# - m3_fig_trends: fitted trends through time (bay/location lines + population band)
# - m3_fig_slopes: bay-level slopes (days/year), computed from posterior predictions
#                  using observed endpoint years per location, then averaged
#                  across locations WITHIN each draw (correct uncertainty)
# - Simple tidyverse + tidybayes approach (no RE-name parsing)
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tidybayes)
library(viridisLite)

# ------------------------------------------------------------
# 1) Exact data used to fit the model (avoids new-level problems)
# ------------------------------------------------------------
m3_fit_df <- model.frame(oyster_first_last)

# ------------------------------------------------------------
# 2) Fixed bay order + colour palette
# ------------------------------------------------------------
m3_bay_levels <- m3_fit_df %>%
  distinct(bay) %>%
  arrange(bay) %>%
  pull(bay)

m3_bay_colors <- setNames(
  viridisLite::viridis(length(m3_bay_levels), option = "plasma"),
  m3_bay_levels
)

# ------------------------------------------------------------
# 3) Year centering constant (must match how n_year.m was created)
# ------------------------------------------------------------
m3_year_center <- mean(m3_dat$n_year, na.rm = TRUE)

# ------------------------------------------------------------
# 4) Prediction grid for BAY + LOCATION trends
#    - real bay/location groups only
#    - smooth across observed n_year.m range per group
# ------------------------------------------------------------
m3_pred_grid_bay <- m3_fit_df %>%
  distinct(bay, location_clean, n_year.m) %>%
  group_by(bay, location_clean) %>%
  summarise(
    n_year.m = seq(
      min(n_year.m, na.rm = TRUE),
      max(n_year.m, na.rm = TRUE),
      length.out = 120
    ),
    .groups = "drop"
  ) %>%
  mutate(
    n_year = n_year.m + m3_year_center,  # natural year for plotting
    bay    = factor(bay, levels = m3_bay_levels)
  )

# ------------------------------------------------------------
# 5) Fitted values INCLUDING bay/location effects
# ------------------------------------------------------------
m3_bay_fit <- fitted(
  oyster_first_last,
  newdata    = m3_pred_grid_bay,
  re_formula = NULL,          # include bay + location effects
  probs      = c(0.10, 0.90)
) %>%
  as_tibble() %>%
  bind_cols(m3_pred_grid_bay)

# ------------------------------------------------------------
# 6) Population-level trend (no random effects)
# ------------------------------------------------------------
m3_pop_grid <- tibble(
  n_year.m = seq(
    min(m3_fit_df$n_year.m, na.rm = TRUE),
    max(m3_fit_df$n_year.m, na.rm = TRUE),
    length.out = 150
  )
) %>%
  mutate(n_year = n_year.m + m3_year_center)

m3_pop_fit <- fitted(
  oyster_first_last,
  newdata    = m3_pop_grid,
  re_formula = NA,            # population only
  probs      = c(0.10, 0.90)
) %>%
  as_tibble() %>%
  bind_cols(m3_pop_grid)

# ------------------------------------------------------------
# 7) FINAL TRENDS FIGURE
# ------------------------------------------------------------
m3_fig_trends <- ggplot() +
  geom_ribbon(
    data = m3_pop_fit,
    aes(x = n_year, ymin = Q10, ymax = Q90),
    alpha = 0.30
  ) +
  geom_line(
    data = m3_pop_fit,
    aes(x = n_year, y = Estimate),
    colour = "black",
    linewidth = 1.2
  ) +
  geom_line(
    data = m3_bay_fit,
    aes(
      x = n_year,
      y = Estimate,
      group  = interaction(bay, location_clean),
      colour = bay
    ),
    linewidth = 0.7,
    alpha = 0.85
  ) +
  scale_colour_manual(values = m3_bay_colors) +
  labs(
    x = "Monitoring year",
    y = "Difference in days (first → max)",
    colour = "Bay"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

m3_fig_trends

# ============================================================
# BAY-LEVEL SLOPES (simple + correct)
# - Compute slopes at each location using endpoint years (min/max) only
# - Predictions include bay/location effects (re_formula = NULL)
# - Average slopes across locations WITHIN each draw -> bay-level
# ============================================================

# ------------------------------------------------------------
# 8) Build endpoint-only grid per bay/location (observed years only)
# ------------------------------------------------------------
m3_slope_grid_loc <- m3_fit_df %>%
  distinct(bay, location_clean, n_year.m) %>%
  group_by(bay, location_clean) %>%
  summarise(
    n_year.m_min = min(n_year.m, na.rm = TRUE),
    n_year.m_max = max(n_year.m, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(n_year.m_min, n_year.m_max),
    names_to = "endpoint",
    values_to = "n_year.m"
  ) %>%
  mutate(
    n_year = n_year.m + m3_year_center,
    bay    = factor(bay, levels = m3_bay_levels)
  )

# ------------------------------------------------------------
# 9) Posterior expected draws at endpoints (INCLUDING bay/location effects)
# ------------------------------------------------------------
m3_slope_draws_loc <- oyster_first_last %>%
  add_epred_draws(
    newdata    = m3_slope_grid_loc,
    re_formula = NULL
  ) %>%
  rename(epred = .epred)

# ------------------------------------------------------------
# 10) Location-level slopes WITHIN each draw (Δy/Δyear using endpoints)
# ------------------------------------------------------------
m3_draw_slopes_loc <- m3_slope_draws_loc %>%
  group_by(.draw, bay, location_clean) %>%
  arrange(n_year) %>%
  summarise(
    slope = (last(epred) - first(epred)) / (last(n_year) - first(n_year)),
    .groups = "drop"
  ) %>%
  filter(is.finite(slope))

# ------------------------------------------------------------
# 11) Bay-level slopes: average locations WITHIN each draw
# ------------------------------------------------------------
m3_draw_slopes_bay <- m3_draw_slopes_loc %>%
  group_by(.draw, bay) %>%
  summarise(
    slope = mean(slope),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 12) Summarise bay slopes across draws: mean + 50% and 90% CrI
# ------------------------------------------------------------
m3_slope_summ_bay <- m3_draw_slopes_bay %>%
  group_by(bay) %>%
  summarise(
    estimate = mean(slope),
    lower50  = quantile(slope, 0.25),
    upper50  = quantile(slope, 0.75),
    lower90  = quantile(slope, 0.05),
    upper90  = quantile(slope, 0.95),
    .groups = "drop"
  ) %>%
  mutate(
    bay = factor(bay, levels = rev(m3_bay_levels))
  )

# ------------------------------------------------------------
# 13) Population-level slope summary (fixed effect only)
# ------------------------------------------------------------
m3_slope_summ_pop <- oyster_first_last %>%
  spread_draws(b_n_year.m) %>%
  summarise(
    estimate = mean(b_n_year.m),
    lower50  = quantile(b_n_year.m, 0.25),
    upper50  = quantile(b_n_year.m, 0.75),
    lower90  = quantile(b_n_year.m, 0.05),
    upper90  = quantile(b_n_year.m, 0.95)
  )

# ------------------------------------------------------------
# 14) FINAL BAY-LEVEL SLOPE FIGURE
# ------------------------------------------------------------
m3_fig_slopes <- ggplot(
  m3_slope_summ_bay,
  aes(x = bay, y = estimate, colour = bay)
) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = m3_slope_summ_pop$lower90,
    ymax = m3_slope_summ_pop$upper90,
    alpha = 0.12
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = m3_slope_summ_pop$lower50,
    ymax = m3_slope_summ_pop$upper50,
    alpha = 0.25
  ) +
  geom_hline(
    yintercept = m3_slope_summ_pop$estimate,
    linewidth  = 1.2,
    colour     = "black"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour   = "grey40"
  ) +
  geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.9) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.2) +
  geom_point(size = 2.8) +
  coord_flip() +
  scale_colour_manual(values = m3_bay_colors) +
  labs(
    x = NULL,
    y = "Slope (days per year)",
    title = "Temporal trends in first–max larval timing",
    subtitle = "Bay-level slopes computed from posterior predictions at observed endpoint years\nPoints = posterior means; thick bars = 50% CrI; thin bars = 90% CrI"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# ------------------------------------------------------------
# 15) Combine plots
# ------------------------------------------------------------
m3_fig_trends + m3_fig_slopes


