# how to clean your global environment
rm(list = ls())


library(tidyverse)
library(ggplot2)
library(viridis)
library(lubridate)
library(brms)
library(patchwork)
library(tidybayes)
library(bayesplot)




# set your own personal working directory below
# Emma MBP
setwd("~/Dropbox/_Projects/PEI Oysters/Data/OMP")
# Maddy
setwd("~/Data/OMP/")


omp_dat <- read.csv("OMP_clean_2025.csv", header= TRUE)
head(omp_dat)


# lets have a look at our work
omp_dat %>% select( location_clean, bay, location, area) %>% distinct() %>% arrange( location_clean, area)

# look at the headers (top 6 rows)
head(omp_dat)

# take the max larvae above 250 microns for every location and year
m1_max_l_dat <- omp_dat %>% group_by(bay, location_clean, f_year) %>% # group by location and year
  filter( larvae_above_250_microns == max(larvae_above_250_microns))  %>% # take max of larvae for each group specified above
  filter(larvae_above_250_microns > 0) %>%
  select(location_clean, area, bay, f_year, n_year, month_day, larvae_above_250_microns, parsed_date, julian_date, water_temp, salinity) %>% # select some columns
  # create numeric month day
  mutate(n_month_day = as.numeric(month_day)) %>%
  mutate(bay = as.factor(bay)) %>% 
  mutate(n_month_day = as.numeric(n_month_day)) %>% 
  group_by(bay) %>%
  #center continuous variables for modelling
  mutate(water_temp.m = water_temp - mean(water_temp, na.rm = TRUE),
         salinity.m = salinity - mean(salinity, na.rm = TRUE),
         n_year.m = n_year - mean(n_year, na.rm = TRUE)
  ) %>% #filter(bay != "Tracadie Bay") %>%
  mutate(bay_loc = interaction(bay, location_clean, drop = TRUE)) 

summary(m1_max_l_dat)
view(m1_max_l_dat)
# ============================================================
# MODEL 1 (MAX larvae >250 µm): 
# ============================================================


# --- Fit model (keep these names exactly as requested) ---
 max_temp_time_sal <- brm(
  julian_date ~ water_temp.m * n_year.m * salinity.m +
    (1 + n_year.m | bay/location_clean),
  data    = m1_max_l_dat,
  iter    = 5000,
  warmup  = 1000,
  family  = gaussian(),
  control = list(adapt_delta = 0.999, max_treedepth = 20)
)


# ----------------------------------------------------------


# Emma's paths
save(max_temp_time_sal, file = "~/Dropbox/_Projects/PEI Oysters/Model_Fits/OMP/max_temp_time_sal.Rdata")
load("~/Dropbox/_Projects/PEI Oysters/Model_Fits/OMP/max_temp_time_sal.Rdata")

# Maddy's path
save(max_temp_time_sal, file = "~/Data/Model_fits/OMP/max_temp_time_sal.Rdata")
load("~/Data/Model_fits/OMP/max_temp_time_sal.Rdata")

# sanity checks 
summary(max_temp_time_sal)
conditional_effects(max_temp_time_sal)

pp_check(max_temp_time_sal)

color_scheme_set("darkgray")
fig_pp_max <- pp_check(max_temp_time_sal) +   
  xlab( "Julian date") + ylab("Density") + 
  ggtitle('Date of max oyster larvae \n> 250 μm (Julian date)')+
  labs(subtitle= "b)")+
  theme_classic() +  theme( plot.title=element_text(size=18, hjust=0.5), legend.position= "none")# predicted vs. observed values
fig_pp_max

#patchwork
fig_pp_start + fig_pp_max

# ============================================================
# FIGURE WORKFLOW FOR max_temp_time_sal
# ============================================================


# ============================================================
# 0) DEFINE SALINITY SLICES (CENTERED SCALE)
# ------------------------------------------------------------
# We define three salinity conditions for visualization:
#  - Low:   −1 SD
#  - Mean:   0
#  - High:  +1 SD
#
# These are NOT new data points in salinity space,
# just representative slices for conditional visualization.
# ============================================================

m1_max_sal_sd <- sd(m1_max_l_dat$salinity.m, na.rm = TRUE)

m1_max_sal_df <- tibble(
  salinity.m = c(-m1_max_sal_sd, 0, m1_max_sal_sd),
  sal_label  = factor(
    c("Low salinity (−1 SD)", "Mean salinity", "High salinity (+1 SD)"),
    levels = c("Low salinity (−1 SD)", "Mean salinity", "High salinity (+1 SD)")
  )
)

# ============================================================
# 1) IDENTIFY YEAR LEVELS FROM THE FIT DATA
# ------------------------------------------------------------
# We explicitly pull year levels from the dataset used to fit
# the model so colour ordering is stable and reproducible.
# ============================================================

m1_max_year_levels <- m1_max_l_dat %>%
  distinct(n_year) %>%
  arrange(n_year) %>%
  pull(n_year)

# ============================================================
# 2) BUILD PREDICTION GRID *ONLY FROM OBSERVED DATA*
# ------------------------------------------------------------
# This is the key restriction you asked for.
#
# We:
#   - take only observed (n_year, water_temp) combinations
#   - cross those with salinity slices
#   - center predictors exactly as in the model
#
# No sequences. No extrapolation.
# ============================================================

m1_max_pred_grid <- m1_max_l_dat %>%
  distinct(n_year, water_temp) %>%      # observed combinations only
  crossing(m1_max_sal_df) %>%
  mutate(
    # Center predictors exactly as used in the model
    water_temp.m = water_temp - mean(m1_max_l_dat$water_temp, na.rm = TRUE),
    n_year.m     = n_year     - mean(m1_max_l_dat$n_year,     na.rm = TRUE),
    
    # SAFE factor creation (no duplicated levels)
    year_group   = factor(n_year, levels = sort(unique(n_year))),
    
    # Row ID needed to map predictions back to grid rows
    row_id       = row_number()
  )


# ============================================================
# 3) POSTERIOR EXPECTED VALUES (POPULATION-LEVEL)
# ------------------------------------------------------------
# posterior_epred():
#   - gives expected julian_date (on response scale)
#   - includes residual variance
#   - we set re_formula = NA → population-level effects only
#
# Output is a draws × rows matrix.
# ============================================================

m1_max_ep <- posterior_epred(
  max_temp_time_sal,
  newdata    = m1_max_pred_grid,
  re_formula = NA
)

# ============================================================
# 4) CONVERT POSTERIOR MATRIX → LONG FORMAT
# ------------------------------------------------------------
# This lets us:
#   - compute uncertainty using tidyverse verbs
#   - avoid apply/sapply/lapply entirely
# ============================================================

m1_max_ep_long <- as_tibble(m1_max_ep) %>%
  mutate(.draw = row_number()) %>%       # draw index
  pivot_longer(
    cols      = - .draw,
    names_to  = "row_id",
    values_to = "epred"
  ) %>%
  mutate(row_id = readr::parse_number(row_id)) %>%
  left_join(
    m1_max_pred_grid %>%
      select(row_id, water_temp, n_year, year_group, sal_label),
    by = "row_id"
  )

# ============================================================
# 5) SUMMARISE POSTERIOR:
#    MEAN + 90% CREDIBLE INTERVAL
# ------------------------------------------------------------
# We summarise ACROSS DRAWS for each prediction point.
# ============================================================

m1_max_summ <- m1_max_ep_long %>%
  group_by(row_id, water_temp, n_year, year_group, sal_label) %>%
  summarise(
    estimate = mean(epred, na.rm = TRUE),
    lower90  = quantile(epred, 0.05, na.rm = TRUE),
    upper90  = quantile(epred, 0.95, na.rm = TRUE),
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
m1_max_fig1_panel <- ggplot(
  m1_max_summ,
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
  scale_x_continuous(
    breaks = c(20, 25)
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 6),
    labels = function(x) {
      format(as.Date(x - 1, origin = "2000-01-01"), "%b %d")
    }
  ) +
  labs(
    x = "Surface water temperature (°C)",
    y = "Date of peak larval abundance (>250 μm)"
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

m1_max_fig1_panel


# ============================================================
# 7) MEAN-SALINITY PANEL ONLY (FOR MAIN TEXT FIGURE)
# ------------------------------------------------------------
# This extracts the central salinity slice for a cleaner figure.
# ============================================================

m1_max_summ_mid <- m1_max_summ %>%
  filter(sal_label == "Mean salinity")

m1_max_fig1_mid_max <- ggplot(
  m1_max_summ_mid,
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
  labs(
    x = "Surface water temperature (°C)",
    y = "Date of peak larval \nabundance (>250 μm)",
    subtitle = "a)"
  ) +
  scale_x_continuous(
    breaks = c(20, 25)
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 6),
    labels = function(x) format(as.Date(x - 1, origin = "2000-01-01"), "%b %d")
  ) +
  coord_cartesian() +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

m1_max_fig1_mid_max

# Figure caption: For the clean salinity, the model shows how the predicted timing of peak larval 
# abundance (>250 μm) varies with surface water temperature (phenology) and how this relationship
# differs among monitoring years.
# AKA (model estimates when the annual peak in larval abundance (>250 μm) occurs as a function of temp
# and shows how that temp timing relationships varies across years with avg salinity conditions)
# Salinity is fixed at the mean (0 SD) = typical (standardized) salinity conditions relative to each bay
# Shaded ribbons represent 90% Bayesian credible intervals around the mean predicted date of peak larval abundance (above 250 microns)
# (90% ribbon = 5th to 95th pct, lines = mean prediction)

# ------------------------------------------------------------
# 5) Spaghetti (mean salinity; central 80% region)
# ------------------------------------------------------------
# 
# # ------------------------------------------------------------
# 
# set.seed(123)  # reproducibility: same draws are chosen each run
# 
# # ------------------------------------------------------------
# # 5a) Restrict to mean-salinity slice
# # ------------------------------------------------------------
# # We are using the long-format posterior predictions created earlier:
# #   m1_max_ep_long has columns: .draw, row_id, epred, water_temp, year_group, sal_label, ...
# #
# # Here we keep only the mean-salinity panel, because we want spaghetti
# # for the main “mean salinity” story figure.
# # ------------------------------------------------------------
# m1_max_ep_mid_long <- m1_max_ep_long %>%
#   filter(sal_label == "Mean salinity")
# 
# # ------------------------------------------------------------
# # 5b) Sample a pool of draws (limit to 2000)
# # ------------------------------------------------------------
# # Plotting ALL posterior draws can be huge and slow.
# # We randomly select a pool of draws, up to 2000.
# #
# # Why 2000?
# #   - Enough to estimate quantiles reliably
# #   - Small enough to keep memory/plotting manageable
# # ------------------------------------------------------------
# m1_max_ci_pool_size <- min(2000, n_distinct(m1_max_ep_mid_long$.draw))
# 
# m1_max_ci_pool_ids <- m1_max_ep_mid_long %>%
#   distinct(.draw) %>%
#   slice_sample(n = m1_max_ci_pool_size) %>%  # randomly choose draw IDs
#   pull(.draw)
# 
# # Filter the mean-salinity posterior predictions down to only those draws
# m1_max_ep_ci_long <- m1_max_ep_mid_long %>%
#   filter(.draw %in% m1_max_ci_pool_ids)
# 
# # ------------------------------------------------------------
# # 5c) Define the "central 80% region" at EACH x-point
# # ------------------------------------------------------------
# # For each prediction point (row_id, i.e., each water_temp × year combination),
# # we compute:
# #   q10 = 10th percentile of epred across the pool of draws
# #   q90 = 90th percentile of epred across the pool of draws
# #
# # This creates a pointwise band that covers the "central 80%" of predictions.
# # ------------------------------------------------------------
# m1_max_bounds80 <- m1_max_ep_ci_long %>%
#   group_by(row_id) %>%
#   summarise(
#     q10 = quantile(epred, 0.10),
#     q90 = quantile(epred, 0.90),
#     .groups = "drop"
#   )
# 
# # ------------------------------------------------------------
# # 5d) Keep only draws that mostly fall inside the central band
# # ------------------------------------------------------------
# # For each draw:
# #   - mark each point as "in80" if epred is between q10 and q90 for that row_id
# #   - compute prop_in80 = proportion of points within the band
# #   - keep draws with prop_in80 >= 0.95
# #
# # Interpretation:
# #   These are draws whose trajectories are not "globally extreme".
# #   They still vary, but they avoid rare huge shifts that dominate plots.
# # ------------------------------------------------------------
# m1_max_draw_keep <- m1_max_ep_ci_long %>%
#   left_join(m1_max_bounds80, by = "row_id") %>%
#   mutate(
#     in80 = epred >= q10 & epred <= q90
#   ) %>%
#   group_by(.draw) %>%
#   summarise(
#     prop_in80 = mean(in80),
#     .groups = "drop"
#   ) %>%
#   filter(prop_in80 >= 0.95)
# 
# 
# # ------------------------------------------------------------
# # 5e) Sample a fixed number of "good" draws for plotting spaghetti
# # ------------------------------------------------------------
# 
# m1_max_n_spaghetti <- 150  # <-- add this
# 
# # how many draws are actually available to sample from?
# m1_max_n_keep <- min(m1_max_n_spaghetti, nrow(m1_max_draw_keep))
# 
# m1_max_spaghetti_draw_ids <- m1_max_draw_keep %>%
#   slice_sample(n = m1_max_n_keep) %>%
#   pull(.draw)
# 
# 
# 
# # Pull the long epred values for only those spaghetti draws
# m1_max_ep_spag_long <- m1_max_ep_ci_long %>%
#   filter(.draw %in% m1_max_spaghetti_draw_ids)
# 
# # ------------------------------------------------------------
# # 5f) Build the spaghetti plot
# # ------------------------------------------------------------
# # We plot:
# #   - faint spaghetti lines: each posterior draw × year trajectory
# #   - thicker posterior mean line per year (from m1_max_summ_mid)
# #
# # Note:
# #   group = interaction(.draw, year_group) ensures each draw-year is a separate line
# #   colour = year_group keeps year colouring consistent with earlier figures
# # ------------------------------------------------------------
# m1_max_fig1_spag <- ggplot() +
#   geom_line(
#     data = m1_max_ep_spag_long,
#     aes(
#       x = water_temp,
#       y = epred,
#       group = interaction(.draw, year_group),
#       colour = year_group
#     ),
#     linewidth = 0.4,
#     alpha = 0.15
#   ) +
#   geom_line(
#     data = m1_max_summ_mid,
#     aes(
#       x = water_temp,
#       y = estimate,
#       group = year_group,
#       colour = year_group
#     ),
#     linewidth = 0.7
#   ) +
#   scale_colour_viridis_d(option = "viridis") +
#   scale_x_continuous(
#     breaks = c(20, 25)
#   ) +
#   scale_y_continuous(
#     breaks = scales::pretty_breaks(n = 6),
#     labels = function(x) format(as.Date(x - 1, origin = "2000-01-01"), "%b %d")
#   ) +
#   labs(
#     color= "Monitoring Year",
#     x = "Surface water temperature (°C)",
#     y = "Date of peak larval \nabundance (>250 μm)",
#     subtitle = "a)"
#   ) +
#   theme_bw(base_size = 18) +
#   theme(
#     panel.grid.minor = element_blank(), legend.position = "bottom"
#   )
# 
# m1_max_fig1_spag
# 

# ------------------------------------------------------------
# 6) Intercepts-by-year (temperature-averaged)
#    3 salinity panels + mean-salinity-only
# ------------------------------------------------------------
# GOAL:
#   For each monitoring year, estimate the "temperature-averaged"
#   phenology timing (Julian date) at each salinity slice.
#
# INTERPRETATION:
#   This is like a year-specific intercept, but averaged across the
#   observed temperature distribution in that year (NOT extrapolated).
#
# KEY REQUIREMENT (your request):
#   Use ONLY the observed temperatures for each year from the dataset
#   used to fit the model (m1_max_l_dat). No artificial temp sequences.
#
# UNCERTAINTY:
#   We do the temperature-average within each posterior draw first,
#   then summarise across draws → correct propagation of uncertainty.
# ------------------------------------------------------------

set.seed(123)

# ------------------------------------------------------------
# 6a) Build a year × temperature grid ONLY from observed data
# ------------------------------------------------------------
# We take distinct observed (n_year, water_temp) pairs from the fit data,
# then cross with salinity slices. This ensures:
#   - temperature averaging is based on observed temps per year
#   - no extrapolation to unobserved temperature values
# ------------------------------------------------------------
m1_max_year_sal_temp_grid <- m1_max_l_dat %>%
  distinct(n_year, water_temp) %>%      # <-- observed year-temp combos only
  crossing(m1_max_sal_df) %>%           # <-- add salinity slices
  mutate(
    # center exactly as in model
    n_year.m     = n_year     - mean(m1_max_l_dat$n_year,     na.rm = TRUE),
    water_temp.m = water_temp - mean(m1_max_l_dat$water_temp, na.rm = TRUE),
    
    # stable row id for joining back after pivoting
    row_id       = row_number()
  )

# ------------------------------------------------------------
# 6b) Posterior expected values for each grid row (population-level)
# ------------------------------------------------------------
m1_max_ep2 <- posterior_epred(
  max_temp_time_sal,
  newdata    = m1_max_year_sal_temp_grid,
  re_formula = NA
)

# ------------------------------------------------------------
# 6c) Convert draws × rows matrix → long tibble
# ------------------------------------------------------------
# This gives one row per (.draw, row_id) with epred and the grid metadata.
# We use readr::parse_number() to convert "V123" column names to row_id.
# ------------------------------------------------------------
m1_max_ep2_long <- as_tibble(m1_max_ep2) %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols      = - .draw,
    names_to  = "row_id",
    values_to = "epred"
  ) %>%
  mutate(row_id = readr::parse_number(row_id)) %>%
  left_join(
    m1_max_year_sal_temp_grid %>% select(row_id, n_year, sal_label, water_temp),
    by = "row_id"
  )

# ------------------------------------------------------------
# 6d) Temperature-average WITHIN EACH DRAW
# ------------------------------------------------------------
# Critical point:
#   Do NOT average after summarising across draws, because that would
#   understate uncertainty.
#
# Instead:
#   For each draw, year, and salinity slice:
#     compute mean epred across observed water_temp values
# ------------------------------------------------------------
m1_max_draw_means <- m1_max_ep2_long %>%
  group_by(.draw, n_year, sal_label) %>%
  summarise(
    epred_mean = mean(epred, na.rm = TRUE),   # mean over observed temps in that year
    .groups = "drop"
  )

# ------------------------------------------------------------
# 6e) Summarise across draws: mean + 50% and 90% CrI
# ------------------------------------------------------------
m1_max_int_summ <- m1_max_draw_means %>%
  group_by(n_year, sal_label) %>%
  summarise(
    estimate = mean(epred_mean, na.rm = TRUE),
    lower50  = quantile(epred_mean, 0.25, na.rm = TRUE),
    upper50  = quantile(epred_mean, 0.75, na.rm = TRUE),
    lower90  = quantile(epred_mean, 0.05, na.rm = TRUE),
    upper90  = quantile(epred_mean, 0.95, na.rm = TRUE),
    .groups  = "drop"
  )
View(m1_max_int_summ)

# ------------------------------------------------------------
# 6f) Figure: 3 salinity panels
# ------------------------------------------------------------
m1_max_fig_intercepts_3panel <- ggplot(
  m1_max_int_summ,
  aes(x = n_year, y = estimate, colour = factor(n_year), group  = factor(n_year))
) +
  # 90% interval (thin)
  geom_linerange(aes(ymin = lower90, ymax = upper90),
                 linewidth = 0.8, alpha = 0.55) +
  # 50% interval (thick)
  geom_linerange(aes(ymin = lower50, ymax = upper50),
                 linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.4) +
  facet_wrap(~ sal_label, nrow = 1) +
  scale_colour_viridis_d(option = "viridis", name = "Monitoring year") +
  scale_x_continuous(breaks = sort(unique(m1_max_l_dat$n_year))[seq(1, length(unique(m1_max_l_dat$n_year)), by = 2)]) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 6),
    labels = function(x) format(as.Date(x - 1, origin = "2000-01-01"), "%b %d")
  ) +
  labs(
    x = "Year",
    y = "Date of peak larval \nabundance (>250 μm)",
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

m1_max_fig_intercepts_3panel

# ------------------------------------------------------------
# 6g) Mean-salinity-only intercept figure (styled to pair with mid panel)
# ------------------------------------------------------------
m1_max_fig_intercepts_mean <- m1_max_int_summ %>%
  filter(sal_label == "Mean salinity") %>%
  mutate(
    # keep year ordering stable for colouring
    year_group = factor(n_year, levels = sort(unique(m1_max_l_dat$n_year)))
  ) %>%
  ggplot(aes(x = n_year, y = estimate, colour = year_group)) +
  geom_linerange(aes(ymin = lower90, ymax = upper90, colour = year_group),
                 linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50, colour = year_group),
                 linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.4) +
  scale_colour_viridis_d(option = "viridis", name = "Monitoring year") +
  scale_x_continuous(breaks = sort(unique(m1_max_l_dat$n_year))) +
  labs(
    x = "Year",
    y = "Date of peak larval \nabundance (>250 μm)",
    subtitle = "b)"
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 6),
    labels = function(x) format(as.Date(x - 1, origin = "2000-01-01"), "%b %d")
  ) +
  scale_x_continuous(
    breaks = c(2012, 2014, 2016, 2018, 2020, 2022, 2024),
    labels = c("2012", "2014", "2016", "2018", "2020", "2022", "2024")
  ) +
  coord_cartesian() +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

m1_max_fig_intercepts_mean

# Figure caption:



#slope
# ============================================================
# 7a) Draw-wise slopes per year and salinity
# ------------------------------------------------------------
# Within each (.draw, year, salinity slice), compute slope of
# epred vs water_temp using covariance/variance (OLS slope).
# ============================================================

m1_max_draw_slopes <- m1_max_ep2_long %>%
  select(.draw, n_year, sal_label, water_temp, epred) %>%
  group_by(.draw, n_year, sal_label) %>%
  summarise(
    n_pts     = n_distinct(water_temp),
    wt_var    = var(water_temp, na.rm = TRUE),
    wt_mean   = mean(water_temp, na.rm = TRUE),
    yy_mean   = mean(epred, na.rm = TRUE),
    wt_yy_cov = mean((water_temp - wt_mean) * (epred - yy_mean), na.rm = TRUE),
    slope     = wt_yy_cov / wt_var,
    .groups   = "drop"
  ) %>%
  filter(n_pts >= 2, wt_var > 0, is.finite(slope)) %>%
  select(.draw, n_year, sal_label, slope)

# ============================================================
# 7b) Summarise slopes across posterior draws
# ============================================================

m1_max_slope_summ <- m1_max_draw_slopes %>%
  group_by(n_year, sal_label) %>%
  summarise(
    estimate = mean(slope),
    lower50  = quantile(slope, 0.25),
    upper50  = quantile(slope, 0.75),
    lower90  = quantile(slope, 0.05),
    upper90  = quantile(slope, 0.95),
    .groups  = "drop"
  )


# ============================================================
# 7c) Slope figure: 3 salinity panels
# ============================================================

m1_max_fig_slopes_3panel <- ggplot(
  m1_max_slope_summ,
  aes(x = n_year, y = estimate, colour = factor(n_year), group = factor(n_year))
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.5) +
  geom_linerange(aes(ymin = lower90, ymax = upper90),
                 linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50),
                 linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.4) +
  facet_wrap(~ sal_label, nrow = 1) +
  scale_colour_viridis_d(option = "viridis", name = "Monitoring year") +
  scale_x_continuous(
    breaks = sort(unique(m1_max_l_dat$n_year))[
      seq(1, length(unique(m1_max_l_dat$n_year)), by = 2)
    ]
  ) +
  labs(
    x = "Year",
    y = "Change in date of peak larval \nabundance (>250 μm) per 1°C",
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

m1_max_fig_slopes_3panel

# ============================================================
# 7d) Mean-salinity-only slope figure
# ============================================================

m1_max_fig_slopes_mean <- m1_max_slope_summ %>%
  filter(sal_label == "Mean salinity") %>%
  mutate(year_group = factor(n_year, levels = sort(unique(m1_max_l_dat$n_year)))) %>%
  ggplot(aes(x = n_year, y = estimate, colour = year_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.5) +
  geom_linerange(aes(ymin = lower90, ymax = upper90),
                 linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50),
                 linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.4) +
  scale_colour_viridis_d(option = "viridis", name = "Monitoring year") +
  scale_x_continuous(
    breaks = c(2012, 2014, 2016, 2018, 2020, 2022, 2024),
    labels = c("2012", "2014", "2016", "2018", "2020", "2022", "2024")
  ) +
  labs(
    x = "Year",
    y = "Change in date of peak larval \nabundance (>250 μm) per 1°C",
    subtitle = "c)"
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

m1_max_fig_slopes_mean

# Figure caption:



# killing the legends in graph b and c
m1_max_fig_intercepts_mean <- m1_max_fig_intercepts_mean +
guides(colour = "none", fill = "none") +
  theme(legend.position = "none")

m1_max_fig_slopes_mean <- m1_max_fig_slopes_mean +
  guides(colour = "none", fill = "none") +
  theme(legend.position = "none")


# combo graph abc
fig_5_combo_abc <- (m1_max_fig1_mid_max +
                      m1_max_fig_intercepts_mean + m1_max_fig_slopes_mean
                    ) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center")

fig_5_combo_abc


# ============================================================
# ALT spaghetti fig: year × temperature (salinity held at mean)
# Namespaced to avoid collisions with earlier objects
# ============================================================

set.seed(123)

# --- references for centering + temp quantiles (UNCENTERED vars) ---
m1_max_alt_ref <- m1_max_l_dat %>%
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
m1_max_alt_temp_levels <- tibble(
  temp_level = factor(c("Cool", "Median", "Warm"),
                      levels = c("Cool", "Median", "Warm")),
  water_temp = c(m1_max_alt_ref$temp_cool,
                 m1_max_alt_ref$temp_med,
                 m1_max_alt_ref$temp_warm)
)

# --- prediction grid (year sequence × temp levels), salinity fixed at mean ---
m1_max_alt_newdat_year_temp <- crossing(
  n_year = seq(min(m1_max_l_dat$n_year, na.rm = TRUE),
               max(m1_max_l_dat$n_year, na.rm = TRUE),
               length.out = 200),
  temp_level = m1_max_alt_temp_levels$temp_level
) %>%
  left_join(m1_max_alt_temp_levels, by = "temp_level") %>%
  mutate(
    salinity     = m1_max_alt_ref$sal_mean,
    salinity.m   = salinity   - m1_max_alt_ref$sal_mean,
    n_year.m     = n_year     - m1_max_alt_ref$n_year_mean,
    water_temp.m = water_temp - m1_max_alt_ref$water_temp_mean
  )

# --- posterior expected draws (population-level) ---
m1_max_alt_draws <- max_temp_time_sal %>%
  add_epred_draws(newdata = m1_max_alt_newdat_year_temp, re_formula = NA) %>%
  ungroup() %>%
  rename(m1_max_alt_epred = .epred)

# --- summarise across draws (median + 80% interval) ---
m1_max_alt_summ <- m1_max_alt_draws %>%
  group_by(n_year, temp_level) %>%
  median_qi(m1_max_alt_epred, .width = 0.80) %>%
  ungroup()

# --- thin spaghetti: sample N draws total (not per temp; consistent across temps) ---
m1_max_alt_n_keep <- 60

m1_max_alt_keep_draws <- m1_max_alt_draws %>%
  distinct(.draw) %>%
  slice_sample(n = min(m1_max_alt_n_keep, nrow(.)))

m1_max_alt_draws_thin <- m1_max_alt_draws %>%
  semi_join(m1_max_alt_keep_draws, by = ".draw")

# --- plot ---
m1_max_alt_fig_year3temp_spag_thin <- ggplot() +
  geom_line(
    data = m1_max_alt_draws_thin,
    aes(
      x = n_year, y = m1_max_alt_epred,
      group = interaction(.draw, temp_level),
      colour = temp_level
    ),
    linewidth = 0.35,
    alpha = 0.10
  ) +
  geom_line(
    data = m1_max_alt_summ,
    aes(
      x = n_year, y = m1_max_alt_epred,
      group = temp_level,
      colour = temp_level
    ),
    linewidth = 1.0
  ) +
  scale_colour_manual(
    name = "Surface water \ntemperature (°C)",
    values = c(
      "Cool"   = "#1B9E77",
      "Median" = "#7570B3",
      "Warm"   = "#D95F02"
    ),
    labels = c(
      "Cool",
      "Median",
      "Warm"
    )
  ) +
  coord_cartesian(
    ylim = range(m1_max_alt_draws_thin$m1_max_alt_epred, na.rm = TRUE) + c(-2, 2)
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 3),
    labels = function(x) format(as.Date(x - 1, origin = "2000-01-01"), "%b %d")
  ) +
  scale_x_continuous(
    breaks = seq(
      from = 2016,
      to   = max(m1_max_l_dat$n_year, na.rm = TRUE),
      by   = 2
    ),
    labels = scales::label_number(accuracy = 1),
    limits = c(
      2016,
      max(m1_max_l_dat$n_year, na.rm = TRUE)
    )
  ) +
  labs(
    x = "Year",
    y = "Date of peak larval \nabundance (>250 μm)",
    subtitle = "d)"
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

m1_max_alt_fig_year3temp_spag_thin

# Figure caption:


# ============================================================
# m1 MAX: Ribbon fig + slopes + intercepts (matched to m2 START style)
# For: max_temp_time_sal (data = m1_max_l_dat)
# - no custom functions
# - tidyverse + tidybayes
# - ribbons = posterior uncertainty (80% + 50%)
# ============================================================

set.seed(123)

# ------------------------------------------------------------
# m1_max_01) References for centering + temp quantiles (UNCENTERED vars)
# ------------------------------------------------------------
m1_max_alt_ref <- m1_max_l_dat %>%
  ungroup() %>%
  summarise(
    n_year_mean     = mean(n_year, na.rm = TRUE),
    water_temp_mean = mean(water_temp, na.rm = TRUE),
    sal_mean        = mean(salinity, na.rm = TRUE),
    temp_cool       = quantile(water_temp, 0.10, na.rm = TRUE),
    temp_med        = quantile(water_temp, 0.50, na.rm = TRUE),
    temp_warm       = quantile(water_temp, 0.90, na.rm = TRUE)
  )

# ------------------------------------------------------------
# m1_max_02) Temperature level lookup table
# ------------------------------------------------------------
m1_max_alt_temp_levels <- tibble(
  temp_level = factor(c("Cool", "Median", "Warm"),
                      levels = c("Cool", "Median", "Warm")),
  water_temp = c(m1_max_alt_ref$temp_cool,
                 m1_max_alt_ref$temp_med,
                 m1_max_alt_ref$temp_warm)
)

# ------------------------------------------------------------
# m1_max_03) Prediction grid (year sequence × temp levels), salinity fixed at mean
# ------------------------------------------------------------
m1_max_alt_newdat_year_temp <- crossing(
  n_year = seq(min(m1_max_l_dat$n_year, na.rm = TRUE),
               max(m1_max_l_dat$n_year, na.rm = TRUE),
               length.out = 200),
  temp_level = m1_max_alt_temp_levels$temp_level
) %>%
  left_join(m1_max_alt_temp_levels, by = "temp_level") %>%
  mutate(
    salinity     = m1_max_alt_ref$sal_mean,
    salinity.m   = salinity   - m1_max_alt_ref$sal_mean,
    n_year.m     = n_year     - m1_max_alt_ref$n_year_mean,
    water_temp.m = water_temp - m1_max_alt_ref$water_temp_mean
  )

# ------------------------------------------------------------
# m1_max_04) Posterior expected draws (population-level)
# ------------------------------------------------------------
m1_max_alt_draws <- max_temp_time_sal %>%
  add_epred_draws(newdata = m1_max_alt_newdat_year_temp, re_formula = NA) %>%
  ungroup() %>%
  rename(m1_max_alt_epred = .epred)

# ------------------------------------------------------------
# m1_max_05) Summaries for ribbons:
#   90% interval for wide ribbon + 50% interval for inner ribbon + median line
# ------------------------------------------------------------
m1_max_alt_summ90 <- m1_max_alt_draws %>%
  group_by(n_year, temp_level) %>%
  median_qi(m1_max_alt_epred, .width = 0.90) %>%
  ungroup() %>%
  rename(
    m1_max_alt_med   = m1_max_alt_epred,
    m1_max_alt_low90 = .lower,
    m1_max_alt_up90  = .upper
  )

m1_max_alt_summ50 <- m1_max_alt_draws %>%
  group_by(n_year, temp_level) %>%
  median_qi(m1_max_alt_epred, .width = 0.50) %>%
  ungroup() %>%
  rename(
    m1_max_alt_low50 = .lower,
    m1_max_alt_up50  = .upper
  )

m1_max_alt_ribbon_dat <- m1_max_alt_summ90 %>%
  left_join(
    m1_max_alt_summ50 %>%
      select(n_year, temp_level, m1_max_alt_low50, m1_max_alt_up50),
    by = c("n_year", "temp_level")
  )

# ------------------------------------------------------------
# m1_max_06) Ribbon figure (year × temp; salinity held at mean)
# ------------------------------------------------------------
m1_max_alt_fig_year3temp_ribbons <- ggplot(m1_max_alt_ribbon_dat) +
  geom_ribbon(
    aes(x = n_year, ymin = m1_max_alt_low90, ymax = m1_max_alt_up90, fill = temp_level),
    alpha = 0.18,
    colour = NA
  ) +
  geom_ribbon(
    aes(x = n_year, ymin = m1_max_alt_low50, ymax = m1_max_alt_up50, fill = temp_level),
    alpha = 0.32,
    colour = NA
  ) +
  # median line
  geom_line(
    aes(x = n_year, y = m1_max_alt_med, colour = temp_level),
    linewidth = 1.0
  ) +
  scale_colour_manual(
    name = "Surface water \ntemperature (°C)",
    values = c(
      "Cool"   = "#1B9E77",  # green (Dark2)
      "Median" = "#7570B3",  # purple (Dark2)
      "Warm"   = "#D95F02"   # orange (Dark2)
    ),
    labels = c(
      "Cool",
      "Median",
      "Warm"
    ) 
  ) +
  scale_fill_manual(
    name = "Surface water \ntemperature (°C)",
    values = c(
      "Cool"   = "#1B9E77",
      "Median" = "#7570B3",
      "Warm"   = "#D95F02"
    ),
    labels = c(
      "Cool",
      "Median",
      "Warm"
    )
  ) +
  # scale_x_continuous(
  #   breaks = seq(
  #     from = 2016,
  #     to   = max(m1_max_l_dat$n_year, na.rm = TRUE),
  #     by   = 2
  #   ),
  #   labels = scales::label_number(accuracy = 1),
  #   limits = c(
  #     2016,
  #     max(m1_max_l_dat$n_year, na.rm = TRUE)
  #   )
  # ) +
  coord_cartesian(
    ylim = quantile(
      m1_max_l_dat$julian_date,
      probs = c(0.01, 0.99),
      na.rm = TRUE
    )
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 3),
    labels = function(x) format(as.Date(x - 1, origin = "2000-01-01"), "%b %d")
  ) +
  scale_x_continuous(
    breaks = c(2012, 2014, 2016, 2018, 2020, 2022, 2024),
    labels = c("2012", "2014", "2016", "2018", "2020", "2022", "2024")
  ) +
  labs(
    x = "Year",
    y = "Date of peak larval \nabundance (>250 μm)",
    subtitle = "d)"
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

m1_max_alt_fig_year3temp_ribbons

# Figure caption: 


# ============================================================
# m1 MAX: SLOPE FIGURE (matched to ribbon figure)
# - draw-by-draw slopes within temp_level
# - summarises across draws (mean + 50%/90% CrI)
# ============================================================

m1_max_slope_draws <- m1_max_alt_draws %>%
  group_by(.draw, temp_level) %>%
  summarise(
    slope = cov(n_year, m1_max_alt_epred) / var(n_year),
    .groups = "drop"
  )

m1_max_slope_summ <- m1_max_slope_draws %>%
  group_by(temp_level) %>%
  summarise(
    estimate = mean(slope),
    lower50  = quantile(slope, 0.25),
    upper50  = quantile(slope, 0.75),
    lower90  = quantile(slope, 0.05),
    upper90  = quantile(slope, 0.95),
    .groups  = "drop"
  )


m1_max_fig_year3temp_slopes <- ggplot(
  m1_max_slope_summ,
  aes(x = temp_level, y = estimate, colour = temp_level)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_linerange(
    aes(ymin = lower90, ymax = upper90),
    linewidth = 0.9,
    alpha = 0.6
  ) +
  geom_linerange(
    aes(ymin = lower50, ymax = upper50),
    linewidth = 2.2
  ) +
  geom_point(size = 3) +
  scale_colour_manual(
    values = c(
      "Cool"   = "#1B9E77",
      "Median" = "#7570B3",
      "Warm"   = "#D95F02"
    )
  ) +
  scale_x_discrete(
    limits = c("Cool", "Median", "Warm")
  ) +
  scale_y_continuous(
    breaks = c(-2, -1, 0),
    limits = c(
      min(m1_max_slope_summ$lower90, -2, na.rm = TRUE),
      0
    )
  ) +
  labs(
    x = NULL,
    y = "Change in date of peak larval \nabundance (>250 μm) per year",
    subtitle = "f)"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_line(),
    legend.position = "none"
  )

m1_max_fig_year3temp_slopes



# ============================================================
# m1 MAX: INTERCEPT FIGURE (matched to ribbons + slopes style)
# - intercept = predicted value at reference year (mean year)
# - summarises across draws (mean + 50%/90% CrI)
# ============================================================

m1_max_ref_year <- m1_max_alt_ref$n_year_mean

m1_max_intercept_draws <- m1_max_alt_draws %>%
  mutate(.dist = abs(n_year - m1_max_ref_year)) %>%
  group_by(.draw, temp_level) %>%
  slice_min(.dist, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    .draw,
    temp_level,
    intercept = m1_max_alt_epred
  )

m1_max_intercept_summ <- m1_max_intercept_draws %>%
  group_by(temp_level) %>%
  summarise(
    estimate = mean(intercept),
    lower50  = quantile(intercept, 0.25),
    upper50  = quantile(intercept, 0.75),
    lower90  = quantile(intercept, 0.05),
    upper90  = quantile(intercept, 0.95),
    .groups  = "drop"
  )


m1_max_fig_year3temp_intercepts <- ggplot(
  m1_max_intercept_summ,
  aes(x = temp_level, y = estimate, colour = temp_level)
) +
  geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.9, alpha = 0.6) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.2) +
  geom_point(size = 3) +
  scale_colour_manual(
    values = c(
      "Cool"   = "#1B9E77",
      "Median" = "#7570B3",
      "Warm"   = "#D95F02"
    )
  ) +
  scale_x_discrete(limits = c("Cool", "Median", "Warm")) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = function(x) format(as.Date(x - 1, origin = "2000-01-01"), "%b %d")
  ) +
  labs(
    x = NULL,
    y = "Date of peak larval \nabundance (>250 μm)",
    subtitle = "e)"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_line(),
    legend.position = "none"
  )

m1_max_fig_year3temp_intercepts


 
 
 
 # Kill legends in the other two panels
 m1_max_fig_year3temp_intercepts <- m1_max_fig_year3temp_intercepts +
   guides(colour = "none", fill = "none") +
   theme(legend.position = "none")
 
 m1_max_fig_year3temp_slopes <- m1_max_fig_year3temp_slopes +
   guides(colour = "none", fill = "none") +
   theme(legend.position = "none")
 
 # Combine and center the legend
 fig_5_combo_def <- (m1_max_alt_fig_year3temp_ribbons +
     m1_max_fig_year3temp_intercepts +
     m1_max_fig_year3temp_slopes) +
   plot_layout(ncol = 3, guides = "collect") &
   theme(
     legend.position = "bottom",
     legend.justification = "center",
     legend.box = "horizontal"
   )

 fig_5_combo_def
 
# FINAL FIG 5
 
fig_5 <- (fig_5_combo_abc / fig_5_combo_def)

fig_5

