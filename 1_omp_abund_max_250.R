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


la_me_dat <- omp_dat %>%  # create numeric month day
  group_by(bay, location_clean, f_year) %>% # group by location and year 
  filter( larvae_above_250_microns == max(larvae_above_250_microns))  %>% # take max of larvae for each group specified above
 # filter(larvae_total > 0) %>%      
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

View(la_me_dat)
summary(la_me_dat)
nrow(la_me_dat)

# ==================================================
# MODEL
# ==================================================


# la_me_mod <- brm(
#   larvae_above_250_microns ~ water_temp.m * n_year.m * salinity.m +
#     #s(julian_date,  k = 10) +
#     (1 + n_year.m | bay/location_clean),
#   data    = la_me_dat,
#   family  = lognormal(),
#   iter    = 5000, warmup = 1000,
#   control = list(adapt_delta = 0.999, max_treedepth = 20)
# )


# ================================================

# Emma's paths
save(la_me_mod, file = "~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/la_me_mod.Rdata")
load("~/Dropbox/_Projects/PEI Oysters/Model_fits/OMP/la_me_mod.Rdata")

# Maddy's path
save(la_me_mod, file = "~/Data/Model_fits/OMP/la_me_mod.Rdata")
load("~/Data/Model_fits/OMP/la_me_mod.Rdata")


summary(la_me_mod)
pp_check(la_me_mod) + xlim(1,1500) 
conditional_effects(la_me_mod)

# ============================================================
# ALIGNMENT SETUP (shared constants)
# ============================================================

set.seed(123)

la_me_ref <- la_me_dat %>%
  ungroup() %>%
  summarise(
    wt_mean   = mean(water_temp, na.rm = TRUE),
    yr_mean   = mean(n_year,     na.rm = TRUE),
    sal_m_q25 = quantile(salinity.m, 0.25, na.rm = TRUE),
    sal_m_q50 = quantile(salinity.m, 0.50, na.rm = TRUE),
    sal_m_q75 = quantile(salinity.m, 0.75, na.rm = TRUE),
    wt_ref    = median(water_temp, na.rm = TRUE),
    #jd_ref    = median(julian_date, na.rm = TRUE)   # <-- NEW (raw julian_date)
  )

# salinity slices on centered scale (matches model input)
la_me_sal_df <- tibble(
  salinity.m = c(la_me_ref$sal_m_q25, la_me_ref$sal_m_q50, la_me_ref$sal_m_q75),
  sal_label  = factor(
    c("Low salinity (25th pct)", "Median salinity", "High salinity (75th pct)"),
    levels = c("Low salinity (25th pct)", "Median salinity", "High salinity (75th pct)")
  )
)

# ============================================================
# la_me_W) TEMP ON X — prediction grid (year-specific temp support)
# ============================================================

la_me_grid <- la_me_dat %>%
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
  crossing(la_me_sal_df) %>%
  mutate(
    water_temp.m  = water_temp - la_me_ref$wt_mean,
    n_year.m      = n_year     - la_me_ref$yr_mean,
    #julian_date = 0,     
    year_group    = factor(n_year),
    row_id        = row_number()
  )

# posterior draws (population-level)
la_me_ep_mat <- posterior_epred(
  la_me_mod,
  newdata    = la_me_grid,
  re_formula = NA
)

# draw pool (memory-safe)
la_me_pool_size <- min(2000, nrow(la_me_ep_mat))
la_me_draw_ids  <- sample(seq_len(nrow(la_me_ep_mat)), size = la_me_pool_size)

# to long + join grid meta
la_me_ep_long <- la_me_ep_mat[la_me_draw_ids, ] %>%
  as_tibble() %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols = - .draw,
    names_to = "row_id",
    values_to = "epred"
  ) %>%
  mutate(row_id = readr::parse_number(row_id)) %>%
  left_join(
    la_me_grid %>% select(row_id, n_year, water_temp, year_group, sal_label),
    by = "row_id"
  )

# ribbons: 50% + 90% + median
la_me_summ <- la_me_ep_long %>%
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
la_me_y_max_obs  <- max(la_me_dat$larvae_total, na.rm = TRUE)
la_me_y_min_plot <- max(1, min(la_me_summ$low90, na.rm = TRUE))
la_me_y_max_plot <- min(la_me_y_max_obs * 1.1, max(la_me_summ$up90, na.rm = TRUE))

# plot
la_me_fig_panel <- ggplot(
  la_me_summ,
  aes(x = water_temp, y = med, group = year_group, colour = year_group)
) +
  geom_ribbon(aes(ymin = low90, ymax = up90), alpha = 0.10, colour = NA) +
  geom_ribbon(aes(ymin = low50, ymax = up50), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ sal_label, nrow = 1) +
  scale_colour_viridis_d(option = "viridis", name = "Monitoring year" ) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30, 35)) +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(4, 8, 16, 32, 64, 128, 256, 512, 1024),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(la_me_y_min_plot, la_me_y_max_plot)) +
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

la_me_fig_panel
# ============================================================
# ltA_WMED) Median salinity only — filtered draws pipeline
#   - slopes/intercepts use ALL filtered draws (stable)
# ============================================================

la_me_med_df <- tibble(
  salinity.m = la_me_ref$sal_m_q50,
  sal_label  = factor("Median salinity", levels = "Median salinity")
)

la_me_grid_med <- la_me_dat %>%
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
  crossing(la_me_med_df) %>%
  mutate(
    water_temp.m  = water_temp - la_me_ref$wt_mean,
    n_year.m      = n_year     - la_me_ref$yr_mean,
    #julian_date   = la_me_ref$jd_ref,          # <-- REPLACE julian_date.m = 0
    year_group    = factor(n_year),
    row_id        = row_number()
  )


la_me_ep_mat_med <- posterior_epred(
  la_me_mod,
  newdata    = la_me_grid_med,
  re_formula = NA
)

la_me_pool_size_med <- min(2000, nrow(la_me_ep_mat_med))
la_me_draw_ids_med  <- sample(seq_len(nrow(la_me_ep_mat_med)), size = la_me_pool_size_med)

la_me_ep_long_med <- la_me_ep_mat_med[la_me_draw_ids_med, ] %>%
  as_tibble() %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols = - .draw,
    names_to = "row_id",
    values_to = "epred"
  ) %>%
  mutate(row_id = readr::parse_number(row_id)) %>%
  left_join(
    la_me_grid_med %>% select(row_id, n_year, water_temp, year_group, sal_label),
    by = "row_id"
  )

# envelope 10–90 at each x-point (row_id), then keep draws with >=95% inside
la_me_bounds80_med <- la_me_ep_long_med %>%
  group_by(row_id) %>%
  summarise(
    q10 = quantile(epred, 0.10, na.rm = TRUE),
    q90 = quantile(epred, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

la_me_draw_keep_med <- la_me_ep_long_med %>%
  left_join(la_me_bounds80_med, by = "row_id") %>%
  mutate(in80 = epred >= q10 & epred <= q90) %>%
  group_by(.draw) %>%
  summarise(prop_in80 = mean(in80, na.rm = TRUE), .groups = "drop") %>%
  filter(prop_in80 >= 0.95)

# ALL kept draws (for slopes/intercepts + stable ribbons/medians)
la_me_ep_keep_med <- la_me_ep_long_med %>%
  semi_join(la_me_draw_keep_med, by = ".draw")
# ============================================================
# la_me_W) TEMP ON med
# ============================================================

la_me_med <- la_me_summ %>% filter(sal_label == "Median salinity")

la_me_fig_med <- ggplot(
  la_me_med,
  aes(x = water_temp, y = med, group = year_group, colour = factor(n_year))
) +
  geom_ribbon(aes(ymin = low90, ymax = up90), alpha = 0.10, colour = NA) +
  geom_ribbon(aes(ymin = low50, ymax = up50), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) +
  #facet_wrap(~ sal_label, nrow = 1) +
  scale_colour_viridis_d(option = "viridis", name = "Monitoring year" ) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30, 35)) +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(4, 8, 16, 32, 64, 128, 256, 512, 1024),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(la_me_y_min_plot, la_me_y_max_plot)) +
  labs(
    x = "Surface water temperature (°C)",
    y = "Total oyster larvae\n >250 microns",
    #title = "Temperature–abundance relationships vary through time and with salinity",
    subtitle = "a)",
    colour= "Year"
    
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

la_me_fig_med

# ============================================================
# INTERCEPTS + SLOPES (ALIGNED: computed from ALL kept draws)
# ============================================================

# intercept-like: predicted larvae at reference temperature (median observed water temp)
la_me_temp_ref <- la_me_ref$wt_ref

la_me_intercept_draws <- la_me_ep_keep_med %>%
  group_by(.draw, n_year, year_group) %>%
  mutate(dist_to_ref = abs(water_temp - la_me_temp_ref)) %>%
  slice_min(order_by = dist_to_ref, n = 1, with_ties = FALSE) %>%
  ungroup()

la_me_intercept_summ <- la_me_intercept_draws %>%
  group_by(n_year, year_group) %>%
  summarise(
    estimate = median(epred),
    lower50  = quantile(epred, 0.25),
    upper50  = quantile(epred, 0.75),
    lower90  = quantile(epred, 0.05),
    upper90  = quantile(epred, 0.95),
    .groups = "drop"
  )

la_me_fig_intercepts <- ggplot(la_me_intercept_summ, aes(x = n_year, y = estimate, colour = factor(n_year))) +
  geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.6) +
  scale_colour_viridis_d(option = "viridis", name = "Monitoring year") +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(4, 8, 16, 32, 64, 128, 256, 512, 1024),
    labels = scales::label_number()
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  labs(
    x = "Year",
    y = "Predicted total larvae\n at reference temperature",
    #title = "Intercept-like differences through time (median salinity)",
    subtitle = "b)"
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")

la_me_fig_intercepts


# slope: within-year temperature sensitivity on log10 scale (stable for lognormal)
la_me_slope_draws <- la_me_ep_keep_med %>%
  mutate(log10_epred = log10(epred)) %>%
  group_by(.draw, n_year, year_group) %>%
  summarise(
    slope_log10 = cov(water_temp, log10_epred) / var(water_temp),
    .groups = "drop"
  )

la_me_slope_summ <- la_me_slope_draws %>%
  group_by(n_year, year_group) %>%
  summarise(
    estimate = median(slope_log10),
    lower50  = quantile(slope_log10, 0.25),
    upper50  = quantile(slope_log10, 0.75),
    lower90  = quantile(slope_log10, 0.05),
    upper90  = quantile(slope_log10, 0.95),
    .groups = "drop"
  )

la_me_fig_slopes <- ggplot(la_me_slope_summ, aes(x = n_year, y = estimate, colour = factor(n_year))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(ymin = lower90, ymax = upper90), linewidth = 0.8, alpha = 0.55) +
  geom_linerange(aes(ymin = lower50, ymax = upper50), linewidth = 2.0, alpha = 0.95) +
  geom_point(size = 2.6) +
  scale_colour_viridis_d(option = "viridis",name = "Monitoring year") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  labs(
    x = "Year",
    y = "Slope: d(log10(larvae))/d(°C)",
    #title = "Temperature sensitivity through time (median salinity)",
    subtitle = "c)" #"Computed from the same kept-draw set used for the spaghetti plot. Thick = 50% CrI; thin = 90% CrI."
  ) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

la_me_fig_slopes

# Optional combined view (patchwork)
la_me_fig_med + la_me_fig_intercepts + la_me_fig_slopes

# extra legends removed
la_me_fig_intercepts <- la_me_fig_intercepts +
  guides(colour = "none") +
  theme(legend.position = "none")

la_me_fig_slopes <- la_me_fig_slopes +
  guides(colour = "none") +
  theme(legend.position = "none")


# combined first panel (graphs a, b, c)
fig_1_combo_abc <- (la_me_fig_med + la_me_fig_intercepts + la_me_fig_slopes) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  )

fig_1_combo_abc


# ============================================================
# la_me_T) YEAR ON X × TEMP LEVELS @ MEDIAN SALINITY (aligned)
# ============================================================

la_me_T_temp_df <- la_me_dat %>%
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


la_me_T_year_seq <- seq(
  floor(min(la_me_dat$n_year, na.rm = TRUE)),
  ceiling(max(la_me_dat$n_year, na.rm = TRUE)),
  by = 1
)

la_me_T_grid <- crossing(
  n_year = la_me_T_year_seq,
  la_me_T_temp_df,
  tibble(salinity.m = la_me_ref$sal_m_q50, sal_label = factor("Median salinity"))
) %>%
  mutate(
    n_year.m      = n_year     - la_me_ref$yr_mean,
    water_temp.m  = water_temp - la_me_ref$wt_mean,
    #julian_date = 0,             # <-- REPLACE julian_date.m = 0
    row_id        = row_number()
  )

la_me_T_ep_mat <- posterior_epred(
  la_me_mod,
  newdata    = la_me_T_grid,
  re_formula = NA
)

la_me_T_pool_size <- min(2000, nrow(la_me_T_ep_mat))
la_me_T_draw_ids  <- sample(seq_len(nrow(la_me_T_ep_mat)), size = la_me_T_pool_size)

la_me_T_ep_long <- la_me_T_ep_mat[la_me_T_draw_ids, ] %>%
  as_tibble() %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(cols = - .draw, names_to = "row_id", values_to = "epred") %>%
  mutate(row_id = readr::parse_number(row_id)) %>%
  left_join(
    la_me_T_grid %>% select(row_id, n_year, temp_label, sal_label),
    by = "row_id"
  )

# unfiltered summary ribbons (50% + 90%)
la_me_T_summ <- la_me_T_ep_long %>%
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
la_me_T_y_max_obs  <- max(la_me_dat$larvae_total, na.rm = TRUE)
la_me_T_y_min_plot <- max(1, min(la_me_T_summ$low90, na.rm = TRUE))
la_me_T_y_max_plot <- min(la_me_T_y_max_obs * 1.1, max(la_me_T_summ$up90, na.rm = TRUE))

la_me_T_fig_time_temp <- ggplot(
  la_me_T_summ,
  aes(x = n_year, y = med, colour = temp_label, group = temp_label)
) +
  geom_ribbon(aes(ymin = low90, ymax = up90, fill = temp_label), alpha = 0.10, colour = NA) +
  geom_ribbon(aes(ymin = low50, ymax = up50, fill = temp_label), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.1) +
  # facet_wrap(~ sal_label, nrow = 1) +
  # scale_colour_viridis_d(option = "viridis") +
  #scale_fill_viridis_d(option = "viridis") +
  # --- Dark2 palette with explicit mapping ---
  scale_colour_manual(
    name   = "Temperature",
    values = c(
      "Cool temp (25th pct)" = "#1B9E77",
      "Median temp"          = "#7570B3",
      "Warm temp (75th pct)" = "#D95F02"
    ),
    labels = c(
      "Cool (25th pct)",
      "Median",
      "Warm (75th pct)"
    )
  ) +
  scale_fill_manual(
    name   = "Temperature",
    values = c(
      "Cool temp (25th pct)" = "#1B9E77",
      "Median temp"          = "#7570B3",
      "Warm temp (75th pct)" = "#D95F02"
    ),
    labels = c(
      "Cool (25th pct)",
      "Median",
      "Warm (75th pct)"
    )
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(4, 8, 16, 32, 64, 128, 256, 512, 1024),
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(la_me_T_y_min_plot, la_me_T_y_max_plot)) +
  labs(
    x = "Year",
    y = "Max oyster larvae\n > 250 microns",
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

la_me_T_fig_time_temp


# ============================================================
# la_me_T) INTERCEPTS + SLOPES for COOL / MEDIAN / WARM temps
# - uses the SAME la_me_T_ep_long object from the year-on-x workflow
# - tidyverse only + no custom functions + no if/else
# - uncertainty: 50% + 90% CrI everywhere
# ============================================================

# ------------------------------------------------------------
# 0) Use your existing temp levels table (must have: water_temp, temp_label)
#     and your year sequence grid already used for la_me_T_ep_long
# ------------------------------------------------------------
# Assumes these already exist from the aligned la_me_T workflow:
#   la_me_T_temp_df, la_me_T_summ, la_me_T_ep_long, la_me_ref, la_me_T_year_seq

# ------------------------------------------------------------
# 1) INTERCEPTS: predicted larvae at a reference year (mean year)
#    - within each draw × temp_label, take epred at year closest to ref
# ------------------------------------------------------------
la_me_T_ref_year <- la_me_dat %>%
  ungroup() %>%
  summarise(ref_year = mean(n_year, na.rm = TRUE)) %>%
  pull(ref_year)

la_me_T_intercept_draws <- la_me_T_ep_long %>%
  group_by(.draw, temp_label) %>%
  mutate(dist_to_ref_year = abs(n_year - la_me_T_ref_year)) %>%
  slice_min(order_by = dist_to_ref_year, n = 1, with_ties = FALSE) %>%
  ungroup()

la_me_T_intercept_summ <- la_me_T_intercept_draws %>%
  group_by(temp_label) %>%
  summarise(
    estimate = median(epred),
    lower50  = quantile(epred, 0.25),
    upper50  = quantile(epred, 0.75),
    lower90  = quantile(epred, 0.05),
    upper90  = quantile(epred, 0.95),
    .groups  = "drop"
  )

la_me_T_fig_intercepts <- ggplot(
  la_me_T_intercept_summ,
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
    y = "Predicted larvae \n at reference year",
    #title = "Baseline larvae abundance by temperature (median salinity)",
    subtitle = "e)"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x  = element_blank(),
    legend.position = "none"
  )

la_me_T_fig_intercepts


# ------------------------------------------------------------
# 2) SLOPES: change over time within each temp level
#    - compute slope within each draw × temp_label using cov/var
#    - (optional) do on log10 scale for stability/interpretability (recommended)
# ------------------------------------------------------------
la_me_T_slope_draws <- la_me_T_ep_long %>%
  mutate(log10_epred = log10(epred)) %>%
  group_by(.draw, temp_label) %>%
  summarise(
    slope_log10 = cov(n_year, log10_epred) / var(n_year),
    .groups = "drop"
  )

la_me_T_slope_summ <- la_me_T_slope_draws %>%
  group_by(temp_label) %>%
  summarise(
    estimate = median(slope_log10),
    lower50  = quantile(slope_log10, 0.25),
    upper50  = quantile(slope_log10, 0.75),
    lower90  = quantile(slope_log10, 0.05),
    upper90  = quantile(slope_log10, 0.95),
    .groups  = "drop"
  )

view(la_me_T_slope_summ)

la_me_T_fig_slopes <- ggplot(
  la_me_T_slope_summ,
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
    y = "Slope: (log10(larvae))/(year)",
    #title = "Temporal trends by temperature (median salinity)",
    subtitle = "f)"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x  = element_blank(),
    legend.position = "none"
  )

la_me_T_fig_slopes


# Optional: combine with your year-on-x ribbon plot (patchwork)
la_me_T_fig_time_temp + la_me_T_fig_intercepts + la_me_T_fig_slopes

# legends cleaned and centered

la_me_T_fig_intercepts <- la_me_T_fig_intercepts +
  guides(colour = "none") +
  theme(legend.position = "none")

la_me_T_fig_slopes <- la_me_T_fig_slopes +
  guides(colour = "none") +
  theme(legend.position = "none")

# Combined bottom panel

fig_1_combo_def <- (la_me_T_fig_time_temp + la_me_T_fig_intercepts + la_me_T_fig_slopes) +
    plot_layout(ncol = 3, guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.box.just = "center"
    )

fig_1_combo_def

# FINAL FIGURE 1 ==============

fig_1 <- (fig_1_combo_abc / fig_1_combo_def)

fig_1


