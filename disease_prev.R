


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
setwd("~/Dropbox/_Projects/PEI Oysters/Data/MSX and Dermo")
# Maddy
setwd("~/Data/MSX and Dermo/")


# load dataset you want
# mon_2013_2024 <- read.csv("Oyster Monitoring Results 2013- present.csv", header= TRUE)
disease_dat <- read.csv("DiseaseSurveillance2025_0.csv", header= TRUE)



head(disease_dat)

disease_dat %>% select(Location, `Sample.Source`, `Sampling.Season`) %>% distinct() %>% 
  arrange(Location, `Sample.Source`, `Sampling.Season`)


disease_dat %>% select(Location, `Sample.Source`, `Sampling.Season`,`Site.Status`) %>% distinct() %>% 
  arrange(Location, `Sample.Source`, `Sampling.Season`, `Site.Status`)

colnames(disease_dat)

disease_prep <- disease_dat %>%
  mutate(
    samp_source = `Sample.Source`,
    samp_seas   = factor(`Sampling.Season`,
                         levels = c("Spring", "Fall")),  # <— order here
    msx_p  = (`MSX.Prevalence....` / 100),
    dermo_p = (`Dermo.Prevalance....` / 100)
  )


head(disease_prep)


msx_p_mod <- brm( msx_p ~  samp_source * samp_seas + ( samp_source * samp_seas | Location) ,
                  data = disease_prep, iter = 10000, warmup = 1000, family = zero_one_inflated_beta(),
                  #control = list(max_treedepth = 15)
                  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

save(msx_p_mod, file = "~/Dropbox/_Projects/PEI Oysters/Model_fits/Disease/msx_p_mod.Rdata")
load("~/Dropbox/_Projects/PEI Oysters/Model_fits/Disease/msx_p_mod.Rdata")

pp_check(msx_p_mod)
conditional_effects(msx_p_mod)


msx_p_c <- conditional_effects(msx_p_mod, effects = 'samp_source:samp_seas', re_formula = NA, method = 'fitted')  # conditional effects

head(msx_p_c)

#alpha_dat$site_status <- factor(alpha_dat$site_status  , levels=c("old field","never-ploughed"))


#head(alpha_c)

pd_pts <- position_jitterdodge(
  jitter.width  = 0.08,
  jitter.height = 0.05,
  dodge.width   = 0.6
)

pd_bar <- position_dodge(width = 0.6)

msx_fig <-   ggplot() + 
  geom_point(data = disease_prep,
             aes(x = samp_source, y = msx_p, colour = samp_seas, group=samp_seas), 
             size = 2, alpha = 0.7,  position = pd_pts) +
  geom_point(data = msx_p_c$samp_source,
             aes(x = samp_source, y = estimate__, colour = samp_seas, group=samp_seas), size = 3,
             position = pd_bar) +
  geom_errorbar(data = msx_p_c$samp_source,
                aes(x = samp_source, ymin = lower__, ymax = upper__, colour = samp_seas, group=samp_seas),
                size = 1, width = 0,  position = pd_bar) +
  labs(x = '',
       y='') +
  #ylim(0,60)+
  scale_color_manual(values =  c(	 "#72aeb6", "#134b73" , "#C0C0C0"),
                     breaks = c("Spring", "Fall"),   # controls order
                     labels = c("Spring 2025", 
                                "Fall 2025"),
                     name = "Sampling Season" )  + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  #ggtitle((expression(paste(italic(alpha), '-scale', sep = ''))))+
  # scale_colour_viridis_d(option = "viridis") +
  theme_bw(base_size=18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               plot.margin= margin(t = 0.2, r = 0.2, b = -0.2, l = 0.2, unit = "cm"),
                               plot.title=element_text(size=18, hjust=0.5),
                               strip.background = element_blank(),legend.position="bottom") +
  labs( subtitle= 'a)'
  ) + ylab( "Prevelence of MSX \nin Oyster Samps (%)") 


msx_fig


# dermo
view(disease_prep)
disease_prep %>%
  summarise(
    n = n(),
    n0 = sum(dermo_p == 0, na.rm = TRUE),
    n1 = sum(dermo_p == 1, na.rm = TRUE),
    p0 = mean(dermo_p == 0, na.rm = TRUE),
    p1 = mean(dermo_p == 1, na.rm = TRUE)
  )

dermo_p_mod <- brm( dermo_p ~  samp_source * samp_seas + ( 1 | Location) ,
                    data = disease_prep, iter = 3000, warmup = 1000, family = zero_one_inflated_beta(),
                    #control = list(max_treedepth = 15)
                    control = list(adapt_delta = 0.95, max_treedepth = 12),
                    refresh = 100
)


save(dermo_p_mod, file = "~/Dropbox/_Projects/PEI Oysters/Model_fits/Disease/dermo_p_mod.Rdata")
load("~/Dropbox/_Projects/PEI Oysters/Model_fits/Disease/dermo_p_mod.Rdata")


conditional_effects(dermo_p_mod)



dermo_p_c <- conditional_effects(dermo_p_mod, effects = 'samp_source:samp_seas', re_formula = NA, method = 'fitted')  # conditional effects

head(dermo_p_c)

#alpha_dat$site_status <- factor(alpha_dat$site_status  , levels=c("old field","never-ploughed"))


#head(alpha_c)

dermo_fig <-   ggplot() + 
  geom_point(data = disease_prep,
             aes(x = samp_source, y = dermo_p, colour = "#C0C0C0", group=samp_seas), 
             size = 2, alpha = 0.7, position = pd_pts) +
  geom_point(data = dermo_p_c$samp_source,
             aes(x = samp_source, y = estimate__, colour = samp_seas, group=samp_seas), size = 3, position = pd_bar) +
  geom_errorbar(data = dermo_p_c$samp_source,
                aes(x = samp_source, ymin = lower__, ymax = upper__, colour = samp_seas, group=samp_seas),
                size = 1, width = 0, position = pd_bar) +
  labs(x = '',
       y='') +
  ylim(0,60)+
  scale_color_manual(values =  c(	"#72aeb6","#134b73" , "#C0C0C0"),
                     breaks = c("Spring", "Fall"),   # controls order
                     labels = c("Spring 2025", 
                                "Fall 2026"),
                     name = "Sampling Season" )  + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  #ggtitle((expression(paste(italic(alpha), '-scale', sep = ''))))+
  # scale_colour_viridis_d(option = "viridis") +
  theme_bw(base_size=18) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 plot.margin= margin(t = 0.2, r = 0.2, b = -0.2, l = 0.2, unit = "cm"),
                                 plot.title=element_text(size=18, hjust=0.5),
                                 strip.background = element_blank(),legend.position="none") +
  labs( subtitle= 'a)'
  ) + ylab( "Prevelence of dermo \nin 2025 Oyster Samps (%)") 


dermo_fig




# site prev

head(disease_prep)


d_dat <- disease_prep %>% #select(Location, samp_source, samp_seas, Site.Status) %>% distinct() %>%
  arrange(Location, samp_source, samp_seas, Site.Status) %>%
  mutate(status = `Site.Status`) %>% select(-`Site.Status`) %>%
  mutate(
    pres = case_when(
      status %in% c("MSX Detected", "MSX & Dermo Detected") ~ 1,
      status == "Not Detected"                              ~ 0,
      TRUE                                                  ~ NA_real_
    )
  )


head(d_dat)
site_prev_mod <- brm( pres ~  samp_source * samp_seas + ( samp_source * samp_seas | Location) ,
                      data = d_dat, iter = 10000, warmup = 1000, family = bernoulli(),
                      #control = list(max_treedepth = 15)
                      control = list(adapt_delta = 0.99, max_treedepth = 15)
)




save(site_prev_mod, file = "~/Dropbox/_Projects/PEI Oysters/Model_fits/Disease/site_prev_mod.Rdata")
load("~/Dropbox/_Projects/PEI Oysters/Model_fits/Disease/site_prev_mod.Rdata")

summary(site_prev_mod)
pp_check(site_prev_mod)
conditional_effects(site_prev_mod)



site_p_c <- conditional_effects(site_prev_mod, effects = 'samp_source:samp_seas', re_formula = NA, method = 'fitted')  # conditional effects

head(site_p_c)

#alpha_dat$site_status <- factor(alpha_dat$site_status  , levels=c("old field","never-ploughed"))


#head(alpha_c)


site_prev_fig <-   ggplot() + 
  geom_point(
    data = d_dat,
    aes(x = samp_source, y = pres, colour = samp_seas),
    size = 1, alpha = 0.7,
    position = pd_pts
  ) +
  # summary points: dodge
  geom_point(
    data = site_p_c$samp_source,
    aes(x = samp_source, y = estimate__, colour = samp_seas),
    size = 3,
    position = pd_bar
  ) +
  # error bars: same dodge
  geom_errorbar(
    data = site_p_c$samp_source,
    aes(x = samp_source, ymin = lower__, ymax = upper__, colour = samp_seas),
    size = 1,
    width = 0,
    position = pd_bar
  )+
  labs(color="Season", x = '',
       y='') +
  ylim(0,1)+
  coord_cartesian()+
  scale_color_manual(values =  c(	"#72aeb6","#134b73" , "#C0C0C0"))  + 
  #ggtitle((expression(paste(italic(alpha), '-scale', sep = ''))))+
  # scale_colour_viridis_d(option = "viridis") +
  theme_bw(base_size=18) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 plot.margin= margin(t = 0.2, r = 0.2, b = -0.2, l = 0.2, unit = "cm"),
                                 plot.title=element_text(size=18, hjust=0.5),
                                 strip.background = element_blank(),legend.position="bottom") +
  labs( subtitle= 'a)'
  ) + ylab( "Prevelence of disease at sites") 


site_prev_fig



