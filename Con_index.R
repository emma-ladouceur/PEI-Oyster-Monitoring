# testing out git with students!
library(tidyverse)
library(ggplot2)
library(viridis)
library(lubridate)
library(brms)
library(patchwork)


# set your own personal working directory below
# Emma MBP
setwd("~/Dropbox/_Projects/PEI Oysters/Data/OMP")
# Maddy
setwd("~/Data/Growth/")

install.packages("readxl")   # only once
library(readxl)
library(dplyr)

con_2025 <- read_excel("Growth 2025 conditional oyster data.xlsx")
con_2024 <- read_excel("Growth 2024 conditiona;.xlsx")

# load dataset you want
con_2025 <- read.csv("Growth 2025 conditional oyster data.xlsx", header= TRUE)

head(con_2025)
summary(con_2025)
view(con_2025)
