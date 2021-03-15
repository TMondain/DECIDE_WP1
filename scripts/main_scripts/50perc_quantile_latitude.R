
### getting 50% quantiles of phenology

library(tidyverse)
library(readr)


mdf <- read_csv("data/edited_insect_data/moth/DayFlyingMoths.csv")

bdf <- read_csv("data/edited_insect_data/butterfly/BNM_NRMS_No_Duplicates_east_north.csv")