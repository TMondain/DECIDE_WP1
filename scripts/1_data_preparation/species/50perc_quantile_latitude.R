
### getting 50% quantiles of phenology

library(tidyverse)
library(readr)
library(lubridate)
library(ggridges)

# load data
# moths
mdf <- read_csv("data/edited_insect_data/moth/DayFlyingMoths.csv")

# butterflies
bdf <- read_csv("data/edited_insect_data/butterfly/BNM_NRMS_No_Duplicates_east_north.csv")


### start with moths
mdf_phenol <- mdf %>% 
  mutate(date = dmy(TO_STARTDATE),
         day = yday(date),
         species = sp_n,
         taxa = 'moth',
         lat = round(lat, 0)) %>% 
  dplyr::select(taxa, species, date, day, lat)

head(mdf_phenol)

ggplot(data = mdf_phenol, aes(x = day, y = factor(lat))) +
  geom_density_ridges() +
  xlim(50, 300) +
  theme_bw() +
  ylab('Rounded latitude') + xlab('Day of year')

summar_moth <- mdf_phenol %>% 
  group_by(taxa, species,lat) %>% 
  summarise(lwr = quantile(day, probs = 0.25),
            upr = quantile(day, probs = 0.75))


### butterflies
bdf_phenol <- bdf %>% 
  mutate(date = dmy(TO_STARTDATE),
         day = yday(date),
         species = sp_n,
         taxa = 'butterfly',
         lat = round(LATITUDE, 0)) %>% 
  dplyr::select(taxa, species, date, day, lat)

head(bdf_phenol)

ggplot(data = bdf_phenol, aes(x = day, y = factor(lat), fill = species)) +
  geom_density_ridges(show.legend = FALSE) +
  # xlim(50, 300) +
  theme_bw() +
  ylab('Rounded latitude') + xlab('Day of year')

summar_butt <- bdf_phenol %>% 
  group_by(taxa, species, lat) %>% 
  summarise(lwr = quantile(day, probs = 0.25),
            upr = quantile(day, probs = 0.75))


## combine moth and butterfly
comb_phenol <- rbind(summar_moth, summar_butt)

# write.csv(comb_phenol,
#           file = "C:/Users/thoval/OneDrive - NERC/Documents/DECIDE/WP2/outputs/25_75_flight_periods_inverts.csv")

#####     Some models     #####

# does the flight period differ across latitudes?
m1 <- glmer(day ~ lat + (1|species),
            family = 'poisson',
            data = bdf_phenol)
summary(m1)

# what about species/latitude
m2 <- glm(day ~ lat*species,
          family = 'poisson',
          data = bdf_phenol)
summary(m2)
