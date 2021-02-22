## Header ----
## Script name: 03_spring_emergence.R
## Purpose of script: Estimate degree-days required for springtime oviposition
## Author: Tyson Wepprich
## Date Created: 2021-02-04
## License: GNU GPLv3
## Email: tyson.wepprich@gmail.com
## ---
## Notes:
## Script makes Figure A1 to show overwinter emergence and oviposition timing
## Both psyllid biotypes tested at 14C and 21C constant after winter diapause
## Number of eggs is key response variable, but bottom of leaf in petri dish
## unobservable as frequently as top, so total eggs per day estimated
## ---

library(dplyr)
library(ggplot2)
library(lubridate)
library(viridis)
library(magrittr)

dat <- read.csv("data/springtime_oviposition.csv", header = TRUE)

# what to do about unknown starting # of adults in each tray?
# we know the # moved and # dead when dishes transferred for 14C trays on 8/2,
# we know the # moved when dishes transferred for 21C trays on 7/23
start14 <- dat %>% 
  filter(Date == as.Date("2018-08-02"), Temperature == 14) %>% 
  group_by(Population, Photoperiod) %>% 
  summarise(StartCount = Dead[1] + StartCount[2])
# assign to correct places in original dataset
dat[which(dat$Temperature == 14 & dat$EggAccumStartDate == as.Date("2018-07-12")), ] %<>%
  mutate(StartCount = rep(start14$StartCount, 9))
start21 <- start14 %>% 
  group_by(Population) %>% 
  mutate(StartCount = round(mean(StartCount)))
dat[which(dat$Temperature == 21 & dat$EggAccumStartDate == as.Date("2018-07-12")), ] %<>%
  mutate(StartCount = rep(start21$StartCount, 4))


# Interpolate # eggs each day ----
# sometimes no adults at end because they all died, but can't have 0 because that messes up division
# estimate eggs on bottom based on cumulative proportion of top eggs
# new eggs can be negative b/c egg det prob means counts can go up or down

dat <- dat %>% 
  group_by(Population, Temperature, Photoperiod, EggAccumStartDate) %>% 
  arrange(Date) %>%
  mutate(maxtop = max(Eggs_top, na.rm = TRUE),
         maxbot = max(Eggs_bottom, na.rm = TRUE),
         CumulPropEggs = ifelse(maxtop > 0,
                                (Eggs_top / maxtop),
                                0),
         Eggs_bottom_imputed = ifelse(is.na(Eggs_bottom), 
                                      round(CumulPropEggs * maxbot),
                                      Eggs_bottom),
         Eggs = Eggs_top + Eggs_bottom_imputed,
         NewEggs = c(NA, diff(Eggs)),
         NumAdultsRaw = StartCount - ifelse(is.na(Dead), 0, Dead),
         NumAdults = c(NumAdultsRaw[1], (NumAdultsRaw[2:length(NumAdultsRaw)] + NumAdultsRaw[1:(length(NumAdultsRaw)-1)])/2)) %>%
  mutate(NewEggs = ifelse(NewEggs < 0, 0, NewEggs)) %>% 
  mutate(OvipRate = NewEggs / 
           as.numeric(difftime(Date, lag(Date, 1), units = "days")) / 
           NumAdults,
         OvipRateDD = NewEggs / 
           c(NA, diff(AccumDD)) / 
           NumAdults,
         DDelapsed = c(NA, diff(AccumDD)),
         RateOffset = log(DDelapsed * NumAdults)) %>% 
  filter(is.na(OvipRate) == FALSE)


theme_set(theme_bw(base_size = 16)) 


# Figure A1 ----
# can change y to be based on day or degree-day
plt <- ggplot(dat, aes(x = AccumDD, y = OvipRateDD, group = Photoperiod, color = Photoperiod)) +
  geom_point(size = 1, alpha = .9) +
  geom_path(size = 2, alpha = .5) +
  scale_color_viridis() +
  facet_wrap(Population~Temperature, scales = "free_y") +
  xlab("Accumulated degree-days") +
  ylab("Oviposition rate per degree-day per adult") +
  ggtitle("Oviposition rate (new eggs / # deg-days / # adults)")
plt

# # alternative based on day
# plt <- ggplot(dat, aes(x = Photoperiod, y = OvipRate, color = AccumDD)) +
#   geom_point(size = 2.5, alpha = .5) +
#   scale_color_viridis() +
#   facet_wrap(Population~Temperature) +
#   theme_bw() +
#   ggtitle("Oviposition rate (new eggs / # days / # adults)")
# plt


ggsave("figures/figA1.png", plot = plt, device = "png", height = 9, width = 12, units = "in")




