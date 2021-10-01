## Header ---------------------------
## Script name: 01_development_rates.R
## Purpose of script: Estimate physiological rates from growth chamber experiments
## Author: Tyson Wepprich
## Date Created: 2021-02-04
## License: GNU GPLv3
## Email: tyson.wepprich@gmail.com
## ---
## Notes: 
## This script makes Figure 2 from experimental development time data
## Linear regression used to estimate Table 1 parameters:
## lower developmental threshold and degree-days required from egg to adult
## ---


library(dplyr)
library(ggplot2)
library(Hmisc)

# Experimental data ----
# Myint et al. (2012) data
# copied from Table 1
myint_dat <- read.csv("data/dev_myint.csv", header = TRUE)

# group male and females together (combine mean/sd of two groups)
# only analyze full development from egg to emerging adult
myint_dev <- myint_dat %>% 
  group_by(sex, temperature) %>% 
  summarise(n = n, mean = total_mean, sd = total_se * sqrt(n)) %>% 
  ungroup() %>% 
  group_by(temperature) %>% 
  mutate(tn = n[1] + n[2],
         tmean = (n[1]*mean[1] + n[2]*mean[2]) / (n[1] + n[2]),
         tsd = sqrt(((n[1]-1)*sd[1]^2 + (n[2]-1)*sd[2]^2 + n[1] * n[2] / (n[1] + n[2]) * (mean[1]^2 + mean[2]^2 - 2 * mean[1] * mean[2])) / (n[1] + n[2] -1)),
         tse = tsd/sqrt(tn))

mdev <- myint_dev[1:4,-1]
myint_dev1 <- mdev %>% 
  mutate(n = tn,
         tlo = tmean - 1.96*tse,
         thi = tmean + 1.96*tse,
         rmean = 1/tmean,
         rlo = 1/tlo,
         rhi = 1/thi,
         experiment = "Myint et al. (2012)") %>% 
  dplyr::select(-mean, -sd, -tn, -tsd, -tse)
myint_dev1$experiment[4] <- "Myint et al. (excluded)"

# Experiment 1

shaw_dev <- read.csv("data/dev_exp1.csv", header = TRUE)
shaw_dev1 <- shaw_dev %>% 
  group_by(temperature) %>% 
  summarise(n = n(),
            tmean = smean.cl.normal(time1)[1],
            tlo = smean.cl.normal(time1)[2],
            thi = smean.cl.normal(time1)[3]) %>% 
  mutate(rmean = 1/tmean,
         rlo = 1/tlo,
         rhi = 1/thi,
         experiment = "Experiment 1")

# Experiment 2
trans_dev <- read.csv("data/dev_exp2.csv", header = TRUE)

# Get mean/CI from duration before converting to 1/duration for growth rate
# Average growth rate at 15C (temperature used to complete development) 
# is used as baseline for development rate observed at lower temperatures
mean15 <- trans_dev %>% 
  filter(temperature == 15) %>% 
  mutate(time1 = time1 + time2,
         devRate = 1 / time1)
mean15 <- mean(mean15$devRate)

trans_dev1 <- trans_dev %>% 
  group_by(temperature) %>% 
  summarise(n = n(),
            tmean = smean.cl.normal(time2)[1],
            tlo = smean.cl.normal(time2)[2],
            thi = smean.cl.normal(time2)[3]) %>% 
  mutate(rmean = (1 - tmean * mean15)/30,
         rlo = (1 - tlo * mean15)/30,
         rhi = (1 - thi * mean15)/30,
         experiment = "Experiment 2")

# Plot experimental 
alldat <- bind_rows(trans_dev1, shaw_dev1, myint_dev1) %>%
  # alldat <- bind_rows(trans2, shaw2, myint_dev1) %>% 
  mutate(ptshp = case_when(experiment == "Myint et al. (excluded)" ~ "a",
                           # temperature == 12 & experiment == "shaw" ~ "b",
                           experiment == "Myint et al. (2012)" ~ "c",
                           experiment == "Experiment 1" ~ "d",
                           experiment == "Experiment 2" ~ "e")) %>% 
  filter(!(temperature == 12 & experiment == "Experiment 1")) # excluded for poor survival/growth rate outlier
alldat$experiment <- factor(alldat$experiment, levels = c("Experiment 1", "Experiment 2", "Myint et al. (2012)", "Myint et al. (excluded)"))

theme_set(theme_bw(base_size = 20) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

# Plot development data ----
# Plot of mean development rate by experiment and treatment
meanplt <- ggplot(alldat, aes(x = temperature, y = rmean, shape = experiment)) +
  geom_point(size = 2.5) +
  geom_linerange(aes(x = temperature, ymin = rlo, ymax = rhi), size = 1) +
  scale_shape_manual(name = NULL, values = c(17, 15, 16, 1))
meanplt

# Plot of mean development time by experiment and treatment
dayplt <- ggplot(alldat, aes(x = temperature, y = tmean, shape = experiment)) +
  geom_point(size = 3) +
  geom_linerange(aes(x = temperature, ymin = tlo, ymax = thi), size = 1) +
  scale_shape_manual(name = NULL, values = c(17, 15, 16, 1))
dayplt


# Linear model ----
# Exclude outliers to estimate linear region of TPC, shaw 12C and myint 30C
lindat <- alldat %>% 
  filter(temperature < 29) %>% 
  filter(!(temperature == 12 & experiment == "Experiment 1"))
wmod <- lm(rmean ~ temperature, data = lindat, weights = n)
summary(wmod)
# Lower developmental threshold and total degree-days required (from wmod coefficients)
ldt <- -(-0.0126)/0.00183
tot_dd <- 1/0.00183

# Fig 2: Plot regression ----
newdat <- data.frame(temperature = seq(0, 30, .01))
newdat <- cbind(newdat, data.frame(predict.lm(wmod, newdat, se.fit = TRUE)))
newdat <- newdat %>% 
  filter(fit >= 0) %>% 
  mutate(rmean = fit,
         rhi = fit + 1.96*se.fit,
         rlo = fit - 1.96*se.fit)

fig2 <- meanplt +
  # geom_line(data = newdat, aes(x = temperature, y = devRate, group = experiment), inherit.aes = FALSE) +
  geom_path(data = newdat, aes(x = temperature, y = rmean), inherit.aes = FALSE, color = "black", size = 1.5, alpha = .5) +
  scale_x_continuous(limits = c(5, 30), breaks = c(5, 10, 15, 20, 25, 30)) +
  geom_ribbon(data = newdat, aes(x = temperature, ymin = rlo, ymax = rhi), inherit.aes = FALSE, alpha = .1) +
  xlab("Constant temperature (°C)") +
  ylab("Development rate (1/days)") +
  theme(legend.position = c(.23, .85))
# ggtitle("Aphalara daily development rate vs constant temperatures")
fig2
ggsave(filename = "figures/figS1.png", plot = fig2, device = "png", width = 8, height = 7)
