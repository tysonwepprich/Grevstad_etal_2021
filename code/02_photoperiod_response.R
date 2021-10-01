## Header -----
## Script name: 02_photoperiod_response.R
## Purpose of script: Estimate short-day diapause response
## Author: Tyson Wepprich
## Date Created: 2021-02-04
## License: GNU GPLv3
## Email: tyson.wepprich@gmail.com
## ---
## Notes:
## This script makes Figure 3 and 
## estimates photoperiod response (Table 1) from experimental data
## ---

library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)
library(MASS)
dat <- read.csv("data/photo_diapause.csv") %>% 
  mutate(Total = Repro + Diap,
         perc_repro = Repro/Total,
         perc_diap = Diap/Total)

# Probit model ----
# best model has additive Pop by AIC, not interaction or omitted
mod <- glm(cbind(Diap, Repro) ~ Treat + Pop - 1, data = dat, family = binomial(link = "probit"))
ilink <- family(mod)$linkinv

# coefs <- coef(modS)
pred <- data.frame(Treat = rep(seq(12.5, 16.5, length.out = 100), 2), Pop = c(rep("South", 100), rep("North", 100)))
preds <- predict(mod, newdata = pred, type = "link", se.fit = TRUE)      
pred$pred <- ilink(preds$fit)
pred$lwr <- ilink(preds$fit - 1.96*preds$se.fit)
pred$upr <- ilink(preds$fit + 1.96*preds$se.fit)

# Use LD50 function to estimate critical photoperiod and SE from GLM
# Results in Table 1 parameters
modcps <- bind_rows(lapply(c(3, 2), FUN = function(x) {
  tmp <- MASS::dose.p(mod, cf = c(x, 1), p = c(.05, .5, .95))
  data.frame(pred = as.vector(tmp),
             se = as.vector(attr(tmp, "SE")),
             p = as.vector(attr(tmp, "p")))
})) %>% 
  mutate(Pop = c(rep("S", 3), rep("N", 3))) %>%
  pivot_wider(names_from = p, values_from = c(pred, se)) %>% 
  mutate(cp_mean = -coef(mod)[3:2]/coef(mod)[1],
         cp_sd = 1/abs(coef(mod)[1]))

# Figure 3: Plot photoperiod response ----
plt <- ggplot(data = dat, aes(x = Treat, y = perc_repro, group = Pop)) +
  geom_point(data = dat, aes(shape = Pop), size = 4, alpha = .5) +
  scale_shape_discrete(name = NULL) +
  geom_line(data = pred, aes(x = Treat, y = 1 - pred, group = Pop, linetype = Pop)) +
  scale_linetype_manual(name = NULL, values = c("solid", "longdash")) +
  geom_ribbon(data = pred, aes(x = Treat, ymin = 1-lwr, ymax = 1-upr, group = Pop), alpha = .2, inherit.aes = FALSE) +
  xlab("Hours of daylight") +
  ylab("Percent reproductive") + 
  theme_bw(base_size = 16) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = c(.85, .15))
plt

ggsave("figures/fig2.png", dpi = 600,
       plot = plt, device = "png", width = 6, height = 5.5, units = "in")
