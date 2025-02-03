### Dispersion Tests Project
## Melina Leite
# Feb 25

library(DHARMa)
library(tidyverse)
library(here)
library(cowplot);
theme_set(theme_cowplot())


#############----###
##### Binomial #####
##############--###

load(here("data", "7_glmmBin_alternative_power.Rdata")) # simulated data

#time spent + stupid way to convert mins to hours for some times
tspent <- map_dbl(out.binA, "time")
stime <- map(out.binA, "time")
#tspent[map(stime, attr, "units") == "mins"] <- tspent[map(stime, attr, "units") 
 #                                                     == "mins"]/60
mean(tspent) # 1

#simulations
simuls.bin <- map_dfr(out.binA, "simulations", .id="model") %>%
  separate(model, c("ngroups", "sampleSize", "intercept"), sep="_") %>%
  rename("overdispersion" = "controlValues")

# power
p.bin <- simuls.bin %>% 
  group_by(sampleSize, ngroups, intercept, overdispersion) %>%
  summarise(p.sig = sum(Alt.p<0.05,na.rm=T),
            nsim = length(Alt.p[!is.na(Alt.p)]))
p.bin$prop.sig <- p.bin$p.sig/p.bin$nsim
p.bin$intercept <- fct_relevel(p.bin$intercept, "-1.5", "0", "1.5")
p.bin$ngroups <- fct_relevel(p.bin$ngroups, "10", "50", "100")
p.bin$sampleSize <- as.factor(as.numeric(p.bin$sampleSize))

###### figure power ####

## all groups
p.bin %>% 
  ggplot(aes(x=overdispersion, y=prop.sig,linetype=ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Binomial", subtitle = "100 sim; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=2, byrow=TRUE))
#ggsave(here("figures", "7_glmmBin_Alternative_power.jpeg"), width=12, height = 15)

#############----###
##### Poisson #####
##############--###



load(here("data", "7_glmmPois_alternative_power.Rdata")) # simulated data

#time spent + stupid way to convert mins to hours for some times
tspent <- map_dbl(out.poisA, "time")
stime <- map(out.poisA, "time")
#tspent[map(stime, attr, "units") == "mins"] <- tspent[map(stime, attr, "units") 
#                                                     == "mins"]/60
mean(tspent) # 1

#simulations
simuls.pois <- map_dfr(out.poisA, "simulations", .id="model") %>%
  separate(model, c("ngroups", "sampleSize", "intercept"), sep="_") %>%
  rename("overdispersion" = "controlValues")

# power
p.pois <- simuls.pois %>% 
  group_by(sampleSize, ngroups, intercept, overdispersion) %>%
  summarise(p.sig = sum(Alt.p<0.05,na.rm=T),
            nsim = length(Alt.p[!is.na(Alt.p)]))
p.pois$prop.sig <- p.pois$p.sig/p.pois$nsim
p.pois$intercept <- fct_relevel(p.pois$intercept, "-1.5", "0", "1.5")
p.pois$ngroups <- fct_relevel(p.pois$ngroups, "10", "50", "100")
p.pois$sampleSize <- as.factor(as.numeric(p.pois$sampleSize))

###### figure power ####

## all groups
p.pois %>% 
  ggplot(aes(x=overdispersion, y=prop.sig,linetype=ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Poisson", subtitle = "100 sim") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=2, byrow=TRUE))

