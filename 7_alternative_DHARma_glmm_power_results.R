### Dispersion Tests Project
## 
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
summary(tspent) # secs

#simulations
simuls.binA <- map_dfr(out.binA, "simulations", .id="model") %>%
  separate(model, c("ngroups", "sampleSize", "intercept"), sep="_") %>%
  rename("overdispersion" = "controlValues")

# power
p.binA <- simuls.binA %>% 
  group_by(sampleSize, ngroups, intercept, overdispersion) %>%
  summarise(p.sig = sum(Alt.p<0.05,na.rm=T),
            nsim = length(Alt.p[!is.na(Alt.p)]))
p.binA$prop.sig <- p.binA$p.sig/p.binA$nsim
p.binA$intercept <- fct_relevel(p.binA$intercept, "-1.5", "0", "1.5")
p.binA$ngroups <- fct_relevel(p.binA$ngroups, "10", "50", "100")
p.binA$sampleSize <- as.factor(as.numeric(p.binA$sampleSize))

###### figure power ####

p.binA %>% 
  ggplot(aes(x=overdispersion, y=prop.sig,linetype=ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Binomial", subtitle = "100 sim; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=2, byrow=TRUE))
#ggsave(here("figures", "7_glmmBin_Alternative_power.jpeg"), width=12, height = 15)




###### Combining with other tests ####

load(here("data", "5_glmmBin_power.Rdata")) # simulated data
#simulations
simuls.bin <- map_dfr(out.bin, "simulations", .id="model") %>%
  separate(model, c("ngroups", "sampleSize", "intercept"), sep="_") %>%
  rename("overdispersion" = "controlValues")

# power
p.bin <- simuls.bin %>% dplyr::select(Pear.p.val, dhaUN.p.val, dhaCO.p.val,
                                       refCO.p.val, replicate,
                                      ngroups,
                                      overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:5, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize, ngroups, intercept, overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]))
p.bin$prop.sig <- p.bin$p.sig/p.bin$nsim
p.bin$intercept <- fct_relevel(p.bin$intercept, "-3", "-1.5", "0", "1.5", "3")
p.bin$ngroups <- fct_relevel(p.bin$ngroups, "10", "50", "100")
p.bin$sampleSize <- as.factor(as.numeric(p.bin$sampleSize))


bindata <- bind_rows(list(p.binA %>% mutate(test = "approxPear"),
                     p.bin)) 

bindata %>% filter(ngroups==100, intercept %in% c("-1.5", "0", "1.5"),
                   !test %in% c("refCO.p.val", "Pear.p.val"),  
                   sampleSize != 10000) %>%
  ggplot(aes(x=overdispersion, y=prop.sig, col=test, linetype=ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  # scale_color_discrete(
  #   labels=c("Sim-based conditional","Sim-based unconditional",  
  #            "Pearson Chi-squared",
  #            "Pearson ParBoot. conditional",
  #            "Pearson ParBoot. unconditional"))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Binomial", subtitle = "100 sim; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=2, byrow=TRUE))

#ggsave(here("figures", "5_glmmBin_powerALL.jpeg"), width=12, height = 15)

#############----###
##### Poisson #####
##############--###



load(here("data", "7_glmmPois_alternative_power.Rdata")) # simulated data

#time spent + stupid way to convert mins to hours for some times
tspent <- map_dbl(out.poisA, "time")
stime <- map(out.poisA, "time")
#tspent[map(stime, attr, "units") == "mins"] <- tspent[map(stime, attr, "units") 
#                                                     == "mins"]/60
summary(tspent) # 1

#simulations
simuls.poisA <- map_dfr(out.poisA, "simulations", .id="model") %>%
  separate(model, c("ngroups", "sampleSize", "intercept"), sep="_") %>%
  rename("overdispersion" = "controlValues")

# power
p.poisA <- simuls.poisA %>% 
  group_by(sampleSize, ngroups, intercept, overdispersion) %>%
  summarise(p.sig = sum(Alt.p<0.05,na.rm=T),
            nsim = length(Alt.p[!is.na(Alt.p)]))
p.poisA$prop.sig <- p.poisA$p.sig/p.poisA$nsim
p.poisA$intercept <- fct_relevel(p.poisA$intercept, "-1.5", "0", "1.5")
p.poisA$ngroups <- fct_relevel(p.poisA$ngroups, "10", "50", "100")
p.poisA$sampleSize <- as.factor(as.numeric(p.poisA$sampleSize))

###### figure power ####

## all groups
p.poisA %>% 
  ggplot(aes(x=overdispersion, y=prop.sig,linetype=ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Poisson", subtitle = "100 sim") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=2, byrow=TRUE))


###### Combining with other tests ####

load(here("data", "5_glmmPois_power.Rdata")) # simulated data
#simulations
simuls.pois <- map_dfr(out.pois, "simulations", .id="model") %>%
  separate(model, c("ngroups", "sampleSize", "intercept"), sep="_") %>%
  rename("overdispersion" = "controlValues")

# power
p.pois <- simuls.pois %>% dplyr::select(Pear.p.val, dhaUN.p.val, dhaCO.p.val,
                                        refCO.p.val, replicate,
                                        ngroups,
                                        overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:5, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize, ngroups, intercept, overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]))
p.pois$prop.sig <- p.pois$p.sig/p.pois$nsim
p.pois$intercept <- fct_relevel(p.pois$intercept, "-3", "-1.5", "0", "1.5", "3")
p.pois$ngroups <- fct_relevel(p.pois$ngroups, "10", "50", "100")
p.pois$sampleSize <- as.factor(as.numeric(p.pois$sampleSize))


poisdata <- bind_rows(list(p.poisA %>% mutate(test = "approxPear"),
                          p.pois)) 

poisdata %>% filter(ngroups==100, intercept %in% c("-1.5", "0", "1.5"),
                    !test %in% c("refCO.p.val", "Pear.p.val"), 
                   sampleSize != 10000) %>%
  ggplot(aes(x=overdispersion, y=prop.sig, col=test, linetype=ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  # scale_color_discrete(
  #   labels=c("Sim-based conditional","Sim-based unconditional",  
  #            "Pearson Chi-squared",
  #            "Pearson ParBoot. conditional",
  #            "Pearson ParBoot. unconditional"))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Poisson", subtitle = "100 sim;") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=2, byrow=TRUE))

