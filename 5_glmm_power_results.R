### Dispersion Tests Project
## Melina Leite
# Dez 24

library(DHARMa)
library(tidyverse)
library(here)
library(cowplot);
theme_set(theme_cowplot())
library(patchwork)


#############----###
##### Binomial #####
##############--###

load(here("data", "5_glmmBin_power.Rdata")) # simulated data

#time spent + stupid way to convert mins to hours for some times
tspent <- map_dbl(out.bin, "time")
stime <- map(out.bin, "time")
tspent[map(stime, attr, "units") == "mins"] <- tspent[map(stime, attr, "units") 
                                                      == "mins"]/60
mean(tspent) # 1

#simulations
simuls.bin <- map_dfr(out.bin, "simulations", .id="model") %>%
  separate(model, c("ngroups", "sampleSize", "intercept"), sep="_") %>%
  rename("overdispersion" = "controlValues")

# power
p.bin <- simuls.bin %>% dplyr::select(Pear.p.val, dhaUN.p.val, dhaCO.p.val,
                                      refUN.p.val, refCO.p.val, replicate,
                                      ngroups,
                                      overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:5, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize, ngroups, intercept, overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]))
p.bin$prop.sig <- p.bin$p.sig/p.bin$nsim
p.bin$intercept <- fct_relevel(p.bin$intercept, "-3", "-1.5", "0", "1.5", "3")
p.bin$sampleSize <- as.factor(as.numeric(p.bin$sampleSize))

###### figure power ####

# 10 groups
p.bin %>% filter(ngroups == "10") %>%
ggplot(aes(x=overdispersion, y=prop.sig, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson ParBoot. conditional",
             "Pearson ParBoot. unconditional"))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Binomial", subtitle = "100 sim; 10 groups; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=4, byrow=TRUE))
ggsave(here("figures", "5_glmmBin_power_10g.jpeg"), width=12, height = 15)


# 100 groups
p.bin %>% filter(ngroups == "100") %>%
  ggplot(aes(x=overdispersion, y=prop.sig, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson ParBoot. conditional",
             "Pearson ParBoot. unconditional"))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Binomial", subtitle = "100 sim; 10 groups; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=4, byrow=TRUE))
ggsave(here("figures", "5_glmmBin_power_100g.jpeg"), width=12, height = 10)




###### dispersion stat ####
d.bin <- simuls.bin %>% dplyr::select(Pear.stat.dispersion, dhaUN.stat.dispersion,
                                    dhaCO.stat.dispersion, refUN.stat.dispersion, 
                                      refCO.stat.dispersion, replicate, ngroups,
                                      overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:5, names_to = "test", values_to = "dispersion") %>%
group_by(sampleSize,ngroups,intercept,overdispersion, test) %>%
  summarise(mean.stat = mean(dispersion, na.rm=T))
d.bin$intercept <- fct_relevel(d.bin$intercept, "-3", "-1.5", "0", "1.5", "3")
d.bin$sampleSize <- as.factor(as.numeric(d.bin$sampleSize))


d.bin %>% filter(ngroups == "10") %>%
ggplot( aes(x=overdispersion, y=mean.stat, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson Param. Bootstrap. conditional",
             "Pearson Param. Boostrap. unconditional"))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Binomial: dispersion statistics", subtitle = "100 sim; 10 groups; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=4, byrow=TRUE))
ggsave(here("figures", "5_glmmBin_dispersionStats_10g.jpeg"), width=12, height = 15)

d.bin %>% filter(ngroups == "100") %>%
  ggplot( aes(x=overdispersion, y=mean.stat, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson Param. Bootstrap. conditional",
             "Pearson Param. Boostrap. unconditional"))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Binomial: dispersion statistics", subtitle = "100 sim; 100 groups; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=4, byrow=TRUE))
ggsave(here("figures", "5_glmmBin_dispersionStats_100g.jpeg"), width=12, height = 10)




#############----###
##### Poisson #####
##############--###

load(here("data", "5_glmmPois_power.Rdata")) # simulated data

#time spent + stupid way to convert mins to hours for some times
tspent <- map_dbl(out.pois, "time")
stime <- map(out.pois, "time")
tspent[map(stime, attr, "units") == "mins"] <- tspent[map(stime, attr, "units") 
                                                      == "mins"]/60
mean(tspent) # 1

#simulations
simuls.pois <- map_dfr(out.pois, "simulations", .id="model") %>%
  separate(model, c("ngroups", "sampleSize", "intercept"), sep="_") %>%
  rename("overdispersion" = "controlValues")

# power
p.pois <- simuls.pois %>% dplyr::select(Pear.p.val, dhaUN.p.val, dhaCO.p.val,
                                      refUN.p.val, refCO.p.val, replicate,
                                      ngroups,
                                      overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:5, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize, ngroups, intercept, overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]))
p.pois$prop.sig <- p.pois$p.sig/p.pois$nsim
p.pois$intercept <- fct_relevel(p.pois$intercept, "-3", "-1.5", "0", "1.5", "3")
p.pois$sampleSize <- as.factor(as.numeric(p.pois$sampleSize))

###### figure power ####

ggplot(p.pois, aes(x=overdispersion, y=prop.sig, col=test, linetype= ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson ParBoot. conditional",
             "Pearson ParBoot. unconditional"))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Poisson", subtitle = "100 sim; 100 groups; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=4, byrow=TRUE))
ggsave(here("figures", "5_glmmPois_power.jpeg"), width=12, height = 10)



###### dispersion stat ####
d.pois <- simuls.pois %>% dplyr::select(Pear.stat.dispersion, dhaUN.stat.dispersion,
                                      dhaCO.stat.dispersion, refUN.stat.dispersion, 
                                      refCO.stat.dispersion, replicate, ngroups,
                                      overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:5, names_to = "test", values_to = "dispersion") %>%
  group_by(sampleSize,intercept,overdispersion, test) %>%
  summarise(mean.stat = mean(dispersion, na.rm=T))
d.pois$intercept <- fct_relevel(d.pois$intercept, "-3", "-1.5", "0", "1.5", "3")
d.pois$sampleSize <- as.factor(as.numeric(d.pois$sampleSize))


ggplot(d.pois, aes(x=overdispersion, y=mean.stat, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson Param. Bootstrap. conditional",
             "Pearson Param. Boostrap. unconditional"))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Poisson: dispersion statistics", subtitle = "100 sim; 100 groups; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  scale_y_log10() + #ylim(0,3) +
  guides(color=guide_legend(nrow=4, byrow=TRUE))
ggsave(here("figures", "5_glmmPois_dispersionStats.jpeg"), width=12, height = 10)



#### FIGURE ALL ####

pow <- bind_rows(list(Poisson = p.pois, Binomial = p.bin), .id="model") %>%
  ungroup() %>%
  mutate(model= fct_relevel(model, "Poisson", "Binomial"))
disp <- bind_rows(list(Poisson = d.pois, Binomial = d.bin), .id="model") %>%
  ungroup() %>%
  mutate(model= fct_relevel(model, "Poisson", "Binomial"))


fig.pow <- pow %>% filter(intercept == 0, test != "Pear.p.val",
               sampleSize == 1000) %>%
  ggplot(aes(x=overdispersion, y= prop.sig, col=test)) +
  geom_point() + geom_line() +
  facet_grid(~model) +
  annotate("rect", xmin = -0.05, xmax = 0.05, ymin = 0, ymax = 1,
           alpha = .1,fill = "blue")+
  ylab("Power")+
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson Param. Bootstrap. conditional",
             "Pearson Param. Boostrap. unconditional")) +
  theme(panel.background = element_rect(color="black"),
        legend.position = "none")+
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  labs(tag="A)")

fig.disp <- disp %>% filter(intercept == 0, test != "Pear.stat.dispersion",
               sampleSize == 1000) %>%
  ggplot(aes(x=overdispersion, y= mean.stat, col=test)) +
  geom_point() + geom_line() +
  facet_grid(~model) +
  ylab("Dispersion statistics")+
  geom_hline(yintercept = 1, linetype="dotted")+
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson Param. Bootstrap. conditional",
             "Pearson Param. Boostrap. unconditional")) +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom")+
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  labs(tag="B)")


fig.pow + fig.disp + plot_layout(ncol=1)
ggsave(here("figures", "5_glmm_power_dispersion.jpeg"), width = 8, heigh=8)
