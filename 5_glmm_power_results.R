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

load(here("data", "5_glmmBin_power_10.Rdata")) # simulated data
out.bin10 <- out.bin
load(here("data", "5_glmmBin_power_50.Rdata")) # simulated data
out.bin50 <- out.bin
load(here("data", "5_glmmBin_power_100.Rdata")) # simulated data
out.bin100 <- out.bin

out.bin <- flatten(list(out.bin10, out.bin50, out.bin100))

#time spent hours
tspent <- map_dbl(out.bin, "time")


#simulations
simuls.bin <- map_dfr(out.bin, "simulations", .id="model") %>%
  separate(model, c("ngroups", "sampleSize", "intercept"), sep="_") %>%
  rename("overdispersion" = "controlValues")

# power
p.bin <- simuls.bin %>% dplyr::select(Pear.p.val, dhaUN.p.val, dhaCO.p.val,
                                      refCO.p.val, replicate,
                                      ngroups,
                                      overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:4, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize, ngroups, intercept, overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]))
p.bin$prop.sig <- p.bin$p.sig/p.bin$nsim
p.bin$intercept <- fct_relevel(p.bin$intercept, "-3", "-1.5", "0", "1.5", "3")
p.bin$ngroups <- fct_relevel(p.bin$ngroups, "10", "50", "100")
p.bin$sampleSize <- as.factor(as.numeric(p.bin$sampleSize))

###### figure power ####

## all groups
p.bin %>% filter(test != "Pear.p.val") %>%
  ggplot(aes(x=overdispersion, y=prop.sig, col=test, linetype=ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  # scale_color_discrete(
  #   labels=c("Sim-based conditional","Sim-based unconditional",  
  #            "Pearson ParBoot."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Binomial", subtitle = "1000 sim; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=2, byrow=TRUE))
#ggsave(here("figures", "5_glmmBin_powerALL.jpeg"), width=12, height = 15)

# 10 groups
p.bin %>% filter(ngroups == "10") %>%
ggplot(aes(x=overdispersion, y=prop.sig, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson ParBoot"))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Binomial", subtitle = "100 sim; 10 groups; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=4, byrow=TRUE))
#ggsave(here("figures", "5_glmmBin_power_10g.jpeg"), width=12, height = 15)

# 50 groups
p.bin %>% filter(ngroups == "50") %>%
  ggplot(aes(x=overdispersion, y=prop.sig, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson ParBoot."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Binomial", subtitle = "100 sim; 50 groups; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=4, byrow=TRUE))
#ggsave(here("figures", "5_glmmBin_power_50g.jpeg"), width=12, height = 9)


# 100 groups
p.bin %>% filter(ngroups == "100") %>%
  ggplot(aes(x=overdispersion, y=prop.sig, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson ParBoot."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Binomial", subtitle = "100 sim; 100 groups; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=4, byrow=TRUE))
#ggsave(here("figures", "5_glmmBin_power_100g.jpeg"), width=12, height = 9)




###### dispersion stat ####
d.bin <- simuls.bin %>% dplyr::select(Pear.stat.dispersion, dhaUN.stat.dispersion,
                                    dhaCO.stat.dispersion,
                                      refCO.stat.dispersion, replicate, ngroups,
                                      overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:4, names_to = "test", values_to = "dispersion") %>%
group_by(sampleSize,ngroups,intercept,overdispersion, test) %>%
  summarise(mean.stat = mean(dispersion, na.rm=T))
d.bin$intercept <- fct_relevel(d.bin$intercept, "-3", "-1.5", "0", "1.5", "3")
d.bin$ngroups <- fct_relevel(d.bin$ngroups, "10", "50", "100")
d.bin$sampleSize <- as.factor(as.numeric(d.bin$sampleSize))


# all groups
d.bin %>% filter(test != "Pear.stat.dispersion") %>%
  ggplot( aes(x=overdispersion, y=mean.stat, col=test, linetype=ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson Param. Bootstrap."))+
  facet_grid(sampleSize~intercept) +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0.75, ymax = 1,
           alpha = .1,fill = "red")+
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Binomial: dispersion statistics", subtitle = "100 sim; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=2, byrow=TRUE))
#ggsave(here("figures", "5_glmmBin_dispersionStatsALL.jpeg"), width=12, height = 15)


# 10 groups
d.bin %>% filter(ngroups == "10") %>%
ggplot( aes(x=overdispersion, y=mean.stat, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson Param. Bootstrap."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Binomial: dispersion statistics", subtitle = "100 sim; 10 groups; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=4, byrow=TRUE))
#ggsave(here("figures", "5_glmmBin_dispersionStats_10g.jpeg"), width=12, height = 15)

# 50 groups
d.bin %>% filter(ngroups == "50") %>%
  ggplot( aes(x=overdispersion, y=mean.stat, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson Param. Bootstrap."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Binomial: dispersion statistics", subtitle = "100 sim; 50 groups; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=4, byrow=TRUE))
#ggsave(here("figures", "5_glmmBin_dispersionStats_50g.jpeg"), width=12, height = 9)


# 100 groups
d.bin %>% filter(ngroups == "100") %>%
  ggplot( aes(x=overdispersion, y=mean.stat, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson Param. Bootstrap."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Binomial: dispersion statistics", subtitle = "100 sim; 100 groups; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=4, byrow=TRUE))
#ggsave(here("figures", "5_glmmBin_dispersionStats_100g.jpeg"), width=12, height = 9)


###### type 1 error #####

  
p.bin %>% filter(overdispersion == 0, test != "Pear.p.val") %>% ungroup() %>%
  mutate(ngroups = fct_relevel(ngroups, "10", "50", "100")) %>%
  ggplot(aes(x=sampleSize, y=prop.sig, col=intercept))+
  geom_point( position = position_dodge(width = 0.9))+
  geom_hline(yintercept = 0.05, linetype="dotted")+
  geom_line(aes(x=as.numeric(sampleSize)),
            position = position_dodge(width = 0.9))+
  facet_grid(ngroups~test, 
             labeller = as_labeller(c(`dhaCO.p.val` = "Conditional Sim-based" ,
                                      `dhaUN.p.val` = "Unconditional Sim-based" ,
                                      `refCO.p.val` = "Conditional Pear. param. boot.",
                                     `10` = " m = 10",`50` = " m = 50",
                                      `100` = " m = 100")))+
  theme(panel.background = element_rect(color="black"),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position = c(0.01,0.9)) 
ggsave(here("figures", "5_glmmBin_type1.jpeg"), width=12, height = 8)

#############----###
##### Poisson #####
##############--###

load(here("data", "5_glmmPois_power_10.Rdata")) # simulated data
out.pois10 <- out.pois
load(here("data", "5_glmmPois_power_50.Rdata")) # simulated data
out.pois50 <- out.pois
load(here("data", "5_glmmPois_power_100.Rdata")) # simulated data
out.pois100 <- out.pois

out.pois <- flatten(list(out.pois10, out.pois50, out.pois100))

#time spent + stupid way to convert mins to hours for some times
(tspent <- map_dbl(out.pois, "time"))


#simulations
simuls.pois <- map_dfr(out.pois, "simulations", .id="model") %>%
  separate(model, c("ngroups", "sampleSize", "intercept"), sep="_") %>%
  rename("overdispersion" = "controlValues")

# power
p.pois <- simuls.pois %>% dplyr::select(Pear.p.val, dhaUN.p.val, dhaCO.p.val,
                                      refCO.p.val, replicate,
                                      ngroups,
                                      overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:4, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize, ngroups, intercept, overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]))
p.pois$prop.sig <- p.pois$p.sig/p.pois$nsim
p.pois$intercept <- fct_relevel(p.pois$intercept, "-3", "-1.5", "0", "1.5", "3")
p.pois$ngroups <- fct_relevel(p.pois$ngroups, "10", "50", "100")
p.pois$sampleSize <- as.factor(as.numeric(p.pois$sampleSize))

###### figure power ####

p.pois %>% filter(test != "Pear.p.val") %>%
  ggplot(aes(x=overdispersion, y=prop.sig, col=test, linetype=ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson ParBoot. conditional"))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Poisson", subtitle = "100 sim") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=2, byrow=TRUE))
#ggsave(here("figures", "5_glmmPois_powerALL.jpeg"), width=12, height = 15)


# 10 groups
p.pois %>% filter(ngroups == "10") %>%
  ggplot(aes(x=overdispersion, y=prop.sig, col=test, linetype= ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson ParBoot."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Poisson", subtitle = "100 sim; 10 groups") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=4, byrow=TRUE))
#ggsave(here("figures", "5_glmmPois_power_10g.jpeg"), width=12, height = 12)

# 50 groups
p.pois %>% filter(ngroups == "50") %>%
  ggplot(aes(x=overdispersion, y=prop.sig, col=test, linetype= ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson ParBoot."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Poisson", subtitle = "100 sim; 50 groups0") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=4, byrow=TRUE))
#ggsave(here("figures", "5_glmmPois_power_50g.jpeg"), width=12, height = 9)


# 100 groups
p.pois %>% filter(ngroups == "100") %>%
ggplot(aes(x=overdispersion, y=prop.sig, col=test, linetype= ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson ParBoot."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Poisson", subtitle = "100 sim; 100 groups") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=4, byrow=TRUE))
#ggsave(here("figures", "5_glmmPois_power_100g.jpeg"), width=12, height = 9)



###### dispersion stat ####
d.pois <- simuls.pois %>% dplyr::select(Pear.stat.dispersion, dhaUN.stat.dispersion,
                                      dhaCO.stat.dispersion, 
                                      refCO.stat.dispersion, replicate, ngroups,
                                      overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:4, names_to = "test", values_to = "dispersion") %>%
  group_by(sampleSize, ngroups, intercept,overdispersion, test) %>%
  summarise(mean.stat = mean(dispersion, na.rm=T))
d.pois$intercept <- fct_relevel(d.pois$intercept, "-3", "-1.5", "0", "1.5", "3")
d.pois$ngroups <- fct_relevel(d.pois$ngroups, "10", "50", "100")
d.pois$sampleSize <- as.factor(as.numeric(d.pois$sampleSize))


# all groups
d.pois %>% filter(test != "Pear.stat.dispersion",
                  mean.stat <1000) %>%
  ggplot( aes(x=overdispersion, y=mean.stat, col=test, linetype=ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_y_log10()+
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Param. Bootstrap."))+
  facet_grid(sampleSize~intercept, scales="free") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0.5, ymax = 1,
           alpha = .1,fill = "red")+
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Poisson: dispersion statistics", subtitle = "100 sim") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=2, byrow=TRUE))
#ggsave(here("figures", "5_glmmPois_dispersionStatsALL.jpeg"), width=12, height = 15)


# 10 groups
d.pois %>% filter(ngroups == "10") %>%
  ggplot(aes(x=overdispersion, y=mean.stat, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson Param. Bootstrap."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Poisson: dispersion statistics", subtitle = "100 sim; 10 groups") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  scale_y_log10() + #ylim(0,3) +
  guides(color=guide_legend(nrow=4, byrow=TRUE))
#ggsave(here("figures", "5_glmmPois_dispersionStats_10g.jpeg"), width=12, height = 12)

# 50 groups
d.pois %>% filter(ngroups == "50") %>%
  filter(mean.stat <110) %>% # excluding weird results 
  ggplot(aes(x=overdispersion, y=mean.stat, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson Param. Bootstrap."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Poisson: dispersion statistics", subtitle = "100 sim; 50 groups") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  scale_y_log10() + #ylim(0,3) +
  guides(color=guide_legend(nrow=4, byrow=TRUE))
#ggsave(here("figures", "5_glmmPois_dispersionStats_50g.jpeg"), width=12, height = 9)


# 100 groups
d.pois %>% filter(ngroups == "100") %>%
ggplot(aes(x=overdispersion, y=mean.stat, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Chi-squared",
             "Pearson Param. Bootstrap."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Poisson: dispersion statistics", subtitle = "100 sim; 100 groups") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  scale_y_log10() + #ylim(0,3) +
  guides(color=guide_legend(nrow=4, byrow=TRUE))
#ggsave(here("figures", "5_glmmPois_dispersionStats_100g.jpeg"), width=12, height = 9)

###### Type 1 error #####

p.pois %>% filter(overdispersion == 0, test != "Pear.p.val") %>% ungroup() %>%
  mutate(ngroups = fct_relevel(ngroups, "10", "50", "100")) %>%
  ggplot(aes(x=sampleSize, y=prop.sig, col=intercept))+
  geom_point( position = position_dodge(width = 0.9))+
  geom_hline(yintercept = 0.05, linetype="dotted")+
  geom_line(aes(x=as.numeric(sampleSize)),
            position = position_dodge(width = 0.9))+
  facet_grid(ngroups~test, 
             labeller = as_labeller(c(`dhaCO.p.val` = "Conditional Sim-based" ,
                                      `dhaUN.p.val` = "Unconditional Sim-based" ,
                                      `refCO.p.val` = "Conditional Pear. param. boot.",
                                      `10` = " m = 10",`50` = " m = 50",
                                      `100` = " m = 100")))+
  theme(panel.background = element_rect(color="black"),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position = c(0.01,0.9)) 
#ggsave(here("figures", "5_glmmPois_type1.jpeg"), width=12, height = 8)



#### FIGURES ALL ####

# all simulations
pow <- bind_rows(list(Poisson = p.pois, Binomial = p.bin), .id="model") %>%
  ungroup() %>%
  mutate(model= fct_relevel(model, "Poisson", "Binomial"),
         ngroups = fct_relevel(ngroups, "10", "50", "100"))
# add sig test
for (i in 1:nrow(pow)) {
  btest <- binom.test(pow$p.sig[i], n=pow$nsim[i], p=0.05)
  pow$p.bin0.05[i] <- btest$p.value
  pow$conf.low[i] <- btest$conf.int[1]
  pow$conf.up[i] <- btest$conf.int[2]
}



disp <- bind_rows(list(Poisson = d.pois, Binomial = d.bin), .id="model") %>%
  ungroup() %>%
  mutate(model= fct_relevel(model, "Poisson", "Binomial"),
         ngroups = fct_relevel(ngroups, "10", "50", "100"))


###### type 1 error ####

# intercept = 0, sample 1000
type1 <- pow %>% filter(test != "Pear.p.val", overdispersion == 0, 
                        intercept==0) %>%
  ggplot(aes(x=sampleSize, y=prop.sig, col=ngroups))+
  geom_point(position=position_dodge(width=0.8)) + 
  geom_line(aes(x=as.numeric(as.factor(sampleSize))), position=position_dodge(width=0.8))+
  ylab("Type I error")+ xlab("")+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up), width=0.1,
                position=position_dodge(width=0.8))+
  facet_grid(model~test, labeller = as_labeller(c(`dhaCO.p.val` = "Conditional Sim-based" ,
                                                    `dhaUN.p.val` = "Unconditional Sim-based" ,
                                                    `refCO.p.val` = "Pear. param. boot.",
                                                    `Binomial` = "Binomial",
                                                    `Poisson` = "Poisson"))) +
  geom_hline(yintercept = 0.05, linetype="dashed")+
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom",
        axis.text.x = element_text(angle=45, hjust=1)) 
type1
ggsave(here("figures", "5_glmm_type1.jpeg"), width = 12, heigh=6)



###### Power ####

# intercept = 0, sample 1000
fig.pow <- pow %>% filter(intercept == 0, 
                          !test %in% c("Pear.p.val", "refUN.p.val" ),
               sampleSize %in% c(1000)) %>%
  ggplot(aes(x=overdispersion, y= prop.sig, col=test, linetype=sampleSize)) +
  geom_point() + geom_line() +
  #geom_errorbar(aes(ymin=conf.low, ymax=conf.up), width=0.01)+
  facet_grid(model~ngroups, labeller = as_labeller(c(`10`= "m = 10 groups",
                                                     `50`= "m = 50 groups",
                                                     `100`= "m = 100 groups",
                                                     `Binomial` = "Binomial",
                                                     `Poisson` = "Poisson"))) +
  annotate("rect", xmin = -0.05, xmax = 0.05, ymin = 0, ymax = 1,
           alpha = .1,fill = "blue")+
  ylab("Power")+
  scale_color_manual(values = col.tests[4:2],
    labels=c("Sim-based conditional","Sim-based unconditional",
             "Pearson Param. Boostrap")) +
  theme(panel.background = element_rect(color="black"),
        legend.box.background = element_rect(fill = "gray94", color="gray94"),
        legend.position = "none") +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  labs(tag="A)") +
  geom_hline(yintercept = 0.05, linetype="dashed")+
  geom_hline(yintercept = 0.5, linetype="dotted")
fig.pow

###### Dispersion stat ####

# intercept = 0, sample 1000
fig.disp <- disp %>% filter(intercept == 0, 
                            !test %in% c("Pear.stat.dispersion"),
               sampleSize == 1000) %>%
  ggplot(aes(x=overdispersion, y= mean.stat, col=test, )) +
  geom_point() + geom_line() +
  facet_grid(model~ngroups, labeller = as_labeller(c(`10`= "m = 10 groups",
                                                     `50`= "m = 50 groups",
                                                     `100`= "m = 100 groups",
                                                     `Binomial` = "Binomial",
                                                     `Poisson` = "Poisson"))) +
  ylab("Dispersion statistics")+
  geom_hline(yintercept = 1, linetype="dotted") +
  scale_color_manual(values = col.tests[4:2],
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Param. Bootstrap.")) +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom",
        legend.box.background = element_rect(fill = "gray94", color="gray94"))+
  guides(color=guide_legend(nrow=1, byrow=TRUE)) +
  labs(tag="B)") +
  scale_y_log10()+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0.75, ymax = 1,
           alpha = .1,fill = "red")
fig.disp

###### Combined power & dispersion ####

# intercept = 0, sample 1000
fig.pow + fig.disp + plot_layout(ncol=1)+
  plot_annotation(title="Alternative dispersion test for GLMMs",
                  theme = theme(plot.title = element_text(hjust=0.5)))
ggsave(here("figures", "5_glmm_results.jpeg"), width = 10, heigh=12)















## ALTERNATIVE ####
## 

library(ggh4x)

pow %>% filter(intercept == 0, test != "Pear.p.val",
               sampleSize %in% c(200,500,1000)) %>%
  ggplot(aes(x=overdispersion, y= prop.sig, col=test)) +
  geom_point(alpha=0.6) + geom_line(alpha=0.6) +
  facet_nested(sampleSize~model + ngroups) +
  ylab("Power") +  
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Param. Bootstrap.")) +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") +
  guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  ggtitle("GLMM Power")


pow %>% filter(intercept == 0, test != "Pear.p.val",
               overdispersion >= 0.2 & overdispersion <=0.5,
               sampleSize %in% c(200,500,1000)) %>%
  ggplot(aes(x=as.numeric(as.character(ngroups)), y= prop.sig, col=test)) +
  geom_point(alpha=0.6) + geom_line(alpha=0.6) +
  facet_nested(overdispersion ~ model + sampleSize) +
  scale_x_continuous(breaks=c(10,50,100))+
  ylab("Power") +  xlab("number of groups (m)")+
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Param. Bootstrap.")) +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") +
  guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  ggtitle("GLMM Power ")




pow %>% filter(intercept == 0, test != "Pear.p.val",
               overdispersion %in% c(0.2,0.5),
               sampleSize %in% c(200,500,1000)) %>%
  ggplot(aes(x=as.numeric(as.character(ngroups)), y= prop.sig, col=test, linetype=as.factor(overdispersion))) +
  geom_point(alpha=0.6) + geom_line(alpha=0.6) +
  facet_grid(model~sampleSize) +
  scale_x_continuous(breaks=c(10,50,100))+
  ylab("Power") +  xlab("number of groups (m)")+
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Param. Bootstrap.")) +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") +
  guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  ggtitle("GLMM Power for 0.5 of overdispersion")

pow %>% filter(intercept == 0, test != "Pear.p.val",
               overdispersion <=0.5,
               sampleSize %in% c(200,500,1000)) %>%
  ggplot(aes(x=overdispersion, y= prop.sig, col=test, linetype=ngroups)) +
  geom_point() + geom_line() +
  facet_grid(model~sampleSize) +
  ylab("Power") +
  scale_color_discrete(
    labels=c("Sim-based conditional","Sim-based unconditional",  
             "Pearson Param. Bootstrap.")) +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("GLMM Power")
