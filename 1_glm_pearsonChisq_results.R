### Dispersion Tests Project
## Melina Leite
# Dec 24

library(DHARMa)
library(tidyverse)
library(here)
library(cowplot);
theme_set(theme_cowplot())
library(patchwork)

####################
##### Binomial #####
####################

load(here("data", "1_glmBin_pearsonChisq.Rdata"))


# manipulaing results
names(final.res.bin) <- 1:length(final.res.bin)

final.bin <- bind_rows(final.res.bin, .id="sim") %>% group_by(controlValues,intercept) %>%
  summarise(ks.sig = sum(ks.p<0.05))

for (i in 1:nrow(final.bin)) {
  confs <- binom.test(final.bin$ks.sig[i], 100)$conf.int
  final.bin$conf.low[i] <-  round(confs[1],2)
  final.bin$conf.up[i] <-  round(confs[2],2)
}

# shape distributions
sims.bin <- bind_rows(final.sims.bin, .id="sim")

# Mean distribution of 100 simulations
sims.mean.bin <- sims.bin %>% group_by(sim, controlValues,intercept) %>%
  arrange(sim, controlValues, intercept, Pear.stat) %>%
  mutate(n = 1:n()) %>%
  group_by(controlValues,intercept,n) %>%
  summarise(mean.pearson = mean(Pear.stat))



##### Figures #####

pbin <- ggplot(final.bin, aes(y=ks.sig/100, x=as.factor(controlValues),
                              col=as.factor(intercept))) +
  geom_point(position = position_dodge(width=0.4)) +
  geom_line(aes(x=as.numeric(as.factor(controlValues))),
            position = position_dodge(width=0.4))+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up),width = 0.1,
                position = position_dodge(width=0.4)) +
  geom_hline(yintercept = 0.05, linetype="dashed") +
  ylim(0,1)+
  scale_color_discrete("intercept")+
  xlab("sampleSize") + ylab("Prop of significant KS test") +
  labs(title="Binomial", tag = "B)") +
  theme(panel.border  = element_rect(color = "black"))
pbin
ggsave(here("figures", "1_glmBin_pearsonChisq.jpeg"), width = 6, height = 4)


# Figure all simulations

sims %>% filter(controlValues == 10) %>%
  ggplot(aes(x=Pear.stat, group=sim)) +
  geom_density(col=2, alpha=0.5) +
  facet_grid(controlValues~intercept)+
  stat_function(fun = dchisq, args = list(df = 8), col="black")+
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none") +
  
  sims %>% filter(controlValues == 50) %>%
  ggplot(aes(x=Pear.stat, group=sim)) +
  geom_density(col=2, alpha=0.5) +
  facet_grid(controlValues~intercept)+
  stat_function(fun = dchisq, args = list(df = 48), col="black")+
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  
  sims %>% filter(controlValues == 100) %>%
  ggplot(aes(x=Pear.stat, group=sim)) +
  geom_density(col=2, alpha=0.5) +
  facet_grid(controlValues~intercept)+
  stat_function(fun = dchisq, args = list(df = 98), col="black")+
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  
  sims %>% filter(controlValues == 500) %>%
  ggplot(aes(x=Pear.stat, group=sim)) +
  geom_density(col=2, alpha=0.5) +
  facet_grid(controlValues~intercept)+
  stat_function(fun = dchisq, args = list(df = 498), col="black")+
  xlab("Pearson Stats") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  plot_layout(ncol=1) + plot_annotation(title="Binomial: Pearson Stats distribution")
ggsave(here("figures", "1_glmBin_pearsonChisq_distrib.jpeg"), width=15,
       height=15)

# figure mean distribution

sims.mean.bin %>% filter(controlValues == 10) %>%
  ggplot(aes(x=mean.pearson)) +
  facet_grid(controlValues~intercept)+
  geom_density(col=2) +
  stat_function(fun = dchisq, args = list(df = 8), col="black") +
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none") +
  
  sims.mean.bin %>% filter(controlValues == 50) %>%
  ggplot(aes(x=mean.pearson)) +
  facet_grid(controlValues~intercept)+
  geom_density(col=2) +
  stat_function(fun = dchisq, args = list(df = 48), col="black") +
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  
  sims.mean.bin %>% filter(controlValues == 100) %>%
  ggplot(aes(x=mean.pearson)) +
  facet_grid(controlValues~intercept)+
  geom_density(col=2) +
  stat_function(fun = dchisq, args = list(df = 98), col="black") +
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  
  sims.mean.bin %>% filter(controlValues == 500) %>%
  ggplot(aes(x=mean.pearson)) +
  facet_grid(controlValues~intercept)+
  geom_density(col=2) +
  stat_function(fun = dchisq, args = list(df = 498), col="black") +
  xlab("Pearson Stats") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  
  plot_layout(ncol=1) + 
  plot_annotation(title="Binomial: Mean Pearson Stats distribution")

ggsave(here("figures", "1_glmBin_pearsonChisq_distrib_MEAN.jpeg"), width=15,
       height=15)



###################
##### Poisson #####
###################


load(here("data", "1_glmPois_pearsonChisq.Rdata"))

# manipulaing results
names(final.res.pois) <- 1:length(final.res.pois)

final.pois <- bind_rows(final.res.pois, .id="sim") %>% 
  group_by(controlValues,intercept) %>%
  summarise(ks.sig = sum(ks.p<0.05))

for (i in 1:nrow(final.pois)) {
  confs <- binom.test(final.pois$ks.sig[i], 100)$conf.int
  final.pois$conf.low[i] <-  round(confs[1],2)
  final.pois$conf.up[i] <-  round(confs[2],2)
}


# shape distributions
sims.pois <- bind_rows(final.sims.pois, .id="sim")

sims.mean.pois <- sims.pois %>% group_by(sim, controlValues,intercept) %>%
  arrange(sim, controlValues, intercept, Pear.stat) %>%
  mutate(n = 1:n()) %>%
  group_by(controlValues,intercept,n) %>%
  summarise(mean.pearson = mean(Pear.stat))




##### Figures #####

ggplot(final.pois, aes(y=ks.sig/100, x=as.factor(controlValues),
                       col=as.factor(intercept))) +
  geom_point(position = position_dodge(width=0.4)) +
  geom_line(aes(x=as.numeric(as.factor(controlValues))),
            position = position_dodge(width=0.4))+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up),width = 0.1,
                position = position_dodge(width=0.4)) +
  geom_hline(yintercept = 0.05, linetype="dashed") +
  scale_color_discrete("intercept")+
  xlab("sampleSize") + ylab("Prop of significant KS test") +
  labs(title="Poisson", tag = "A)") +
  theme(panel.border  = element_rect(color = "black")) +
  pbin + plot_layout(ncol=1, guides="collect")
ggsave(here("figures", "1_glmBOTH_pearsonChisq.jpeg"), width = 6, height = 8)



ggsave(here("figures", "1_glmPois_pearsonChisq.jpeg"), width = 6, height = 8)


# trick to include the binomial plot in the same figure / whe the other code is also loaded
# +
#  pbin + plot_layout(ncol=1, guides="collect")
#ggsave(here("figures", "1_glmBOTH_pearsonChisq.jpeg"), width = 6, height = 8)


# all simulations
sims.pois %>% filter(controlValues == 10) %>%
  ggplot(aes(x=Pear.stat, group=sim)) +
  geom_density(col=2, alpha=0.5) +
  facet_grid(controlValues~intercept)+
  stat_function(fun = dchisq, args = list(df = 8), col="black")+
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none") +
  
sims.pois %>% filter(controlValues == 50) %>%
  ggplot(aes(x=Pear.stat, group=sim)) +
  geom_density(col=2, alpha=0.5) +
  facet_grid(controlValues~intercept)+
  stat_function(fun = dchisq, args = list(df = 48), col="black")+
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  
sims.pois %>% filter(controlValues == 100) %>%
  ggplot(aes(x=Pear.stat, group=sim)) +
  geom_density(col=2, alpha=0.5) +
  facet_grid(controlValues~intercept)+
  stat_function(fun = dchisq, args = list(df = 98), col="black")+
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  
sims.pois %>% filter(controlValues == 500) %>%
  ggplot(aes(x=Pear.stat, group=sim)) +
  geom_density(col=2, alpha=0.5) +
  facet_grid(controlValues~intercept)+
  stat_function(fun = dchisq, args = list(df = 498), col="black")+
  xlab("Pearson Stats") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  
  plot_layout(ncol=1) + plot_annotation(title="Poisson: Pearson Stats distribution")
ggsave(here("figures", "1_glmPois_pearsonChisq_distrib.jpeg"), width=15,
       height=15)




sims.mean %>% filter(controlValues == 10) %>%
  ggplot(aes(x=mean.pearson, col=model)) +
  facet_grid(controlValues~intercept)+
  stat_density(geom="line",size=1.1, position="identity")+
  stat_function(fun = dchisq, args = list(df = 8), aes(color = "Chi-squared")) +
  scale_color_manual(values = c("binomial" ="olivedrab4", 
                                "poisson" = "mediumpurple3", 
                                "Chi-squared" = "black"))+
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = c(0.87,0.85),
        legend.key = element_rect()) +
  
  sims.mean %>% filter(controlValues == 50) %>%
  ggplot(aes(x=mean.pearson, col=model)) +
  facet_grid(controlValues~intercept)+
  geom_density(size=1.1) +
  scale_color_manual(values = c("olivedrab4", "mediumpurple3"))+
  stat_function(fun = dchisq, args = list(df = 48), col="black") +
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  
  
  sims.mean %>% filter(controlValues == 100) %>%
  ggplot(aes(x=mean.pearson, col=model)) +
  facet_grid(controlValues~intercept)+
  geom_density(size=1.1) +
  scale_color_manual(values = c("olivedrab4", "mediumpurple3"))+
  stat_function(fun = dchisq, args = list(df = 98), col="black") +
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  
  sims.mean %>% filter(controlValues == 500) %>%
  ggplot(aes(x=mean.pearson, col=model)) +
  facet_grid(controlValues~intercept)+
  geom_density(size=1.1) +
  scale_color_manual(values = c("olivedrab4", "mediumpurple3"))+
  stat_function(fun = dchisq, args = list(df = 498), col="black") +
  xlab("Pearson Stats") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  
  plot_layout(ncol=1) + 
  plot_annotation(title="Pearson Statistics X Chi-squared distribution")

ggsave(here("figures", "1_glm_pearsonChisq_distrib_MEAN.jpeg"), width=12,
       height=10)














