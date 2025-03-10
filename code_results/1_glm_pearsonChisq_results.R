### Dispersion Tests Project
## Melina Leite
# Dec 24

library(DHARMa)
library(tidyverse)
library(here)
library(cowplot);
theme_set(theme_cowplot())
library(patchwork)

source(here("code_results", "plotColors.R"))

##############-----#
##### Binomial #####
###############----#

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

#################-#
##### Poisson #####
#################-#


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

# both distributions
sims.mean <- bind_rows(list(Poisson = sims.mean.pois, 
                            Binomial = sims.mean.bin),
                       .id = "model")


#################-#
##### Figures #####
#################-#

pbin <- ggplot(final.bin, aes(y=ks.sig/100, x=as.factor(controlValues),
                              col=as.factor(intercept))) +
  geom_point(position = position_dodge(width=0.4)) +
  geom_line(aes(x=as.numeric(as.factor(controlValues))),
            position = position_dodge(width=0.4))+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up),width = 0.1,
                position = position_dodge(width=0.4)) +
  geom_hline(yintercept = 0.05, linetype="dashed") +
  ylim(0,1)+
  scale_color_manual("intercept",values = col.intercept)+
  xlab("sampleSize") + ylab("Prop of significant KS test") +
  labs(title="Binomial", tag = "B)") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.7,0.8),
        legend.box.background = element_rect(fill="gray94",color="gray94"))
pbin

ggplot(final.pois, aes(y=ks.sig/100, x=as.factor(controlValues),
                       col=as.factor(intercept))) +
  geom_point(position = position_dodge(width=0.4)) +
  geom_line(aes(x=as.numeric(as.factor(controlValues))),
            position = position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up),width = 0.1,
                position = position_dodge(width=0.4)) +
  geom_hline(yintercept = 0.05, linetype="dashed") +
  scale_color_manual("intercept",values = col.intercept)+
  xlab("sampleSize") + ylab("Prop of significant KS test") +
  labs(title="Poisson", tag = "A)") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none") +
  pbin + plot_layout(ncol=2)
ggsave(here("figures", "1_glmBOTH_pearsonChisq.jpeg"), width = 9, height = 5)


# distributions

sims.mean %>% filter(controlValues == 10) %>%
  ggplot(aes(x=mean.pearson, col=model)) +
  facet_grid(controlValues~intercept)+
  stat_density(geom="line",size=1.5, position="identity")+
  stat_function(fun = dchisq, args = list(df = 8), aes(color = "Chi-squared")) +
  scale_color_manual(values = c("Binomial" ="olivedrab4", 
                                "Poisson" = "mediumpurple3", 
                                "Chi-squared" = "black"))+
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "inside",
        legend.box.background = element_rect(fill="gray94", color="gray94"),
        legend.position.inside  = c(0.87,0.85),
        legend.key = element_rect(),
        plot.margin = margin(-20, 5, -10, 5)) +
  
  sims.mean %>% filter(controlValues == 20) %>%
  ggplot(aes(x=mean.pearson, col=model)) +
  facet_grid(controlValues~intercept)+
  geom_density(size=1.5) +
  scale_color_manual(values = c("olivedrab4", "mediumpurple3"))+
  stat_function(fun = dchisq, args = list(df = 18), col="black") +
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(-20, 5, -10, 5)) +
  
  sims.mean %>% filter(controlValues == 50) %>%
  ggplot(aes(x=mean.pearson, col=model)) +
  facet_grid(controlValues~intercept)+
  geom_density(size=1.5) +
  scale_color_manual(values = c("olivedrab4", "mediumpurple3"))+
  stat_function(fun = dchisq, args = list(df = 48), col="black") +
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(-20, 5, -10, 5)) +
  
  sims.mean %>% filter(controlValues == 100) %>%
  ggplot(aes(x=mean.pearson, col=model)) +
  facet_grid(controlValues~intercept)+
  geom_density(size=1.5) +
  scale_color_manual(values = c("olivedrab4", "mediumpurple3"))+
  stat_function(fun = dchisq, args = list(df = 98), col="black") +
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(-20, 5, -10, 5)) +
  
  sims.mean %>% filter(controlValues == 200) %>%
  ggplot(aes(x=mean.pearson, col=model)) +
  facet_grid(controlValues~intercept)+
  geom_density(size=1.5) +
  scale_color_manual(values = c("olivedrab4", "mediumpurple3"))+
  stat_function(fun = dchisq, args = list(df = 198), col="black") +
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(-20, 5, -10, 5)) +
  
  sims.mean %>% filter(controlValues == 500) %>%
  ggplot(aes(x=mean.pearson, col=model)) +
  facet_grid(controlValues~intercept)+
  geom_density(size=1.5) +
  scale_color_manual(values = c("olivedrab4", "mediumpurple3"))+
  stat_function(fun = dchisq, args = list(df = 498), col="black") +
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(-20, 5, -10, 5)) +
  
  sims.mean %>% filter(controlValues == 1000) %>%
  ggplot(aes(x=mean.pearson, col=model)) +
  facet_grid(controlValues~intercept)+
  geom_density(size=1.5) +
  scale_color_manual(values = c("olivedrab4", "mediumpurple3"))+
  stat_function(fun = dchisq, args = list(df = 998), col="black") +
  xlab("") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(-20, 5, -20, 5)) +
  
  sims.mean %>% filter(controlValues == 10000) %>%
  ggplot(aes(x=mean.pearson, col=model)) +
  facet_grid(controlValues~intercept)+
  geom_density(size=1.5) +
  scale_color_manual(values = c("olivedrab4", "mediumpurple3"))+
  stat_function(fun = dchisq, args = list(df = 9998), col="black") +
  xlab("Pearson Stats") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1),
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(-20, 5, 0, 5)) +
  
  plot_layout(ncol=1) + 
  plot_annotation(title="Pearson Statistics X Chi-squared distribution")

ggsave(here("figures", "1_glm_pearsonChisq_distrib_MEAN.jpeg"), width=12,
       height=16)














