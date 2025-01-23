### Dispersion Tests Project
## Melina Leite
# jan 25

library(DHARMa)
library(here)
library(tidyverse)
library(cowplot); theme_set(theme_cowplot())




# seeing data:
load(here("data","6_alternative_DHARMA.Rdata"))

simpois <- map_dfr(out.pois, "simulations", .id="ngroups")  %>%
  separate(ngroups, c("nSim", "sampleSize", "intercept"), sep="_") %>%
  rename("overdispersion" = "controlValues") %>%
  mutate(nSim = as.factor(as.numeric(nSim)))

# looking at the percentage of zeros in the simulations

simpois %>%
  ggplot(aes(x=as.factor(overdispersion), y=prop.zero, col=as.factor(sampleSize))) +
  geom_boxplot()+
  facet_grid(nSim~intercept, scales="free") +
  labs(title="Poisson", 
       subtitle="Proportion of obs with zero SD estimated for simulated observation") +
  theme(panel.background = element_rect(color="black"))


## evaluating just the results with NO zero SD-obs

## type I error:
p.pois <- simpois %>% filter(prop.zero == 0, overdispersion == 0) %>% 
  select(sampleSize,intercept, nSim, ends_with(".p")) %>%
  pivot_longer(4:6, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize, intercept, nSim, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]))
p.pois$prop.sig <- p.pois$p.sig/p.pois$nsim
p.pois$nSim <- as.factor(p.pois$nSim)
# add sig test
for (i in 1:nrow(p.pois)) {
  btest <- binom.test(p.pois$p.sig[i], n=p.pois$nsim[i], p=0.05)
  p.pois$p.bin0.05[i] <- btest$p.value
  p.pois$conf.low[i] <- btest$conf.int[1]
  p.pois$conf.up[i] <- btest$conf.int[2]
}



p.pois %>%
  ggplot(aes(x=nSim, y=prop.sig, col=intercept))+
  facet_grid(sampleSize~test)+
  geom_point(position = position_dodge(width = 0.2)) + 
  geom_line(position = position_dodge(width = 0.2), aes(x=as.numeric(nSim))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up), width=0.02,
                position = position_dodge(width = 0.2)) +
  scale_y_sqrt(breaks=c(0,0.01, 0.05,0.25,0.5, 0.75)) +
  geom_hline(yintercept = 0.05, linetype = "dotted")+
  labs(title= "Type I error Poisson") +
  theme(panel.background = element_rect(color="black"))



#dispersion stat for zeo overdisp
dis <- simpois %>%  filter(prop.zero ==0, overdispersion==0) %>%
  select(sampleSize, intercept, nSim, ends_with("dispersion")) %>%
  pivot_longer(4:6,names_to = "test", values_to = "disp.statistics") %>%
  mutate(nSim = as.factor(nSim))

dis %>% 
  ggplot(aes(x = nSim, y = disp.statistics, col = intercept))+
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype="dashed")  +
  scale_y_log10() +
  facet_grid(sampleSize~test, scales="free")+
  labs(title= "Dispersion statistics") +
  theme(panel.background = element_rect(color="black"))

## differences between dispersions stats and Pearson
dis2 <- simpois %>% filter(prop.zero ==0, overdispersion==0) %>%
  mutate(P.DHA = Pear.stat.dispersion - DHA.stat.dispersion,
                            P.Alt = Pear.stat.dispersion - Alt.stat.dispersion) %>% 
  select(sampleSize,intercept,nSim, P.DHA, P.Alt) %>%
  pivot_longer(4:5,names_to = "test", values_to = "dif")

dis2 %>% 
  ggplot(aes(x = as.factor(nSim), y = dif, col = intercept))+
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype="dashed")  +
  facet_grid(sampleSize~test, scales="free")+
  labs(title= "Dispersion statistics") +
  theme(panel.background = element_rect(color="black"))


dis2 %>% filter(nSim %in% c(250,1000)) %>%
  ggplot(aes(x = as.factor(nSim), y = dif, col = intercept))+
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype="dashed")  +
  facet_grid(sampleSize~test, scales="free")+
  labs(title= "Dispersion statistics") +
  theme(panel.background = element_rect(color="black"))



dis2 %>% filter(nSim %in% c(250,1000)) %>%
  ggplot(aes(x = as.factor(nSim), y = dif, col = test))+
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype="dashed")  +
  facet_grid(sampleSize~intercept, scales="free")+
  labs(title= "Dispersion statistics") +
  theme(panel.background = element_rect(color="black"))







# only sim with zeros# only sim with z0_figure1_explanation.Reros
dis0 <- simpois %>% filter(zeroSD != 0) %>%
  select(sampleSize,intercept,nSim, ends_with("dispersion"),
                                               ends_with(".stat")) %>%
  pivot_longer(4:9,names_to = "test", values_to = "disp.statistics")

dis0 %>% #filter(disp.statistics <10) %>%
  ggplot(aes(x = as.factor(nSim), y = disp.statistics, col = test))+
  geom_hline(yintercept = 1, linetype="dashed")+
  geom_boxplot() +
  facet_grid(sampleSize ~ intercept, scales = "free") +
  scale_y_sqrt() +
  theme(panel.background = element_rect(color="black"))

dis0 %>% filter(disp.statistics <10) %>%
  ggplot(aes(x = as.factor(nSim), y = disp.statistics, col = test))+
  geom_hline(yintercept = 1, linetype="dashed")+
  geom_boxplot() +
  facet_grid(sampleSize ~ intercept, scales = "free") +
  scale_y_sqrt() +
  theme(panel.background = element_rect(color="black"))








# all simulations, including the ones with zero SD

#dispersion stat
dis <- simpois %>% select(sampleSize,intercept,nSim, ends_with("dispersion"),
                       ends_with(".stat")) %>%
  pivot_longer(4:9,names_to = "test", values_to = "disp.statistics")

dis %>% #filter(disp.statistics <10) %>%
  ggplot(aes(x = as.factor(nSim), y = disp.statistics, col = test))+
  geom_hline(yintercept = 1, linetype="dashed")+
  geom_boxplot() +
  facet_grid(sampleSize ~intercept, scales = "free") +
  scale_y_sqrt() +
  theme(panel.background = element_rect(color="black"))
 
dis %>% filter(disp.statistics <10) %>%
  ggplot(aes(x = as.factor(nSim), y = disp.statistics, col = sampleSize))+
  geom_hline(yintercept = 1, linetype="dashed")+
  geom_boxplot() +
  facet_grid(intercept ~ test, scales = "free") +
  scale_y_sqrt()+
  theme(panel.background = element_rect(color="black"))

dis %>% filter(nSim %in% c(250,1000)) %>%
  ggplot(aes(x = as.factor(nSim), y = disp.statistics, col = sampleSize))+
  geom_boxplot() +
  facet_grid(intercept ~ test, scales = "free") +
  geom_hline(yintercept = 1, linetype="dashed")+
  scale_y_sqrt()+
  theme(panel.background = element_rect(color="black"))



## absolute difference with Pearson dispersion 

simpois <- simpois %>% mutate(P.DHA = Pear.stat.dispersion - DHA.stat.dispersion,
                   P.A1 = Pear.stat.dispersion - A1.stat,
                   P.A2 = Pear.stat.dispersion - A2.stat,
                   P.A3 = Pear.stat.dispersion - A3.stat,
                   P.A4 = Pear.stat.dispersion - A4.stat)
dis2 <- simpois %>% select(sampleSize,intercept,nSim, P.DHA, P.A1, P.A2, P.A3,
                           P.A4) %>%
  pivot_longer(4:8,names_to = "test", values_to = "dif")

dis2 %>%
  ggplot(aes(x = as.factor(nSim), y = dif, col = test))+
  geom_boxplot() +
  facet_grid(intercept ~ sampleSize, scales = "free") +
  theme(panel.background = element_rect(color="black"))

dis2 %>% filter(nSim %in% c(250,1000)) %>%
  ggplot(aes(x = as.factor(nSim), y = dif, col = test))+
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype="dashed")+
  facet_grid(sampleSize~intercept, scales = "free") +
  theme(panel.background = element_rect(color="black"))

