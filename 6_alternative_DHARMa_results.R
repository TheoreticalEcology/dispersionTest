### Dispersion Tests Project
## Melina Leite
# jan 25

library(DHARMa)
library(here)
library(tidyverse)
library(cowplot); theme_set(theme_cowplot())

# plot Colors
source(here("plotColors.R"))

# Simulated data ####
load(here("data","6_approxPear_pois.Rdata"))

simpois <- map_dfr(out.pois, "simulations", .id="ngroups")  %>%
  separate(ngroups, c("nSim", "sampleSize", "intercept"), sep="_") %>%
  rename("overdispersion" = "controlValues") %>%
  mutate(nSim = as.factor(as.numeric(nSim)))

load(here("data","6_approxPear_bin.Rdata"))

simbin <- map_dfr(out.bin, "simulations", .id="ngroups")  %>%
  separate(ngroups, c("nSim", "sampleSize", "intercept"), sep="_") %>%
  rename("overdispersion" = "controlValues") %>%
  mutate(nSim = as.factor(as.numeric(nSim)))





# Percentage of zeros ####
# looking at the percentage of zeros in the simulations

simpois %>%
  ggplot(aes(x=as.factor(overdispersion), y=prop.zero, col=as.factor(sampleSize))) +
  geom_boxplot()+
  facet_grid(nSim~intercept, scales="free") +
  labs(title="Poisson: Prop obs with zero SD simulated data") +
  theme(panel.background = element_rect(color="black"))
#ggsave(here('figures', '6_propzero.jpeg'), height=10, width=10)

simbin %>%
  ggplot(aes(x=as.factor(overdispersion), y=prop.zero, col=as.factor(sampleSize))) +
  geom_boxplot()+
  facet_grid(nSim~intercept, scales="free") +
  labs(title="Binomial: Prop obs with zero SD simulated data") +
  theme(panel.background = element_rect(color="black"))
ggsave(here('figures', '6_propzero_bin.jpeg'), height=10, width=10)


# evaluating just the results with NO zero SD-obs ####

## type I error ####

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
#ggsave(here('figures', '6_type1.jpeg'), height=6, width=10)

p.bin <- simbin %>% filter(prop.zero == 0, overdispersion == 0) %>% 
  select(sampleSize,intercept, nSim, ends_with(".p")) %>%
  pivot_longer(4:6, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize, intercept, nSim, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]))
p.bin$prop.sig <- p.bin$p.sig/p.bin$nsim
p.bin$nSim <- as.factor(p.bin$nSim)
# add sig test
for (i in 1:nrow(p.bin)) {
  btest <- binom.test(p.bin$p.sig[i], n=p.bin$nsim[i], p=0.05)
  p.bin$p.bin0.05[i] <- btest$p.value
  p.bin$conf.low[i] <- btest$conf.int[1]
  p.bin$conf.up[i] <- btest$conf.int[2]
}

p.bin %>%
  ggplot(aes(x=nSim, y=prop.sig, col=intercept))+
  facet_grid(sampleSize~test)+
  geom_point(position = position_dodge(width = 0.2)) + 
  geom_line(position = position_dodge(width = 0.2), aes(x=as.numeric(nSim))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up), width=0.02,
                position = position_dodge(width = 0.2)) +
  scale_y_sqrt(breaks=c(0,0.01, 0.05,0.25,0.5, 0.75)) +
  geom_hline(yintercept = 0.05, linetype = "dotted")+
  labs(title= "Type I error Binomial") +
  theme(panel.background = element_rect(color="black"))
#ggsave(here('figures', '6_type1_bin.jpeg'), height=6, width=10)




# dispersion stat for zeo overdisp ####
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
#ggsave(here('figures', '6_dispersion.jpeg'), height=6, width=10)


disbin <- simbin %>%  filter(prop.zero ==0, overdispersion==0) %>%
  select(sampleSize, intercept, nSim, ends_with("dispersion")) %>%
  pivot_longer(4:6,names_to = "test", values_to = "disp.statistics") %>%
  mutate(nSim = as.factor(nSim))

disbin %>% 
  ggplot(aes(x = nSim, y = disp.statistics, col = intercept))+
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype="dashed")  +
  scale_y_log10() +
  facet_grid(sampleSize~test, scales="free")+
  labs(title= "Binomial: Dispersion statistics") +
  theme(panel.background = element_rect(color="black"))
#ggsave(here('figures', '6_dispersionBin.jpeg'), height=6, width=10)



## differences between dispersions stats and Pearson ####
dis2 <- simpois %>% filter(prop.zero ==0, overdispersion==0) %>%
  mutate(P.DHA = Pear.stat.dispersion - DHA.stat.dispersion,
                            P.Alt = Pear.stat.dispersion - Alt.stat.dispersion) %>% 
  select(sampleSize,intercept,nSim, P.DHA, P.Alt) %>%
  pivot_longer(4:5,names_to = "test", values_to = "dif")

#library(patchwork)
dis2 %>% 
  ggplot(aes(x = as.factor(nSim), y = dif, col = intercept))+
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype="dashed")  +
  facet_grid(sampleSize~test, scales="free")+
  labs(title= "Dispersion statistics") +
  theme(panel.background = element_rect(color="black")) #+

dis2 %>% filter(nSim %in% c(250,1000)) %>%
  ggplot(aes(x = as.factor(nSim), y = dif, col = intercept))+
  geom_boxplot() +
  ylim(-0.1,0.1) + ylab("Dif (Pearson - disp)")+
   geom_hline(yintercept = 0, linetype="dashed")  +
  facet_grid(sampleSize~test, scales="free")+
  labs(title= "Dispersion statistics") +
  theme(panel.background = element_rect(color="black")) #+
#  plot_layout(ncol=1)
#ggsave(here('figures', '6_dispersion_dif.jpeg'), height=12, width=10)


#binomial
disbin2 <- simbin %>% filter(prop.zero ==0, overdispersion==0) %>%
  mutate(P.DHA = Pear.stat.dispersion - DHA.stat.dispersion,
         P.Alt = Pear.stat.dispersion - Alt.stat.dispersion) %>% 
  select(sampleSize,intercept,nSim, P.DHA, P.Alt) %>%
  pivot_longer(4:5,names_to = "test", values_to = "dif")

#library(patchwork)
disbin2 %>% 
  ggplot(aes(x = as.factor(nSim), y = dif, col = intercept))+
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype="dashed")  +
  facet_grid(sampleSize~test, scales="free")+
  labs(title= "Binomial: Dispersion statistics") +
  theme(panel.background = element_rect(color="black")) #+

disbin2 %>% filter(nSim %in% c(250,1000)) %>%
  ggplot(aes(x = as.factor(nSim), y = dif, col = intercept))+
  geom_boxplot() +
  ylim(-0.1,0.1) + ylab("Dif (Pearson - disp)")+
  geom_hline(yintercept = 0, linetype="dashed")  +
  facet_grid(sampleSize~test, scales="free")+
  labs(title= "Binomial: Dispersion statistics") +
  theme(panel.background = element_rect(color="black")) #+
 # plot_layout(ncol=1)
ggsave(here('figures', '6_dispersion_difBin.jpeg'), height=12, width=10)



# dipersion with oeverdispersion
disper <- simpois %>%  filter(prop.zero ==0) %>%
  select(sampleSize, intercept, nSim, ends_with("dispersion")) %>%
  pivot_longer(4:6,names_to = "test", values_to = "disp.statistics") %>%
  mutate(nSim = as.factor(nSim))


disper %>% group_by(overdispersion,sampleSize, nSim, test, intercept) %>%
  summarise(mean.disp = mean(disp.statistics)) %>%
  ggplot(aes(x = overdispersion, y = mean.disp, col = test, linetype=intercept))+
  geom_line() +
  geom_hline(yintercept = 1, linetype="dashed")  +
  scale_y_log10() +
  facet_grid(sampleSize~nSim, scales="free")+
  labs(title= "Dispersion statistics") +
  theme(panel.background = element_rect(color="black"))
#ggsave(here('figures', '6_dispersion_overdisp.jpeg'), height=6, width=12)

disperbin <- simbin %>%  filter(prop.zero ==0) %>%
  select(sampleSize, intercept, nSim, ends_with("dispersion")) %>%
  pivot_longer(4:6,names_to = "test", values_to = "disp.statistics") %>%
  mutate(nSim = as.factor(nSim))


disperbin %>% group_by(overdispersion,sampleSize, nSim, test, intercept) %>%
  summarise(mean.disp = mean(disp.statistics)) %>%
  filter(nSim %in% c(250, 1000)) %>%
  ggplot(aes(x = overdispersion, y = mean.disp, col = test, linetype=intercept))+
  geom_line() +
  geom_hline(yintercept = 1, linetype="dashed")  +
  scale_y_log10() +
  facet_grid(sampleSize~nSim, scales="free")+
  labs(title= "Binomial: Dispersion statistics") +
  theme(panel.background = element_rect(color="black"))
ggsave(here('figures', '6_dispersion_overdispBin.jpeg'), height=6, width=8)




# power ####
power <- simpois %>% filter(prop.zero == 0,) %>% 
  select(sampleSize,intercept, nSim, overdispersion, ends_with(".p")) %>%
  pivot_longer(5:7, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize, intercept, nSim, test, overdispersion) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]))
power$prop.sig <- power$p.sig/power$nsim
power$nSim <- as.factor(power$nSim)
# add sig test
for (i in 1:nrow(power)) {
  btest <- binom.test(power$p.sig[i], n=power$nsim[i], p=0.05)
  power$p.bin0.05[i] <- btest$p.value
  power$conf.low[i] <- btest$conf.int[1]
  power$conf.up[i] <- btest$conf.int[2]
}

power %>% 
  ggplot(aes(x = overdispersion, y = prop.sig, col = test, linetype=intercept))+
  geom_line() +
  geom_hline(yintercept = 0)  +
  facet_grid(sampleSize~nSim, scales="free")+
  labs(title= "Power") +
  theme(panel.background = element_rect(color="black"))
#ggsave(here('figures', '6_power.jpeg'), height=6, width=12)



powerbin <- simbin %>% filter(prop.zero == 0,) %>% 
  select(sampleSize,intercept, nSim, overdispersion, ends_with(".p")) %>%
  pivot_longer(5:7, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize, intercept, nSim, test, overdispersion) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]))
powerbin$prop.sig <- powerbin$p.sig/powerbin$nsim
powerbin$nSim <- as.factor(powerbin$nSim)
# add sig test
for (i in 1:nrow(powerbin)) {
  btest <- binom.test(powerbin$p.sig[i], n=powerbin$nsim[i], p=0.05)
  powerbin$p.bin0.05[i] <- btest$p.value
  powerbin$conf.low[i] <- btest$conf.int[1]
  powerbin$conf.up[i] <- btest$conf.int[2]
}

powerbin %>% 
  ggplot(aes(x = overdispersion, y = prop.sig, col = test, linetype=intercept))+
  geom_line() +
  geom_hline(yintercept = 0)  +
  facet_grid(sampleSize~nSim, scales="free")+
  labs(title= "Power") +
  theme(panel.background = element_rect(color="black"))
ggsave(here('figures', '6_powerBin.jpeg'), height=6, width=12)


# lm approximate pearson and the pearson residuals ####

#lm(alterna$resPA - resP ~ resP)
# expecting intercept and slope in zero

# slope
slope <- simpois %>%  filter(prop.zero == 0) %>%
  select(sampleSize, intercept, nSim, overdispersion, resPA.slope, resPA.slopeP) %>%
  mutate(nSim = as.factor(nSim))

slope %>%
  ggplot(aes(x=as.factor(overdispersion), y=resPA.slope, col=intercept))+
  geom_boxplot() +
  facet_grid(sampleSize ~ nSim, scales="free") +
  geom_hline(yintercept = 0, linetype="dotted")+
  theme(panel.background = element_rect(color = "black"))

slope %>% filter(nSim %in% c("250","1000")) %>%
  ggplot(aes(x=as.factor(overdispersion), y=resPA.slope, col=intercept))+
  geom_boxplot() +
  facet_grid(sampleSize ~ nSim) +
  geom_hline(yintercept = 0, linetype="dotted")+
  theme(panel.background = element_rect(color = "black"))

slope %>% filter(nSim %in% c("250","1000")) %>%
  group_by(nSim, sampleSize, intercept, overdispersion) %>%
  summarise(mean.slope = mean(resPA.slope),
            median.slope = median(resPA.slope)) %>%
  ggplot(aes(x=overdispersion, y=mean.slope, col=intercept)) +
  geom_line()+ geom_point()+
  facet_grid(sampleSize ~ nSim)   +
  geom_hline(yintercept = 0, linetype="dotted")+
  theme(panel.background = element_rect(color = "black"))

# binomial
library(patchwork)
slopebin <- simbin %>%  filter(prop.zero == 0) %>%
  select(sampleSize, intercept, nSim, overdispersion, resPA.slope, resPA.slopeP) %>%
  mutate(nSim = as.factor(nSim))

slopebin %>% filter(nSim %in% c("250","1000")) %>%
  ggplot(aes(x=as.factor(overdispersion), y=resPA.slope, col=intercept))+
  geom_boxplot() +
  facet_grid(sampleSize ~ nSim) +
  theme(panel.background = element_rect(color = "black")) 


slopebin %>% filter(nSim %in% c("250","1000")) %>%
  group_by(nSim, sampleSize, intercept, overdispersion) %>%
  summarise(mean.slope = mean(resPA.slope),
            median.slope = median(resPA.slope)) %>%
  ggplot(aes(x=overdispersion, y=mean.slope, col=intercept)) +
  geom_line()+ geom_point()+
  facet_grid(sampleSize ~ nSim)   +
  labs(title="Binomial: mean slope", subtitle = "(App.Pear - Pear) ~ Pear")+
  theme(panel.background = element_rect(color = "black")) 



# intercept
inter <- simpois %>%  filter(prop.zero ==0) %>%
  select(sampleSize, intercept, nSim, overdispersion, resPA.inter, resPA.interP) %>%
  mutate(nSim = as.factor(nSim))

inter %>%
  ggplot(aes(x=as.factor(overdispersion), y=resPA.inter, col=intercept))+
  geom_boxplot() +
  facet_grid(sampleSize ~ nSim, scales="free") +
  geom_hline(yintercept = 0, linetype="dotted")+
  theme(panel.background = element_rect(color = "black"))

inter %>% filter(nSim %in% c("250","1000")) %>%
  ggplot(aes(x=as.factor(overdispersion), y=resPA.inter, col=intercept))+
  geom_boxplot() +
  facet_grid(sampleSize ~ nSim) +
  geom_hline(yintercept = 0, linetype="dotted")+
  theme(panel.background = element_rect(color = "black"))

inter %>% filter(nSim %in% c("250","1000")) %>%
  group_by(nSim, sampleSize, intercept, overdispersion) %>%
  summarise(mean.inter = mean(resPA.inter),
            median.inter = median(resPA.inter)) %>%
  ggplot(aes(x=overdispersion, y=mean.inter, col=intercept)) +
  geom_line()+ geom_point()+
  facet_grid(sampleSize ~ nSim)   +
  geom_hline(yintercept = 0, linetype="dotted")+
  theme(panel.background = element_rect(color = "black"))


# binomial
interbin <- simbin %>%  filter(prop.zero ==0) %>%
  select(sampleSize, intercept, nSim, overdispersion, resPA.inter, resPA.interP) %>%
  mutate(nSim = as.factor(nSim))

inter %>% filter(nSim %in% c("250","1000")) %>%
  group_by(nSim, sampleSize, intercept, overdispersion) %>%
  summarise(mean.inter = mean(resPA.inter),
            median.inter = median(resPA.inter)) %>%
  ggplot(aes(x=overdispersion, y=mean.inter, col=intercept)) +
  geom_line()+ geom_point()+
  facet_grid(sampleSize ~ nSim)   +
  geom_hline(yintercept = 0, linetype="dotted")+
  labs(title="Binomial: mean intercept", subtitle = "(App.Pear - Pear) ~ Pear")+
  theme(panel.background = element_rect(color = "black"))







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

