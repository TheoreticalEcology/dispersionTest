### Dispersion Tests Project
## Melina Leite
# jan 25

library(DHARMa)
library(here)
library(tidyverse)
library(cowplot); theme_set(theme_cowplot())




# seeing data:
load(here("data","7_alternatives_pois.Rdata"))

simpois <- map_dfr(out.pois, "simulations", .id="ngroups")  %>%
  separate(ngroups, c("sampleSize", "intercept")) %>%
  rename("nSim" = "controlValues") 

# looking at the number of simulations with at least one zero SD
simpois %>% group_by(intercept,sampleSize, nSim) %>%
  summarise(pro.zero = sum(zeroSD)/n())

# looking at the percentage of zeros in the simulations

simpois %>%
  ggplot(aes(x=intercept, y=prop.zeroSD, col=as.factor(nSim))) +
  geom_boxplot()+
  facet_wrap(~sampleSize, scales="free")


## evaluating just the results with NO zero SD-obs

p.pois <- simpois %>% filter(zeroSD ==0) %>% 
  select(sampleSize,intercept, nSim, ends_with(".p")) %>%
  select(-c(A2.p, A3.p, A4.p)) %>%
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
  ggplot(aes(x=nSim, y=prop.sig, col=intercept, linetype=sampleSize))+
  facet_grid(~test)+
  geom_point(position = position_dodge(width = 0.2)) + 
  geom_line(position = position_dodge(width = 0.2), aes(x=as.numeric(nSim))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up), width=0.02,
                position = position_dodge(width = 0.2)) +
  scale_y_sqrt(breaks=c(0,0.01, 0.05,0.25,0.5, 0.75)) +
  geom_hline(yintercept = 0.05, linetype = "dotted")





#dispersion stat
dis <- simpois %>%  filter(zeroSD ==0) %>%
  select(sampleSize, intercept, nSim, ends_with("dispersion"), A1.stat) %>%
  pivot_longer(4:6,names_to = "test", values_to = "disp.statistics") %>%
  mutate(nSim = as.factor(nSim))

dis %>% 
  ggplot(aes(x = nSim, y = disp.statistics, col = test))+
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype="dashed")  +
  scale_y_sqrt()


# only simulations with at least one zero obs
sim0 <- simpois %>% filter(zeroSD == 1)

#dispersion stat
dis <- sim0 %>% select(sampleSize,intercept,nSim, ends_with("dispersion"),
                       ends_with(".stat")) %>%
  pivot_longer(4:8,names_to = "test", values_to = "disp.statistics")

dis %>% filter(disp.statistics <10) %>%
  ggplot(aes(x = test, y = disp.statistics, col = sampleSize))+
  geom_hline(yintercept = 1, linetype="dashed")+
  geom_boxplot() +

  facet_grid(nSim ~ intercept, scales = "free") +
  scale_y_sqrt()
