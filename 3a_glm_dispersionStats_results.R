### Dispersion Tests Project
## Melina Leite
# Jan 25

library(DHARMa)
library(tidyverse)
library(here)


##############-###
#### Binomial ####
##############-###

load(here("data","3a_glmBin_dispersionStat.Rdata"))

#simulations
bin  <- map_dfr(slope.bin, "simulations", .id="slope") %>%  
  rename("overdispersion" = "controlValues")



##### p values #####
p.bin <- bin %>% dplyr::select(Pear.p.val,DHA.p.val,Ref.p.val, replicate,
                                     overdispersion, slope) %>%
  pivot_longer(1:3, names_to = "test", values_to = "p.val") %>%
  group_by(slope, overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]) )
p.bin$prop.sig <- p.bin$p.sig/p.bin$nsim
for (i in 1:nrow(p.bin)) {
  btest <- binom.test(p.bin$p.sig[i], n=p.bin$nsim[i], p=0.05)
  p.bin$p.bin0.05[i] <- btest$p.value
  p.bin$conf.low[i] <- btest$conf.int[1]
  p.bin$conf.up[i] <- btest$conf.int[2]
}
p.bin$test <- factor(p.bin$test, levels = c("Pear.p.val", "Ref.p.val","DHA.p.val"))



##### Dispersion statistics ####
stats.bin <- bin %>% dplyr::select(Pear.stat.dispersion,
                                          DHA.stat.dispersion,
                                          Ref.stat.dispersion, 
                                          overdispersion, slope) %>%
  pivot_longer(1:3, names_to = "test", values_to = "Dispersion") %>%
  group_by(slope, overdispersion,test) %>% 
  summarise(mean.stat = mean(Dispersion, na.rm=T))



##############-###
#### Poisson ####
##############-###

load(here("data","3a_glmPois_dispersionStat.Rdata"))

#simulations
pois  <- map_dfr(slope.pois, "simulations", .id="slope") %>%  
  rename("overdispersion" = "controlValues")



##### p values #####
p.pois <- pois %>% dplyr::select(Pear.p.val,DHA.p.val,Ref.p.val, replicate,
                               overdispersion, slope) %>%
  pivot_longer(1:3, names_to = "test", values_to = "p.val") %>%
  group_by(slope, overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]) )
p.pois$prop.sig <- p.pois$p.sig/p.pois$nsim
for (i in 1:nrow(p.pois)) {
  btest <- binom.test(p.pois$p.sig[i], n=p.pois$nsim[i], p=0.05)
  p.pois$p.pois0.05[i] <- btest$p.value
  p.pois$conf.low[i] <- btest$conf.int[1]
  p.pois$conf.up[i] <- btest$conf.int[2]
}
p.pois$test <- factor(p.pois$test, levels = c("Pear.p.val", "Ref.p.val","DHA.p.val"))



##### Dispersion statistics ####
stats.pois <- pois %>% dplyr::select(Pear.stat.dispersion,
                                   DHA.stat.dispersion,
                                   Ref.stat.dispersion, 
                                   overdispersion, slope) %>%
  pivot_longer(1:3, names_to = "test", values_to = "Dispersion") %>%
  group_by(slope, overdispersion,test) %>% 
  summarise(mean.stat = mean(Dispersion, na.rm=T))



## Figures ####

pval <- bind_rows(list(Poisson = p.pois, Binomial = p.bin), .id="model")

pval %>%
  ggplot(aes(x=overdispersion, y=prop.sig, col=test)) +
  geom_point() + geom_line()+
  facet_grid(model~slope)

statval <- bind_rows(list(Poisson = stats.pois, Binomial = stats.bin), .id="model")


statval %>%
  ggplot(aes(x=overdispersion, y=mean.stat, col=test))+
  geom_point() + geom_line() +
  facet_grid(model~slope, scales="free") +
  labs(title="GLM dispersion stats for dif slopes", 
       subtitle  = "n=500, intercept=0")




