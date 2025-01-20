### Dispersion Tests Project
## Melina Leite
# jan 25

library(DHARMa)
library(here)
library(tidyverse)




# seeing data:
load(here("data","7_alternatives_pois.Rdata"))

simpois <- map_dfr(out.pois, "simulations", .id="ngroups")  %>%
  separate(ngroups, c("sampleSize", "nSim")) %>%
  rename("intercept" = "controlValues") %>%
  mutate(nSim = fct_relevel(nSim, "250", "1000"),
         intercept = fct_relevel(as.character(intercept), "-3", "-1.5", "0"))

# looking at the number of simulations with at least one zero SD
simpois %>% group_by(intercept,nSim, sampleSize ) %>%
  summarise(pro.zero = sum(zeroSD)/n())

# looking at the percentage of zeros in the simulations

simpois %>%
  ggplot(aes(x=intercept, y=prop.zeroSD, col=nSim)) +
  geom_boxplot()+
  facet_wrap(~sampleSize, scales="free")




#dispersion stat
dis <- simpois %>% select(sampleSize,intercept,nSim, ends_with("dispersion"),
                       ends_with(".stat")) %>%
  pivot_longer(4:8,names_to = "test", values_to = "disp.statistics")

dis %>% filter(disp.statistics <100)%>%
  ggplot(aes(x = test, y = disp.statistics, col = sampleSize))+
  geom_hline(yintercept = 1, linetype="dashed")+
  geom_boxplot() +
  
  facet_grid(nSim ~ intercept, scales = "free") +
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
