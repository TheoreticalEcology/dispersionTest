### Dispersion Tests Project
## Melina Leite
# Dez 24

library(DHARMa)
library(lme4)
library(tidyverse); 
library(here)





#####################################
##### Instructions & Simulations #####
#####################################

# GLMM Poisson tests of type 1 and power for 5 dispersion tests

# 1) Simulating 100 GLMM Poisson datasets with different sample sizes, intercepts, RE-intercept variance, and different dispersion levels:
#       - sampleSize (3): c(50, 100, 500)
#       - intercept (4):  c(-1, 0, 2, 4)
#       - RE variance (3): c(0.5, 1, 1.5)
#       - dispersion (6): c(0, 0.1, 0.2, 0.3, 0.4, 0.6)
# 2) fitting them to correct GLMM models
# 5) calculating type I error rate and power for: 
#     - Pearson-chisq naive DF 
#     - Bootstrapped pearson with UNCONDITIONAL simulations
#     - Bootstrapped pearson with CONDITIONAL simulations
#     - DHARMa residuals with UNCONDITIONAL simulations
#     - DHARMa residuals with CONDITIONAL simulations


# values
sampleSize <- c(50, 100, 500)
intercept <- c(-1, 0, 2, 4)
REvariance <- c(0.5, 1, 1.5)
dispersion <- c(0, 0.1, 0.2, 0.3, 0.4, 0.6)


## function varying dispersion 
# function to varying sampleSize
calculateStatistics <- function(control = 0){
  # data
  testData <- DHARMa::createData(sampleSize = 50, numGroups = 10,
                                 intercept = -1,
                                 randomEffectVariance = 1,
                                 overdispersion = control,
                                 family = poisson())
  # model
  fittedModel <- lme4::glmer(observedResponse ~ Environment1 + (1|group), 
                            data = testData, family = poisson())
  #results
  out <- list()
  
  # pearson residual
  out$pear <- testDispersion(fittedModel, plot = F, alternative = "greater",
                                   type="PearsonChisq")$p.value
  
  # Bootstrapped pearson with UNCONDITIONAL simulations
  resUN <- simulateResiduals(fittedModel, refit=T)
  out$refitUN <- testDispersion(resUN, type = "DHARMa",plot = F)$p.value
  
  # Bootstrapped pearson with CONDITIONAL simulations
  resCO <- simulateResiduals(fittedModel, refit=T, re.form=NULL)
  out$refitCO <- testDispersion(resCO, plot = F, type = "DHARMa")$p.value
  
  # DHARMa residuals with UNCONDITIONAL simulations
  resUN <- simulateResiduals(fittedModel, refit = F, re.form=NA)
  out$dhaUN <- testDispersion(resUN, plot = F, type = "DHARMa")$p.value
  
  #DHARMa residuals with CONDITIONAL simulations
  resCO <- simulateResiduals(fittedModel, refit = F, re.form=NULL)
  out$dhaCO <- testDispersion(resCO, plot = F, type = "DHARMa")$p.value
  
  
  return(unlist(out))
}
#calculateStatistics()

####################
##### Running #####
###################

# k=0;j=50;i=0.5

output <- list()

for (i in REvariance){ # Varying RE variances

  for (j in sampleSize){ # Varying sample sizes
    
    for (k in intercept){ # Varying intercepts
     
      out <- runBenchmarks(calculateStatistics, controlValues = dispersion,
                           intercept = k,
                           sampleSize = j,
                           REvariance = i,
                           nRep=100, parallel = T)
      output[[length(output) + 1]] <- out

    }
  }
}

nms <- expand.grid(inter = intercept,
            sample = sampleSize,
            revar = REvariance) %>%
        unite("name")

names(output) <- as.vector(nms)$name

save(output, file=here("data", "3_glmmPois_dispersionTests.Rdata"))




####################
##### Results #####
###################

library(cowplot);
theme_set(theme_cowplot() +
           theme(panel.background=element_rect(color="black")))
library(patchwork)

load(here("data", "3_glmmPois_dispersionTests.Rdata"))

props <- list()
for (i in 1:length(output)){
  sims <- strsplit(names(output)[i], "_")[[1]]
  suma <- output[[i]]$summaries$propSignificant
  suma$intercept = sims[1]
  suma$sampleSize = sims[2]
  suma$REvariance = sims[3]
  props <- rbind(props,suma)
  }
propSig <- props %>% pivot_longer(cols = 2:6, names_to = "test",
                                  values_to = "propSig") %>%
  mutate(sampleSize = fct_relevel(sampleSize, "50", "100", "500"))

propSig %>% filter(REvariance == 0.5) %>%
ggplot(aes(x = controlValues, y= propSig, col = test))+
  facet_grid(sampleSize ~ intercept) +
  geom_point() + geom_line() + 
  ylim(0,0.35) +
propSig %>% filter(REvariance == 1) %>%
  ggplot(aes(x = controlValues, y= propSig, col = test))+
  facet_grid(sampleSize ~ intercept) +
  geom_point() + geom_line() +
  ylim(0,0.35) +
propSig %>% filter(REvariance == 1.5) %>%
  ggplot(aes(x = controlValues, y= propSig, col = test))+
  facet_grid(sampleSize ~ intercept) +
  geom_point() + geom_line() +
  ylim(0,0.35) +
plot_layout(ncol=2, guides="collect")


propSig %>% filter(test == "dhaCO") %>%
  ggplot(aes(x = controlValues, y= propSig, col = intercept))+
  facet_grid(sampleSize ~ REvariance) +
  geom_point() + geom_line() + 
  ylim(0,0.35) +
propSig %>% filter(test == "dhaUN") %>%
  ggplot(aes(x = controlValues, y= propSig, col = intercept))+
  facet_grid(sampleSize ~ REvariance) +
  geom_point() + geom_line() + 
  ylim(0,0.35) +
propSig %>% filter(test == "refitCO") %>%
  ggplot(aes(x = controlValues, y= propSig, col = intercept))+
  facet_grid(sampleSize ~ REvariance) +
  geom_point() + geom_line() + 
  ylim(0,0.35) +
propSig %>% filter(test == "refitUN") %>%
  ggplot(aes(x = controlValues, y= propSig, col = intercept))+
  facet_grid(sampleSize ~ REvariance) +
  geom_point() + geom_line() + 
  ylim(0,0.35) +
plot_layout(ncol=2, guides="collect")

#ignoring intercept

propSig %>% filter(intercept == 0) %>%
  ggplot(aes(x = controlValues, y= propSig, col = test))+
  facet_grid(sampleSize ~ REvariance, 
             labeller = as_labeller(c(`0.5` = "REvar = 0.5",
                                    `1`   = "REvar = 1",
                                    `1.5` = "REvar = 1.5",
                                    `50`  = "N=50; gr=10",
                                    `100` = "N=100; gr=10",
                                    `500` = "N=500; gr=20"))) +
  geom_point() + geom_line() + 
  xlab("Overdispersion") + ylab("Prop")+
  ylim(0,0.35) +
  ggtitle("Poisson GLMM: power 5 dispersion tests",
          subtitle = "100 sim; intercept=0") +
  theme(panel.grid.major = element_line(color="gray",linewidth=0.1) )
ggsave(here("figures", "3_glmmPois_dispersion.jpeg"))























