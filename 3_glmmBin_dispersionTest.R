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

# GLMM Binomial tests of type 1 and power for 5 dispersion tests

# 1) Simulating 100 GLMM Binomial datasets with different sample sizes, intercepts, RE-intercept variance, and different dispersion levels. fixing number of trials in 20.
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
                                 binomialTrials = 20,
                                 randomEffectVariance = 1,
                                 overdispersion = control,
                                 family = binomial())
  # model
  fittedModel <- lme4::glmer(cbind(observedResponse1,observedResponse0) ~ 
                               Environment1 + (1|group), 
                             data = testData, family = binomial())
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

save(output, file=here("data", "3_glmmBin_dispersionTests.Rdata"))


