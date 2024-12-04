### Dispersion Tests Project
## Melina Leite
# Dez 24

library(DHARMa)
library(tidyverse)
library(here)

#####################################
##### Instructions & Simulations #####
#####################################

# 1) Simulating 1000 binomial prop datasets with different sample sizes and overdispersion. Fixing the number of trials to 20:
#       - sampleSize: c(10,50,100,500)
#       - overdisperstion: seq(0,1,0.10)
# 2) fitting them to correct GLM models
# 5) calculating power for the dispersion tests 
#     - Pearson-chisq, 
#     - DHARMa default (simulated residuals)
#     - DHARMa refit (boostrapped Pearson residuals)


# function to varying sampleSize
calculateStatistics <- function(control = 10){
  # data
  testData <- DHARMa::createData(sampleSize = 10,
                                 overdispersion = control,
                                 numGroups = 1,
                                 binomialTrials = 20,
                                 family = binomial())
  # model
  fittedModel <- stats::glm(cbind(observedResponse1, observedResponse0)  ~ 
                              Environment1, data = testData, family = binomial()) 
  #results
  out <- list()
  
  # pearson residual
  out$Pear.p.val <- testDispersion(fittedModel, plot = F, 
                                   type="PearsonChisq")$p.value
  # DHARMa default residuals
  res <- simulateResiduals(fittedModel)
  out$DHA.p.val<- testDispersion(res, type = "DHARMa",plot = F)$p.value
  
  # DHARMa refit residuals -> bootstrapped Pearson
  res <- simulateResiduals(fittedModel, refit=T)
  out$Ref.p.val <- testDispersion(res, plot = F, type = "DHARMa")$p.value
  return(unlist(out))
}
#calculateStatistics()

# varying intercept in a loop 
overdispersion <- seq(0,1,0.10)
sampleSize = c(10,50,100,500)

out.out <- list()
for (i in sampleSize){
  
  out <- runBenchmarks(calculateStatistics, controlValues = overdispersion,
                       intercept = 0,
                       nRep=1000, parallel = T)
  out.out[[length(out.out) + 1]] <- out
}

names(out.out) <- sampeSize

# saving sim results
save(out.out,sampleSize,overdispersion, file=here("data", 
                                             "4_glmBin_power.Rdata"))


