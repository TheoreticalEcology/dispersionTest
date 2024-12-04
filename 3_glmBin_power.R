### Dispersion Tests Project
## Melina Leite
# Dez 24

library(DHARMa)
library(tidyverse)
library(here)

#####################################
##### Instructions & Simulations #####
#####################################

# 1) Simulating 1000 binomial prop datasets with different sample sizes, intercepts and overdispersion. Fixing the number of trials to 10:
#       - sampleSize: c(10,50,100,500)
#       - intercept:  c(-3,-1,0,1,3)
#       - overdisperstion: seq(0,1,0.10)
# 2) fitting them to correct GLM models
# 5) calculating power for the dispersion tests 
#     - Pearson-chisq, 
#     - DHARMa default (simulated residuals)
#     - DHARMa refit (boostrapped Pearson residuals)



# varying parameters
overdispersion <- seq(0,1,0.10)
sampleSize = c(10,50,100,500)
intercept <- c(-3,-1,0,1,3)


out.out <- list()
for (k in intercept){
  for (i in sampleSize){
  
  # function to varying sampleSize
  calculateStatistics <- function(control = 0){
    # data
    testData <- DHARMa::createData(overdispersion = control,
                                   sampleSize = k,
                                   intercept = i,
                                   numGroups = 10,
                                   randomEffectVariance = 0,
                                   binomialTrials = 10,
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
    out$DHA.p.val <- testDispersion(res, type = "DHARMa",plot = F)$p.value
    
    # DHARMa refit residuals -> bootstrapped Pearson
    res <- simulateResiduals(fittedModel, refit=T)
    out$Ref.p.val <- testDispersion(res, plot = F, type = "DHARMa")$p.value
    return(unlist(out))
  }
  
  out <- runBenchmarks(calculateStatistics, controlValues = overdispersion,
                       nRep=1000, parallel = T, exportGlobal = T)
  out.out[[length(out.out) + 1]] <- out
}
}
#names(out.out) <- sampeSize

# saving sim results
save(out.out,sampleSize,intercept,overdispersion, file=here("data", 
                                             "3_glmBin_power.Rdata"))


