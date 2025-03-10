### Dispersion Tests Project
## Melina Leite
# Mar 25

library(DHARMa)
library(tidyverse)
library(here)

# parameters
overdispersion <- seq(0, 1, 0.10)
sampleSize = 500
intercept <- 0
slope = 1
ntrials = c(5, 20)
nRep = 1000

#################-##
##### Binomial #####
#################-##


trial.bin <- list()

for (k in ntrials){
  
  # function to varying sampleSize
  calculateStatistics <- function(control = 0){
    # data
    testData <- createData(overdispersion = control,
                           fixedEffects = slope,
                           sampleSize = sampleSize,
                           intercept = intercept,
                           numGroups = 10,
                           randomEffectVariance = 0,
                           binomialTrials = k,
                           family = binomial())
    # model
    fittedModel <-glm(cbind(observedResponse1, observedResponse0)  ~
                        Environment1, data = testData, family = binomial())
    #results
    out <- list()
    
    # pearson residual
    out$Pear.p.val <- testDispersion(fittedModel, plot = F,
                                     type="PearsonChisq")$p.value
    out$Pear.stat <- testDispersion(fittedModel, plot = F,
                                    type="PearsonChisq")$statistic
    # DHARMa default residuals
    res <- simulateResiduals(fittedModel)
    out$DHA.p.val <- testDispersion(res, type = "DHARMa",plot = F)$p.value
    out$DHA.stat  <- testDispersion(res, type = "DHARMa",plot = F)$statistic
    
    # DHARMa refit residuals -> bootstrapped Pearson
    res <- simulateResiduals(fittedModel, refit=T)
    out$Ref.p.val <- testDispersion(res, plot = F, type = "DHARMa")$p.value
    out$Ref.stat <- testDispersion(res, plot = F, type = "DHARMa")$statistic
    
    return(unlist(out))
  }
  
  out <- runBenchmarks(calculateStatistics, controlValues = overdispersion,
                       nRep=nRep, parallel = T, exportGlobal = T)
 trial.bin[[length(trial.bin) + 1]] <- out
  
  names(trial.bin)[length(trial.bin)] <- k
  # saving sim results
  save(trial.bin, file=here("data","3b_glmBin_ntrials.Rdata"))
}
