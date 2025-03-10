### Dispersion Tests Project
## Melina Leite
# Dez 24

library(DHARMa)
library(tidyverse)
library(here)

# parameters
overdispersion <- seq(0, 1, 0.10)
sampleSize = 500
intercept <- 0
slope = c(-2,-1, 0, 1, 2, 3) 
nRep = 1000

#################-##
##### Binomial #####
#################-##

#slope.bin <- list() first sim
#load(here("data","3a_glmBin_dispersionStat.Rdata")) # slope.bin

for (k in slope){

    # function to varying sampleSize
    calculateStatistics <- function(control = 0){
      # data
      testData <- createData(overdispersion = control,
                             fixedEffects = k,
                             sampleSize = sampleSize,
                             intercept = intercept,
                             numGroups = 10,
                             randomEffectVariance = 0,
                             binomialTrials = 10,
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
    slope.bin[[length(slope.bin) + 1]] <- out

    names(slope.bin)[length(slope.bin)] <- k
    # saving sim results
    save(slope.bin, file=here("data","3a_glmBin_dispersionStat.Rdata"))
}



#################-#
##### Poisson #####
################-##


#slope.pois <- list() #first round
#load(here("data","3a_glmPois_dispersionStat.Rdata"))


for (k in slope){

    # function to varying overdispersion
    calculateStatistics <- function(control = 0){
      # data
      testData <- DHARMa::createData(overdispersion = control,
                                     fixedEffects = k,
                                     sampleSize = sampleSize,
                                     intercept = intercept,
                                     numGroups = 10,
                                     randomEffectVariance = 0,
                                     family = poisson())
      # model
      fittedModel <- stats::glm(observedResponse ~ Environment1,
                                data = testData, family = poisson())
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
    slope.pois[[length(slope.pois) + 1]] <- out

    names(slope.pois)[length(slope.pois)] <- k
    # saving sim results
    save(slope.pois, file=here("data","3a_glmPois_dispersionStat.Rdata"))
}



