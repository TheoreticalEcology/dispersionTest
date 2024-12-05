### Dispersion Tests Project
## Melina Leite
# Dez 24

library(DHARMa)
library(tidyverse)
library(here)

####################
##### Binomial #####
####################

# 1) Simulating 1000 binomial prop datasets with different sample sizes, intercepts and overdispersion. Fixing the number of trials to 10
# 2) fitting them to correct GLM models
# 5) calculating power for the dispersion tests 
#     - Pearson-chisq, 
#     - DHARMa default (simulated residuals)
#     - DHARMa refit (boostrapped Pearson residuals)



# varying parameters
overdispersion <- seq(0,1,0.10)
sampleSize = c(10,20,50,100,200,500,1000,10000)
intercept <- c(-3,-1.5,0,1.5,3)


out.bin <- list()
for (k in sampleSize){
  for (i in intercept){
    
    # function to varying sampleSize
    calculateStatistics <- function(control = 0){
      # data
      testData <- createData(overdispersion = control,
                             sampleSize = k,
                             intercept = i,
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
                         nRep=1000, parallel = T, exportGlobal = T)
    out.bin[[length(out.bin) + 1]] <- out
  }
}

names(out.bin) <- as.vector(unite(expand.grid(intercept,sampleSize), "sim"))$sim

# saving sim results
save(out.bin,sampleSize,intercept,overdispersion, file=here("data", 
                                                            "3_glmBin_power.Rdata"))



###################
##### Poisson #####
###################
# 1) Simulating 1000 Poisson datasets with different sample sizes, intercepts and overdispertion.
# 2) fitting them to correct GLM models
# 5) calculating power for the dispersion tests 
#     - Pearson-chisq, 
#     - DHARMa default (simulated residuals)
#     - DHARMa refit (boostrapped Pearson residuals)



# varying parameters
overdispersion <- seq(0,1,0.10)
sampleSize = c(10,20,50,100,200,500,1000,10000)
intercept <- c(-3,-1.5,0,1.5,3)


out.pois <- list()
for (k in sampleSize){
  for (i in intercept){
    
    # function to varying overdispersion
    calculateStatistics <- function(control = 0){
      # data
      testData <- DHARMa::createData(overdispersion = control,
                                     sampleSize = k,
                                     intercept = i,
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
                         nRep=1000, parallel = T, exportGlobal = T)
    out.pois[[length(out.pois) + 1]] <- out
  }
}

names(out.pois) <- as.vector(unite(expand.grid(intercept,sampleSize), "sim"))$sim


# saving sim results
save(out.pois,sampleSize,intercept,overdispersion, file=here("data", 
                                                            "3_glmPois_power.Rdata"))