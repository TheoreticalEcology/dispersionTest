### Dispersion Tests Project
## 
# Dez 24

library(DHARMa)
library(tidyverse)
library(here)
library(lme4)


# parameters
nRep       <- 10000
sampleSize <- 1000
ngroups    <- c(10,20,50,100) 
intercept  <- 0
ntrials    <- 10 # binomial
overdispersion <- seq(0,1,0.10)


#########################
##### Binomial GLMM #####
#########################

out.bin <- list()
for(m in ngroups){
      
  calculateStatistics <- function(control = 0){
    # data
    testData <- createData(overdispersion = control,
                           sampleSize = sampleSize,
                           intercept = intercept,
                           numGroups = m,
                           randomEffectVariance = 1,
                           binomialTrials = ntrials,
                           family = binomial())
    # model
    fittedModel <- glmer(cbind(observedResponse1, observedResponse0)  ~ 
                           Environment1 + (1|group),
                         data = testData, family = binomial()) 
    #results
    out <- list()
    
    # pearson residual
    out$Pear2.p.val <- testDispersion(fittedModel, plot = F, 
                                      alternative = "two.sided",
                                     type="PearsonChisq")$p.value
    out$PearG.p.val <- testDispersion(fittedModel, plot = F, 
                                      alternative = "great",
                                      type="PearsonChisq")$p.value
    out$Pear.stat <- testDispersion(fittedModel, plot = F, 
                                    type="PearsonChisq")$statistic
    return(unlist(out))
  }
  
  out <- runBenchmarks(calculateStatistics, controlValues = overdispersion,
                       nRep = nRep, parallel = T, exportGlobal = T)
  out.bin[[length(out.bin) + 1]] <- out
  names(out.bin)[length(out.bin)] <- m
  # saving sim results
  save(out.bin, file=here("data","4_glmmBin_pearsonChisq.Rdata"))
}


#########################
##### Poisson GLMM #####
#########################

out.pois <- list()
for(m in ngroups){

  calculateStatistics <- function(control = 0){
    # data
    testData <- createData(overdispersion = control,
                           sampleSize = sampleSize,
                           intercept = intercept,
                           numGroups = m,
                           randomEffectVariance = 1,
                           family = poisson())
    # model
    fittedModel <- glmer(observedResponse  ~ 
                           Environment1 + (1|group),
                         data = testData, family = poisson()) 
    #results
    out <- list()
    
    # pearson residual
    out$Pear2.p.val <- testDispersion(fittedModel, plot = F, 
                                      alternative = "two.sided",
                                      type="PearsonChisq")$p.value
    out$PearG.p.val <- testDispersion(fittedModel, plot = F, 
                                      alternative = "great",
                                      type="PearsonChisq")$p.value
    out$Pear.stat <- testDispersion(fittedModel, plot = F, 
                                    type="PearsonChisq")$statistic
    return(unlist(out))
  }
  
  out <- runBenchmarks(calculateStatistics, controlValues = overdispersion,
                       nRep = nRep, parallel = T, exportGlobal = T)
  out.pois[[length(out.pois) + 1]] <- out
  names(out.pois)[length(out.pois)] <-m
  # saving sim results
  save(out.pois, file=here("data","4_glmmPois_pearsonChisq.Rdata"))
}




