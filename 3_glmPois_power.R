### Dispersion Tests Project
## Melina Leite
# Dec 24

library(DHARMa)
library(here)


#####################################
##### Instructions & Simulatipms #####
#####################################

# 1) Simulating 1000 Poisson datasets with different sample sizes, intercepts and overdispertion.
#       - sampleSize: c(10,50,100,500)
#       - intercept:  c(-3,-1,0,2,4)
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
for (k in sampleSize){
  for (i in interept){
  
  # function to varying overdispersion
  calculateStatistics <- function(control = 0){
    # data
    testData <- DHARMa::createData(overdispersion = control,
                                   sampleSize = k,
                                   intercept = i
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
    # DHARMa default residuals
    res <- simulateResiduals(fittedModel)
    out$DHA.p.val<- testDispersion(res, type = "DHARMa",plot = F)$p.value
    
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

#names(out.out) <- sampleSize

# saving sim results
save(out.out,sampleSize,intercept,overdispersion, file=here("data", 
                                             "3_glmPois_power.Rdata"))



# prep data

#load(here("data", "4_glmPois_power.Rdata"))

# props <- list()
# for (i in 1:length(out.out)){
#   suma <- out.out[[i]]$simulations
#   suma$sampleSize <- sampleSize[i]
# }









