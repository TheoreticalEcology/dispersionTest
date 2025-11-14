### Dispersion Tests Project
## 
# Fev 25

library(DHARMa)
library(here)

# testing the diff alternatives for DHARMa disp test

sampleSize <- c(100,1000)
intercept <- c(-1.5,0,1.5)
nSim = c(5,10,50,250,1000)
overdispersion = seq(0,1,0.1)
nRep = 1000


# functions for the alternative
getP <- function(simulated, observed, alternative, plot = FALSE, ...){
  
  if(alternative == "greater") p = mean(simulated >= observed)
  if(alternative == "less") p = mean(simulated <= observed)
  if(alternative == "two.sided") p = min(min(mean(simulated <= observed), mean(simulated >= observed) ) * 2,1)
  
  if(plot == T){
    hist(simulated, xlim = range(simulated, observed), col = "lightgrey", main = "Distribution of test statistic \n grey = simulated, red = observed", ...)
    abline(v = mean(simulated), col = 1, lwd = 2)
    abline(v = observed, col = "red", lwd = 2)
  }
  
  return(p)
}

getApproximatePearson <- function(simulationOutput, alternative = c("two.sided","greater", "less"), plot = T, ...){
  
  expectedSD = apply(simulationOutput$simulatedResponse, MARGIN = 1, sd)
  noSimVar = expectedSD == 0
  expectedSD[noSimVar] = sd(c(rep(0, simulationOutput$nSim - 1),1))
  
  ntrials <- simulationOutput[["fittedModel"]][["prior.weights"]]
  
  spread <- function(x) sum(((x - simulationOutput$fittedPredictedResponse*ntrials) / expectedSD )^2) 
  out = testGeneric(simulationOutput, summary = spread, alternative = alternative, methodName = "DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated", plot = plot, ...)
  names(out$statistic) = "dispersion"
  out$noSimVar = noSimVar
  out$expectedSD = expectedSD
  out$resPA = simulationOutput$fittedResiduals / expectedSD # doesn't make sense for binomial
} 


# poisson

# out.pois<- list()
# #load(here("data","6_approxPear_pois.Rdata"))
# 
# for(s in nSim) {
#   for (k in sampleSize){
#     for (i in intercept){
# 
#     # function to varying sampleSize
#     calculateStatistics <- function(control = 0){
#       # data
#       testData <- createData(overdispersion = control,
#                              sampleSize = k,
#                              intercept = i,
#                              numGroups = 10,
#                              randomEffectVariance = 0,
#                              family = poisson())
#       # model
#       fittedModel <-glm(observedResponse ~
#                           Environment1, data = testData, family = poisson())
#       #results
#       out <- list()
# 
#       simulationOutput <- simulateResiduals(fittedModel, n = s)
# 
#       #Pearson
#       out$Pear.p <- testDispersion(fittedModel, plot = F,
#                                        type="PearsonChisq")$p.value
#       out$Pear.stat <- testDispersion(fittedModel, plot = F,
#                                       type="PearsonChisq")$statistic
#       #DHARMa default
#       out$DHA.p <- testDispersion(simulationOutput, type = "DHARMa",
#                                   plot = F)$p.value
#       out$DHA.stat  <- testDispersion(simulationOutput, type = "DHARMa",
#                                       plot = F)$statistic
# 
#       # Alternative
#       alterna <- getApproximatePearson(simulationOutput, "two.sided", plot = F)
#       out$Alt.p <- alterna$p.value
#       out$Alt.stat <- alterna$statistic
#       out$prop.zero <- sum(alterna$noSimVar)/length(alterna$noSimVar)
# 
#       resP = residuals(fittedModel, type = "pearson")
#       m2 <- summary(lm(alterna$resPA - resP ~ resP))
#       out$resPA.inter = m2$coefficients[1,1]
#       out$resPA.interP = m2$coefficients[1,4]
#       out$resPA.slope = m2$coefficients[2,1]
#       out$resPA.slopeP = m2$coefficients[2,4]
# 
#       return(unlist(out))
#     }
# 
#     out <- runBenchmarks(calculateStatistics, controlValues = overdispersion,
#                          nRep = nRep, parallel = T, exportGlobal = T)
#     out.pois[[length(out.pois) + 1]] <- out
#     names(out.pois)[length(out.pois)] <- paste(s, k, i, sep="_")
#     # saving sim results
#     save(out.pois, file=here("data","6_approxPear_pois.Rdata"))
#   }
# }
# }

# binomial

out.bin <- list()
#load(here("data","6_approxPear_bin.Rdata"))

for(s in nSim) {
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
        
        simulationOutput <- simulateResiduals(fittedModel, n = s)
        
        #Pearson
        out$Pear.p <- testDispersion(fittedModel, plot = F, 
                                     type="PearsonChisq")$p.value
        out$Pear.stat <- testDispersion(fittedModel, plot = F, 
                                        type="PearsonChisq")$statistic
        #DHARMa default
        out$DHA.p <- testDispersion(simulationOutput, type = "DHARMa",
                                    plot = F)$p.value
        out$DHA.stat  <- testDispersion(simulationOutput, type = "DHARMa",
                                        plot = F)$statistic
        
        # Alternative
        alterna <- getApproximatePearson(simulationOutput, "two.sided", plot = F)
        out$Alt.p <- alterna$p.value
        out$Alt.stat <- alterna$statistic
        out$prop.zero <- sum(alterna$noSimVar)/length(alterna$noSimVar)
        
        resP = residuals(fittedModel, type = "pearson")
        m2 <- summary(lm(alterna$resPA - resP ~ resP))
        out$resPA.inter = m2$coefficients[1,1]
        out$resPA.interP = m2$coefficients[1,4]
        out$resPA.slope = m2$coefficients[2,1]
        out$resPA.slopeP = m2$coefficients[2,4]
        
        return(unlist(out))
      }
      
      out <- runBenchmarks(calculateStatistics, controlValues = overdispersion,
                           nRep = nRep, parallel = T, exportGlobal = T)
      out.bin[[length(out.bin) + 1]] <- out
      names(out.bin)[length(out.bin)] <- paste(s, k, i, sep="_")
      # saving sim results
      save(out.bin, file=here("data","6_approxPear_bin.Rdata"))
    }
  }
}


