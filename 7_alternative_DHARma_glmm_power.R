### Dispersion Tests Project
## Melina Leite
# feb 25

library(DHARMa)
library(here)


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
  
  spread <- function(x) sum(((x - simulationOutput$fittedPredictedResponse) / expectedSD )^2) 
  out = testGeneric(simulationOutput, summary = spread, alternative = alternative, methodName = "DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated", plot = plot, ...)
  names(out$statistic) = "dispersion"
  out$noSimVar = noSimVar
  out$expectedSD = expectedSD
  out$resPA = simulationOutput$fittedResiduals / expectedSD
  return(out)
} 



# testing the alternative for DHARMa disp test

# varying parameters
overdispersion <- seq(0,1,0.10)
intercept <- c(-1.5,0,1.5)
## sample size depend on ngroups (see it within the loop)
ngroups <- 100
nRep = 100


out.binA <- list()
#load(here("data", "7_glmmBin_alternative_power.Rdata")) # loading out.bin

for(m in ngroups){
  
  if(m == 10) sampleSize <- c(20,50,100,200,500,1000) # 10,000 excluded
  if(m == 50) sampleSize <- c(100,200,500,1000) # 10,000 excluded
  if(m == 100) sampleSize <- c(200,500,1000) # 10,000 excluded
  
  for (k in sampleSize){
    
    for (i in intercept){
      
      # function to varying sampleSize
      calculateStatistics <- function(control = 0){
        # data
        testData <- createData(overdispersion = control,
                               sampleSize = k,
                               intercept = i,
                               numGroups = m,
                               randomEffectVariance = 1,
                               binomialTrials = 10,
                               family = binomial())
        # model
        fittedModel <- glmer(cbind(observedResponse1, observedResponse0)  ~
                               Environment1 + (1|group),
                             data = testData, family = binomial())
        #results
        out <- list()

        # Alternative DHARMa test
        alterna <- getApproximatePearson(simulationOutput, "two.sided", plot = F)
        out$Alt.p <- alterna$p.value
        out$Alt.stat <- alterna$statistic
        out$prop.zero <- sum(alterna$noSimVar)/length(alterna$noSimVar)
        
        
        return(unlist(out))
      }
      
      out <- runBenchmarks(calculateStatistics, controlValues = overdispersion,
                           nRep = nRep, parallel = T, exportGlobal = T)
      out.binA[[length(out.binA) + 1]] <- out
      names(out.binA)[length(out.binA)] <- paste(m, k, i, sep="_")
      # saving sim results
      save(out.binA, file=here("data","7_glmmBin_alternative_power.Rdata"))
    }
  }
}

######################-#
##### Poisson GLMM #####
#####################-##

out.poisA <- list()
#load(here("data", "7_glmmPois_alternative_power.Rdata")) # loading out.pois

for(m in ngroups){
  
  if(m == 10) sampleSize <- c(50,100,200,500,1000) # 20, 10,000 excluded
  if(m == 50) sampleSize <- c(100,200,500,1000) # 10,000 excluded
  if(m == 100) sampleSize <- c(200,500,1000,10000)
  
  for (k in sampleSize){
    
    for (i in intercept){
      
      # function to varying sampleSize
      calculateStatistics <- function(control = 0){
        # data
        testData <- createData(overdispersion = control,
                               sampleSize = k,
                               intercept = i,
                               numGroups = m,
                               randomEffectVariance = 1,
                               family = poisson())
        # model
        fittedModel <- glmer(observedResponse  ~ 
                               Environment1 + (1|group),
                             data = testData, family = poisson()) 
        #results
        out <- list()
        
        # Alternative DHARMa test
        alterna <- getApproximatePearson(simulationOutput, "two.sided", plot = F)
        out$Alt.p <- alterna$p.value
        out$Alt.stat <- alterna$statistic
        out$prop.zero <- sum(alterna$noSimVar)/length(alterna$noSimVar)
        
        return(unlist(out))
      }
      
      out <- runBenchmarks(calculateStatistics, controlValues = overdispersion,
                           nRep=nRep, parallel = T, exportGlobal = T)
      out.poisA[[length(out.poisA) + 1]] <- out
      names(out.poisA)[length(out.poisA)] <- paste(m, k, i, sep="_")
      # saving sim results
      save(out.poisA, file=here("data","7_glmmPois_alternative_power.Rdata"))
    }
  }
}