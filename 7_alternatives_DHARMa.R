### Dispersion Tests Project
## Melina Leite
# jan 25

library(DHARMa)
library(here)

# testing the diff alternatives for DHARMa disp test

sampleSize <- c(100,1000)
intercept <- c(-3,-1.5,0)
nSim = c(250,1000)
nRep = 10000


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

spread2 <- function(x, simulationOutput, simSD.obs){
  val <- (x - simulationOutput$fittedPredictedResponse)/simSD.obs
  (sum(val^2)/(length(x)-1))
}

dispersion <- function(simulationOutput, simSD.obs){
  observed = spread2(simulationOutput$observedResponse,simulationOutput, simSD.obs)
  simulated = apply(simulationOutput$simulatedResponse, 2, spread2, 
                    simulationOutput = simulationOutput, 
                    simSD.obs = simSD.obs)
  dispersion = observed/mean(simulated)
  pval = getP(simulated, observed, "two.sided")
  return(list(statistic = dispersion, p.value = pval))
}


# poisson

out.pois<- list()
#load(here("data","6_DHARMa_dispersion_pois.Rdata"))


for (k in sampleSize){
  for (i in nSim){
    
    # function to varying sampleSize
    calculateStatistics <- function(control = 0){
      # data
      testData <- createData(overdispersion = 0,
                             sampleSize = k,
                             intercept = control,
                             numGroups = 10,
                             family = poisson())
      # model
      fittedModel <-glm(observedResponse ~ 
                          Environment1, data = testData, family = poisson()) 
      #results
      out <- list()
      
      simulationOutput <- simulateResiduals(fittedModel, n = i)
      
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
      
      # Alternatives
      # Observation-wise SD
      simulationSD.obs <-  apply(simulationOutput$simulatedResponse,1,sd)
      
      # if there is a simOBS with zero sd
      out$zeroSD <- min(simulationSD.obs)==0
      out$prop.zeroSD <- sum(simulationSD.obs==0)/simulationOutput$nObs
      
      # Alternative 1: use te minimum SD
      simSD.obs1 <- simulationSD.obs
      if(min(simulationSD.obs)==0){
        newSD <- min(simulationSD.obs[simulationSD.obs != 0])
        simSD.obs1[simSD.obs1==0] <- newSD
      }
      out$A1.p <- dispersion(simulationOutput, simSD.obs1)$p.value
      out$A1.stat <- dispersion(simulationOutput,simSD.obs1)$statistic
      
      # Alternative 2: min SD divided by the max N unique vals
      n.uniq <- function(x) length(unique(x))
      nuniq <- apply(simulationOutput$simulatedResponse,1,n.uniq)
      
      simSD.obs2 <- simulationSD.obs
      if(min(simulationSD.obs)==0){
        newSD <- min(simulationSD.obs[simulationSD.obs != 0])/max(nuniq)
        simSD.obs2[simSD.obs2==0] <- newSD
      }
      out$A2.p <- dispersion(simulationOutput,simSD.obs2)$p.value
      out$A2.stat <- dispersion(simulationOutput,simSD.obs2)$statistic
      
      # Alternative 3: min SD divided by the nSim
      simSD.obs3 <- simulationSD.obs
      if(min(simulationSD.obs)==0){
        newSD <- min(simulationSD.obs[simulationSD.obs != 0])/
          simulationOutput$nSim
        simSD.obs3[simSD.obs3==0] <- newSD
      }
      out$A3.p <- dispersion(simulationOutput,simSD.obs3)$p.value
      out$A3.stat <- dispersion(simulationOutput,simSD.obs3)$statistic

      return(unlist(out))
    }
    
    out <- runBenchmarks(calculateStatistics, controlValues = intercept,
                         nRep = nRep, parallel = T, exportGlobal = T)
    out.pois[[length(out.pois) + 1]] <- out
    names(out.pois)[length(out.pois)] <- paste(k, i, sep="_")
    # saving sim results
    save(out.pois, file=here("data","7_alternatives_pois.Rdata"))
  }
}


# binomial

# out.bin <- list()
# #load(here("data","6_DHARMa_dispersion_bin.Rdata"))
# 
# for (k in sampleSize){
#   for (i in nSim){
#     
#     # function to varying sampleSize
#     calculateStatistics <- function(control = 0){
#       # data
#       testData <- createData(overdispersion = 0,
#                              sampleSize = k,
#                              intercept = control,
#                              numGroups = 10,
#                              randomEffectVariance = 0,
#                              binomialTrials = 10,
#                              family = binomial())
#       # model
#       fittedModel <-glm(cbind(observedResponse1, observedResponse0)  ~ 
#                           Environment1, data = testData, family = binomial()) 
#       #results
#       out <- list()
#       
#       res <- simulateResiduals(fittedModel, n = i)
#       zeros <- sum( apply(res$simulatedResponse,1,sd)==0)
#       out$zero = zeros
#       out$prop.zero = zeros/res$nObs
#       return(unlist(out))
#     }
#     
#     out <- runBenchmarks(calculateStatistics, controlValues = intercept,
#                          nRep = nRep, parallel = T, exportGlobal = T)
#     out.bin[[length(out.bin) + 1]] <- out
#     names(out.bin)[length(out.bin)] <- paste(k, i, sep="_")
#     # saving sim results
#     save(out.bin, file=here("data","7_alternatives_bin.Rdata"))
#   }
# }






























