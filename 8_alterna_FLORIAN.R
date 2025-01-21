### Dispersion Tests Project
## Melina Leite
# jan 25

library(DHARMa)
library(here)


# testing Florian's idea to check error in the SD calculation for the simulated version of Pearson Residuals when we have low simulations

# 1) Simulating data without zero in the SD of the simulations, 
# 2) Changing only the number of simulations

# function for the "Pearson simulated residuals"
pear.spread <- function(x){
  (x - simulationOutput$fittedPredictedResponse)/simSD
}

# getP from DHARMa
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


# params
sampleSize = 100

calculateStatistics <- function(control = 5){

  testData <- createData(overdispersion = 0,
                         sampleSize = sampleSize,
                         intercept = 3,
                         numGroups = 10,
                         family = poisson())
  # model
  fittedModel <-glm(observedResponse ~ 
                      Environment1, data = testData, family = poisson()) 
  #results
  out <- list()
  
  simulationOutput <- simulateResiduals(fittedModel, n = control)
  
  # Observation-wise SD
  simSD <-  apply(simulationOutput$simulatedResponse,1,sd)
  
  out$zeroSD <- min(simSD)==0

  # pearson sim residual
  obs.spread <- pear.spread(simulationOutput$observedResponse)
  sum.obs.spread2 <- sum(obs.spread^2)
  sim.spread <- apply(simulationOutput$simulatedResponse, 2, pear.spread)
  sim.spread2 <-sim.spread^2 
  sum.sim.spread2 <- apply(sim.spread2,2,sum)
  
  out$pear.spread.stat <- sum.obs.spread2/mean(sum.sim.spread2)
 
  out$pval <- getP(simulated = sum.sim.spread2, observed = sum.obs.spread2,"two.sided")
  
  return(out)
   
}
calculateStatistics(300)
