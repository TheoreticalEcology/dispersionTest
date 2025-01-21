### Dispersion Tests Project
## Melina Leite
# jan 25

library(DHARMa)
library(here)


# testing Florian's idea to check error in the SD calculation for the simulated version of Pearson Residuals when we have low simulations

# 1) Simulating data without zero in the SD of the simulations, 
# 2) Changing only the number of simulations

calculateStatistics <- function(control = 5){

  testData <- createData(overdispersion = 0,
                         sampleSize = 100,
                         intercept = 3,
                         numGroups = 10,
                         family = poisson())
  # model
  fittedModel <-glm(observedResponse ~ 
                      Environment1, data = testData, family = poisson()) 
  #results
  simulationOutput <- simulateResiduals(fittedModel, n = control)
  
  # Observation-wise SD
  simulationSD.obs <-  apply(simulationOutput$simulatedResponse,1,sd)
  min(simulationSD.obs)==0

  }