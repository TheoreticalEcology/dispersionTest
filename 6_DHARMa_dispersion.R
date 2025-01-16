### Dispersion Tests Project
## Melina Leite
# jan 25

library(DHARMa)
library(here)

# Evaluating the frequency of zero variance/sd in simulated observations for:

sampleSize <- c(100,1000)
intercept <- c(-3,-1.5,0)
nSim = c(1000)
nRep = 10000

# binomial

#out.bin <- list()
load(here("data","6_DHARMa_dispersion_bin.Rdata"))

for (k in sampleSize){
  for (i in nSim){
    
    # function to varying sampleSize
    calculateStatistics <- function(control = 0){
      # data
      testData <- createData(overdispersion = 0,
                             sampleSize = k,
                             intercept = control,
                             numGroups = 10,
                             randomEffectVariance = 0,
                             binomialTrials = 10,
                             family = binomial())
      # model
      fittedModel <-glm(cbind(observedResponse1, observedResponse0)  ~ 
                          Environment1, data = testData, family = binomial()) 
      #results
      out <- list()
      
      res <- simulateResiduals(fittedModel, n = i)
      zeros <- sum( apply(res$simulatedResponse,1,sd)==0)
      out$zero = zeros
      out$prop.zero = zeros/res$nObs
      return(unlist(out))
    }
    
    out <- runBenchmarks(calculateStatistics, controlValues = intercept,
                         nRep = nRep, parallel = T, exportGlobal = T)
    out.bin[[length(out.bin) + 1]] <- out
    names(out.bin)[length(out.bin)] <- paste(k, i, sep="_")
    # saving sim results
    save(out.bin, file=here("data","6_DHARMa_dispersion_bin.Rdata"))
  }
}

# poisson

#out.pois<- list()
load(here("data","6_DHARMa_dispersion_pois.Rdata"))


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
      
      res <- simulateResiduals(fittedModel, n = i)
      zeros <- sum( apply(res$simulatedResponse,1,sd)==0)
      out$zero = zeros
      out$prop.zero = zeros/res$nObs
      return(unlist(out))
    }
    
    out <- runBenchmarks(calculateStatistics, controlValues = intercept,
                         nRep = nRep, parallel = T, exportGlobal = T)
    out.pois[[length(out.pois) + 1]] <- out
    names(out.pois)[length(out.pois)] <- paste(k, i, sep="_")
    # saving sim results
    save(out.pois, file=here("data","6_DHARMa_dispersion_pois.Rdata"))
  }
}



# results:
# library(tidyverse)
# 
# # seeing data:
# load(here("data","6_DHARMa_dispersion_bin.Rdata"))
# load(here("data","6_DHARMa_dispersion_pois.Rdata"))
# 
# simbin < map_dfr(out.bin, "simulations", .id="ngroups")  %>%
#   separate(ngroups, c("sampleSize", "nSim")) %>%
#   rename("intercept" = "controlValues")
# simpois <- map_dfr(out.pois, "simulations", .id="ngroups")  %>%
#   separate(ngroups, c("sampleSize", "nSim")) %>%
#   rename("intercept" = "controlValues")
# 
# res <- bind_rows(list(Poisson = simpois, Binomial = simbin), .id="model")
# 
# 
# res %>%
#   mutate(nSim = fct_relevel(nSim, "50","100","250", "1000") ) %>%
#   ggplot(aes(x = nSim, y = prop.zero, col = sampleSize))+
#   geom_boxplot() +
#   #geom_violin() +
#   facet_grid(model ~ intercept, scales = "free") +
#   scale_y_sqrt()
# 



























