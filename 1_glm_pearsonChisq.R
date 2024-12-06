### Dispersion Tests Project
## Melina Leite
# Dec 24

library(DHARMa)
library(tidyverse)
library(here)

####################
##### Binomial #####
####################

## KS GOF test to compare Pearson Stats and Chi-squared distribution

# 1) Simulating 1000 Binomial proportion datasets with different sample sizes and intercepts. Fixing number of trials in 10.
# 2) fitting them to correct GLM model
# 3) calculating pearson statistics of the pearson residuals (sum(resˆ2)) 
# 4) comparing distribution of these statistics with the Chi-squared distribution with the same DF.
# 5) Repeating these steps 100 times to get the proportion of significant results per combination of parameters.


# sampleSizes and intercepts
sampleSize = c(10,20,50,100,200,500,1000,10000)
intercept <- c(-3,-1.5,0,1.5,3)

# KS TEST
ks.p <- function(x, df) ks.test(x,"pchisq", df=df)$p.value

final.res.bin <- list()
final.sims.bin <- list()

for (k in 1:100) { # MANY SIMULATIONS TO HAVE A PROP OF SIG RESULTS
  set.seed(k)  
  result <- list()
  sims <- list()
  
  
  for (i in 1:length(intercept)){
    
    calculateStatistics <- function(control = 10){
      testData <- DHARMa::createData(sampleSize = control, 
                                     intercept = intercept[i],
                                     numGroups = 10,
                                     randomEffectVariance = 0,
                                     binomialTrials = 10,
                                     family = binomial())
      
      fittedModel <- stats::glm(cbind(observedResponse1,observedResponse0) ~
                                  Environment1, data = testData, 
                                family = binomial()) 
      
      pearson <- residuals(fittedModel, "pearson")
      
      out <- list()
      out$Pear.stat <- sum(residuals(fittedModel, "pearson")^2)
      out$rdf <- df.residual(fittedModel)
      return(unlist(out))
    }
    
    out <- runBenchmarks(calculateStatistics, controlValues = sampleSize,
                         nRep=1000, parallel = T, exportGlobal = TRUE)
    #simulations
    stats <- as.data.frame(apply(as.data.frame(out$simulations), 2, unlist))
    stats$intercept <- intercept[i]
    sims <- rbind(stats,sims)
    
    #results
    res <- stats %>% group_by(controlValues) %>%
      summarise(ks.p = ks.p(Pear.stat,rdf))
    res$intercept <- intercept[i]
    result <- rbind(result,res)
  }
  final.res.bin[[k]] <- result
  final.sims.bin[[k]] <- sims
}

# saving results
save(final.res.bin, final.sims.bin, file = here("data", "1_glmBin_pearsonChisq.Rdata"))

rm(ls())


###################
##### Poisson #####
###################

## KS GOF test to compare Pearson Stats and Chi-squared distribution

# 1) Simulating 1000 Binomial proportion datasets with different sample sizes and intercepts. Fixing number of trials in 10.
# 2) fitting them to correct GLM model
# 3) calculating pearson statistics of the pearson residuals (sum(resˆ2)) 
# 4) comparing distribution of these statistics with the Chi-squared distribution with the same DF.
# 5) Repeating these steps 100 times to get the proportion of significant results per combination of parameters.


# sampleSizes and intercepts
sampleSize = c(10,20,50,100,200,500,1000,10000)
intercept <- c(-3,-1.5,0,1.5,3)

# KS TEST
ks.p <- function(x, df) ks.test(x,"pchisq", df=df)$p.value

final.res.pois <- list()
final.sims.pois <- list()

for (k in 1:100) { # MANY SIMULATIONS TO HAVE A PROP OF SIG RESULTS
  set.seed(k)  
  result <- list()
  sims <- list()
  
  
  for (i in 1:length(intercept)){
    
    calculateStatistics <- function(control = 10){
      testData <- DHARMa::createData(sampleSize = control, 
                                     intercept = intercept[i],
                                     numGroups = 10,
                                     randomEffectVariance = 0,
                                     family = poisson())
      
      fittedModel <- stats::glm(observedResponse ~
                                  Environment1, data = testData, 
                                family = poisson()) 
      
      pearson <- residuals(fittedModel, "pearson")
      
      out <- list()
      out$Pear.stat <- sum(residuals(fittedModel, "pearson")^2)
      out$rdf <- df.residual(fittedModel)
      return(unlist(out))
    }
    
    out <- runBenchmarks(calculateStatistics, controlValues = sampleSize,
                         nRep=1000, parallel = T, exportGlobal = TRUE)
    #simulations
    stats <- as.data.frame(apply(as.data.frame(out$simulations), 2, unlist))
    stats$intercept <- intercept[i]
    sims <- rbind(stats,sims)
    
    #results
    res <- stats %>% group_by(controlValues) %>%
      summarise(ks.p = ks.p(Pear.stat,rdf))
    res$intercept <- intercept[i]
    result <- rbind(result,res)
  }
  final.res.pois[[k]] <- result
  final.sims.pois[[k]] <- sims
}

# saving results
save(final.res.pois, final.sims.pois, file = here("data", "1_glmPois_pearsonChisq.Rdata"))
