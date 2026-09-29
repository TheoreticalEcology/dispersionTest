### Dispersion Tests Project
## Power of dispersion tests under UNDERDISPERSION (GLMs)
# Sep 26

library(DHARMa)
library(tidyverse)
library(here)

# createData_ud(): DHARMa::createData() extended with the arguments
#   underdispersion = u in [0,1)  -> target Var(y|eta) ~ (1 - u) * nominal variance
#   underdispersionModel          -> "CMP" (Poisson) or "CMB" (binomial)
# The conditional mean is kept at family$linkinv(eta) in both cases.
#   Poisson : Conway-Maxwell-Poisson (mean parametrisation), nu = 1/(1-u)
#   Binomial: Conway-Maxwell-binomial (mean-matched),          nu = 1/(1-u)
source(here("functions_others", "createData_underdispersion.R"))

# Note: all tests are run two-sided (DHARMa default, as in 3_glm_power.R).
# One-sided p-values (alternative = "less") are stored as well (*.p.less),
# at no extra simulation cost.


####################
##### Binomial #####
####################

# 1) Simulating 1000 binomial prop datasets with different sample sizes, intercepts
#    and underdispersion (Conway-Maxwell-binomial). Fixing the number of trials to 10
# 2) fitting them to correct GLM models
# 3) calculating power for the dispersion tests
#     - Pearson-chisq,
#     - DHARMa default (simulated residuals)
#     - DHARMa refit (boostrapped Pearson residuals)


# varying parameters
underdispersion <- seq(0, 0.9, 0.10)   # nu = 1/(1-u): 1 (= binomial) ... 10
sampleSize = c(10,50,100,200,500,1000)
intercept <- c(-1.5,0,1.5) # no extreme -3 or 3


out.bin <- list()
for (k in sampleSize){
  for (i in intercept){

    # function to varying underdispersion
    calculateStatistics <- function(control = 0){
      # data
      testData <- createData_ud(underdispersion = control,
                                underdispersionModel = "CMB",
                                sampleSize = k,
                                intercept = i,
                                numGroups = 10,
                                randomEffectVariance = 0,
                                binomialTrials = 10,
                                family = binomial())
      # model
      fittedModel <- glm(cbind(observedResponse1, observedResponse0)  ~
                          Environment1, data = testData, family = binomial())
      #results
      out <- list()

      # pearson residual
      t2 <- testDispersion(fittedModel, plot = F, type = "PearsonChisq")
      out$Pear.p.val  <- t2$p.value
      out$Pear.stat   <- t2$statistic
      out$Pear.p.less <- testDispersion(fittedModel, plot = F, type = "PearsonChisq",
                                        alternative = "less")$p.value

      # DHARMa default residuals
      res <- simulateResiduals(fittedModel)
      t2 <- testDispersion(res, type = "DHARMa", plot = F)
      out$DHA.p.val  <- t2$p.value
      out$DHA.stat   <- t2$statistic
      out$DHA.p.less <- testDispersion(res, type = "DHARMa", plot = F,
                                       alternative = "less")$p.value

      # DHARMa refit residuals -> bootstrapped Pearson
      res <- simulateResiduals(fittedModel, refit = T)
      t2 <- testDispersion(res, type = "DHARMa", plot = F)
      out$Ref.p.val  <- t2$p.value
      out$Ref.stat   <- t2$statistic
      out$Ref.p.less <- testDispersion(res, type = "DHARMa", plot = F,
                                       alternative = "less")$p.value

      return(unlist(out))
    }

    out <- runBenchmarks(calculateStatistics, controlValues = underdispersion,
                         nRep = 1000, parallel = T, exportGlobal = T)
    out.bin[[length(out.bin) + 1]] <- out
  }
}

names(out.bin) <- as.vector(unite(expand.grid(intercept,sampleSize), "sim"))$sim

# saving sim results
save(out.bin, sampleSize, intercept, underdispersion,
     file = here("data", "3_glmBin_power_underdisp.Rdata"))



###################
##### Poisson #####
###################
# 1) Simulating 1000 Poisson datasets with different sample sizes, intercepts and
#    underdispersion (Conway-Maxwell-Poisson, mean parametrisation).
# 2) fitting them to correct GLM models
# 3) calculating power for the dispersion tests
#     - Pearson-chisq,
#     - DHARMa default (simulated residuals)
#     - DHARMa refit (boostrapped Pearson residuals)


# varying parameters
underdispersion <- seq(0, 0.9, 0.10)   # nu = 1/(1-u): 1 (= Poisson) ... 10
sampleSize = c(10,50,100,200,500,1000)
intercept <- c(-1.5,0,1.5,3) # no extreme -3


out.pois <- list()
for (k in sampleSize){
  for (i in intercept){

    # function to varying underdispersion
    calculateStatistics <- function(control = 0){
      # data
      testData <- createData_ud(underdispersion = control,
                                underdispersionModel = "CMP",
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
      t2 <- testDispersion(fittedModel, plot = F, type = "PearsonChisq")
      out$Pear.p.val  <- t2$p.value
      out$Pear.stat   <- t2$statistic
      out$Pear.p.less <- testDispersion(fittedModel, plot = F, type = "PearsonChisq",
                                        alternative = "less")$p.value

      # DHARMa default residuals
      res <- simulateResiduals(fittedModel)
      t2 <- testDispersion(res, type = "DHARMa", plot = F)
      out$DHA.p.val  <- t2$p.value
      out$DHA.stat   <- t2$statistic
      out$DHA.p.less <- testDispersion(res, type = "DHARMa", plot = F,
                                       alternative = "less")$p.value

      # DHARMa refit residuals -> bootstrapped Pearson
      res <- simulateResiduals(fittedModel, refit = T)
      t2 <- testDispersion(res, type = "DHARMa", plot = F)
      out$Ref.p.val  <- t2$p.value
      out$Ref.stat   <- t2$statistic
      out$Ref.p.less <- testDispersion(res, type = "DHARMa", plot = F,
                                       alternative = "less")$p.value

      return(unlist(out))
    }

    out <- runBenchmarks(calculateStatistics, controlValues = underdispersion,
                         nRep = 1000, parallel = T, exportGlobal = T)
    out.pois[[length(out.pois) + 1]] <- out
  }
}

names(out.pois) <- as.vector(unite(expand.grid(intercept,sampleSize), "sim"))$sim


# saving sim results
save(out.pois, sampleSize, intercept, underdispersion,
     file = here("data", "3_glmPois_power_underdisp.Rdata"))
