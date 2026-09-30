### Dispersion Tests Project
## Power of dispersion tests under alternative OVERDISPERSION mechanisms (GLMs):
## negative binomial (Poisson) and beta-binomial (binomial) - Rev. 2, Comment 2
# Sep 26

library(DHARMa)
library(tidyverse)
library(here)

# createData_ud(): DHARMa::createData() extended with the arguments
#   phi >= 1             -> Var(y|eta) = phi * nominal variance (exact)
#   overdispersionModel  -> "NB" (Poisson) or "betabinomial" (binomial)
# The conditional mean is kept at family$linkinv(eta) in both cases, unlike the
# default DHARMa mechanism used in 3_glm_power.R (Gaussian noise on the linear
# predictor, which also shifts the marginal mean, e.g. exp(eta + sd^2/2)).
#   Poisson : negative binomial, size = mu/(phi - 1)       -> Var = phi * mu
#   Binomial: beta-binomial, rho = (phi - 1)/(Ntrials - 1)  -> Var = phi * N p (1 - p)
# phi is therefore the same for every observation, intercept and sample size,
# i.e. a common effect-size scale across simulation settings (Rev. 2, Comment 13).
source(here("functions_others", "createData_underdispersion.R"))

# Note: all tests are run two-sided (DHARMa default, as in 3_glm_power.R).
# One-sided p-values (alternative = "greater") are stored as well (*.p.greater),
# at no extra simulation cost.

# Dispersion factor phi: 11 levels (as the 11 overdispersion levels in
# 3_glm_power.R), roughly log-spaced to have resolution close to 1 (large n)
# and a range up to 5 (small n). phi = 1 is exactly the Poisson / binomial null.
# Beta-binomial with Ntrials = 10: phi = 5 <-> rho = 0.44 (phi must be < Ntrials).
phi <- c(1, 1.1, 1.2, 1.35, 1.5, 1.75, 2, 2.5, 3, 4, 5)


####################
##### Binomial #####
####################

# 1) Simulating 1000 binomial prop datasets with different sample sizes, intercepts
#    and overdispersion (beta-binomial). Fixing the number of trials to 10
# 2) fitting them to correct (binomial) GLM models
# 3) calculating power for the dispersion tests
#     - Pearson-chisq,
#     - DHARMa default (simulated residuals)
#     - DHARMa refit (boostrapped Pearson residuals)


# varying parameters (same grid as 3_glm_power_underdisp.R)
sampleSize = c(10,50,100,200,500,1000)
intercept <- c(-1.5,0,1.5)


out.bin <- list()
for (k in sampleSize){
  for (i in intercept){

    # function to varying overdispersion (phi)
    calculateStatistics <- function(control = 1){
      # data
      testData <- createData_ud(phi = control,
                                overdispersionModel = "betabinomial",
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
      out$Pear.p.val     <- t2$p.value
      out$Pear.stat      <- t2$statistic
      out$Pear.p.greater <- testDispersion(fittedModel, plot = F, type = "PearsonChisq",
                                           alternative = "greater")$p.value

      # DHARMa default residuals
      res <- simulateResiduals(fittedModel)
      t2 <- testDispersion(res, type = "DHARMa", plot = F)
      out$DHA.p.val     <- t2$p.value
      out$DHA.stat      <- t2$statistic
      out$DHA.p.greater <- testDispersion(res, type = "DHARMa", plot = F,
                                          alternative = "greater")$p.value

      # DHARMa refit residuals -> bootstrapped Pearson
      res <- simulateResiduals(fittedModel, refit = T)
      t2 <- testDispersion(res, type = "DHARMa", plot = F)
      out$Ref.p.val     <- t2$p.value
      out$Ref.stat      <- t2$statistic
      out$Ref.p.greater <- testDispersion(res, type = "DHARMa", plot = F,
                                          alternative = "greater")$p.value

      return(unlist(out))
    }

    out <- runBenchmarks(calculateStatistics, controlValues = phi,
                         nRep = 1000, parallel = T, exportGlobal = T)
    out.bin[[length(out.bin) + 1]] <- out
  }
}

names(out.bin) <- as.vector(unite(expand.grid(intercept,sampleSize), "sim"))$sim

# saving sim results
save(out.bin, sampleSize, intercept, phi,
     file = here("data", "3_glmBin_power_overdisp2.Rdata"))



###################
##### Poisson #####
###################
# 1) Simulating 1000 Poisson datasets with different sample sizes, intercepts and
#    overdispersion (negative binomial with constant dispersion factor phi).
# 2) fitting them to correct (Poisson) GLM models
# 3) calculating power for the dispersion tests
#     - Pearson-chisq,
#     - DHARMa default (simulated residuals)
#     - DHARMa refit (boostrapped Pearson residuals)


# varying parameters (same grid as 3_glm_power_underdisp.R)
sampleSize = c(10,50,100,200,500,1000)
intercept <- c(-1.5,0,1.5,3)


out.pois <- list()
for (k in sampleSize){
  for (i in intercept){

    # function to varying overdispersion (phi)
    calculateStatistics <- function(control = 1){
      # data
      testData <- createData_ud(phi = control,
                                overdispersionModel = "NB",
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
      out$Pear.p.val     <- t2$p.value
      out$Pear.stat      <- t2$statistic
      out$Pear.p.greater <- testDispersion(fittedModel, plot = F, type = "PearsonChisq",
                                           alternative = "greater")$p.value

      # DHARMa default residuals
      res <- simulateResiduals(fittedModel)
      t2 <- testDispersion(res, type = "DHARMa", plot = F)
      out$DHA.p.val     <- t2$p.value
      out$DHA.stat      <- t2$statistic
      out$DHA.p.greater <- testDispersion(res, type = "DHARMa", plot = F,
                                          alternative = "greater")$p.value

      # DHARMa refit residuals -> bootstrapped Pearson
      res <- simulateResiduals(fittedModel, refit = T)
      t2 <- testDispersion(res, type = "DHARMa", plot = F)
      out$Ref.p.val     <- t2$p.value
      out$Ref.stat      <- t2$statistic
      out$Ref.p.greater <- testDispersion(res, type = "DHARMa", plot = F,
                                          alternative = "greater")$p.value

      return(unlist(out))
    }

    out <- runBenchmarks(calculateStatistics, controlValues = phi,
                         nRep = 1000, parallel = T, exportGlobal = T)
    out.pois[[length(out.pois) + 1]] <- out
  }
}

names(out.pois) <- as.vector(unite(expand.grid(intercept,sampleSize), "sim"))$sim


# saving sim results
save(out.pois, sampleSize, intercept, phi,
     file = here("data", "3_glmPois_power_overdisp2.Rdata"))
