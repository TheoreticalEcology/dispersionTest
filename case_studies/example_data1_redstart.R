#####      Study case 1: redstarts breeding pairs      #####
# Melina Leite
# March 2026

# The data is availabe in the blmeco package
# for details of the dataset and references see the main text and help(redstart)


# Loading Packages ----
library(blmeco) # dataset
library(DHARMa) # residual tests
library(MASS) # glm.nb negative binomial GLM
library(glmmTMB) # spatial GLM
library(tidyverse) # plotting/manipulating data
library(broom); library(broom.mixed) # tidy summary tables


### Dataset ----
data("redstart")

redstart$elevation.z <- as.numeric(scale(redstart$elevation))
redstart$elevation.z2 <- redstart$elevation.z^2
redstart$x.s <- scale(redstart$x)
redstart$y.s <- scale(redstart$y)



#### Poisson model ########

modPois <- glm(counts ~ elevation.z + elevation.z2 + forests,
            data= redstart,
            family=poisson())
summary(modPois)

###### dispersion tests -----

resPois <- simulateResiduals(modPois)
resPoisRefit <- simulateResiduals(modPois, refit = T) 

simPois <- testDispersion(resPois) # simulation-based residual variance
pearPois <- testDispersion(resPois, type = "PearsonChisq") # parametric Pearson
npearPois <- testDispersion(resPoisRefit) # non-parametric Pearson


###### Residual spatial autocorrelation tests -----

spatPois <- testSpatialAutocorrelation(resPois, x=~x, y=~y)




#### Spatial Poisson model ########

# arranging data
redstart$pos <- numFactor(redstart$x.s, redstart$y.s)
redstart$group <- factor(rep(1, nrow(redstart)))

modSPois <- glmmTMB(counts ~ elevation.z + elevation.z2 + forests + 
                      exp(pos + 0 | group),
                    data= redstart, family=poisson())

summary(modSPois)

###### dispersion tests -----

resSPois <- simulateResiduals(modSPois)
resSPoisRefit <- simulateResiduals(modSPois, refit = T) # takes time! 

simSPois <- testDispersion(resSPois) # simulation-based residual variance test
pearSPois <- testDispersion(resSPois, type = "PearsonChisq", 
                            alternative = "greater") # parametric Pearson
npearSPois <- testDispersion(resSPoisRefit) # non-parametric Pearson


###### Residual spatial autocorrelation tests -----

resspatSPois <- simulateResiduals(modSPois, rotation = "estimated")
spatSPois <- testSpatialAutocorrelation(resspatSPois, x=~x.s, y=~y.s)




#### Negative binomial model ########

modNB <- glm.nb(counts ~ elevation.z + elevation.z2 + forests,
            data= redstart)
summary(modNB)

###### dispersion tests -----

resNB <- simulateResiduals(modNB)
resNBRefit <- simulateResiduals(modNB, refit = T) 

simNB <- testDispersion(resNB) 
pearNB <- testDispersion(resNB, type = "PearsonChisq")
npearNB <- testDispersion(resNBRefit) # non-parametric Pearson


###### Residual spatial autocorrelation tests -----

resspatNB <- simulateResiduals(modNB, rotation = "estimated")
spatNB <- testSpatialAutocorrelation(resspatNB, x=~x, y=~y)




#### RESULTS and tables #####

###### Coefficients ------

coefs <- bind_rows(list(Poisson = tidy(modPois),
                        `negative binomial` = tidy(modNB), 
                        `spatial Poisson` = tidy(modSPois, effects="fixed")),
                   .id="model") %>%
  dplyr::select(-effect, -component) %>% mutate_if(is.numeric, round, 3)
coefs
#knitr::kable(coefs) 


###### Dispersion tests ------
tests <- c("sim-based res. variance", "parametric Pearson", "nonparametric Pearson")
dPoisson = data.frame(test=c(tests),
                     Dispersion=c(simPois$statistic, pearPois$statistic, 
                                  npearPois$statistic),
                     p.value=c(simPois$p.value, pearPois$p.value, 
                               npearPois$p.value))
dSPois = data.frame(test=tests,
                    Dispersion=c(simSPois$statistic, pearSPois$statistic, 
                                npearSPois$statistic),
                    p.value=c(simSPois$p.value, pearSPois$p.value, 
                              npearSPois$p.value))
dNB = data.frame(test=tests,
                 Dispersion=c(simNB$statistic, pearNB$statistic, npearNB$statistic),
                 p.value=c(simNB$p.value, pearNB$p.value, npearNB$p.value))


dispersion_tests <- bind_rows(list(Poisson = dPoisson,
                                   `negative binomial` = dNB,
                                   `spatial Poisson` = dSPois),
                             .id="model") %>% 
  mutate(Dispersion= round(Dispersion, 2), p.value=round(p.value, 3))
dispersion_tests
#knitr::kable(dispersion_tests) 


###### Residual spatial autocorrelation tests ------

spatTests <- data.frame(model = c("Poisson", "spatial Poisson", 
                                  "negative binomial"),
          MoransI = c(spatPois$statistic[1], spatSPois$statistic[1],
                    spatNB$statistic[1]),
          Expected = c(spatPois$statistic[2], spatSPois$statistic[2],
                     spatNB$statistic[2]),
          SD = c(spatPois$statistic[3], spatSPois$statistic[3],
               spatNB$statistic[3]),
          p.value=c(spatPois$p.value, spatSPois$p.value, spatNB$p.value)) %>%
  mutate_if(is.numeric, round, 4) %>%
  arrange(model)
spatTests
knitr::kable(spatTests) 

