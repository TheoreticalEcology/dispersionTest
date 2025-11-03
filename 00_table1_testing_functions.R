# Dispersion Tests paper ----
# testing table 1 functions found in R for dispersion tests

# september 2025

# Data/models for testing ----
library(DHARMa)
set.seed(123)
testData1 <- createData(randomEffectVariance = 0, fixedEffects = 0) # no effect
testData2 <- createData(randomEffectVariance = 0, fixedEffects = 1.2) #  effect

glm1 <- glm(observedResponse ~ Environment1, data = testData2,
            family = "poisson")
glmnb1 <- MASS::glm.nb(observedResponse ~ Environment1, data = testData2)
glmm1 <- lme4::glmer(observedResponse ~ Environment1 + (1|group), data = testData2,
            family = "poisson")


# 1. Likelihood ratio tests ----

## 1.1 Likelihood ratio tests ----

pscl::odTest(glmnb1)

DCluster::test.nb.pois(glmnb1, glm1)

anova(glm1, glmnb1, test="LRT")

lmtest::lrtest(glm1, glmnb1)



# 2. Score tests ----

DCluster::DeanB(glm1)

DCluster::DeanB2(glm1)

Rfast2::overdispreg.test(testData2$observedResponse, testData2$Environment1)



# 3. Residual dispersion ----

DHARMa::testDispersion(glm1, type="Pearson")

performance::check_overdispersion(glm1)

msme::P__disp(glm1)

RVAideMemoire:: overdisp.glmer(glmm1)

overdisp::overdisp(testData2, 2, 3)

AER::dispersiontest(glm1)



# 4. Response variance ----

DHARMa::testDispersion(glm1, type="DHARMa")



# 5. OTHERS ----


## 5.1 vcd::goodfit ----
# good fit is not based on model but test the fit of the response variable to 
# the Poisson/binomial/negative binomial models
# decided not to include in the paper because it is not a dispersion test and 
# not based o models.

library(vcd)
?goodfit
gfit1 <- goodfit(testData1$observedResponse)
summary(gfit1)
gfit2 <- goodfit(testData2$observedResponse)
summary(gfit2)


## 5.2 gcc::qcc.overdispesion.test ----
# only for the response variable, it compares the observed variance divided
# by the theoretical variance in the response data * (N-1)
# test it with X2

library(qcc)
qtest1 <- qcc.overdispersion.test(testData1$observedResponse)
qtest1
qtest2 <- qcc.overdispersion.test(testData2$observedResponse)
qtest2




