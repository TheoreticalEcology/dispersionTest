library(DHARMa)
returnStatistics <- function(control = 1){
  testData = createData(sampleSize = 200, family = poisson(),
                        overdispersion = control, fixedEffects = 0,
                        randomEffectVariance = 0)
  fittedModel <- glm(observedResponse ~ Environment1, data = testData, family = poisson())
  x = summary(fittedModel)
  res <- simulateResiduals(fittedModel = fittedModel, n = 250)
  out <- c("Type I error GLM slope" = x$coefficients[2,4],
           "DHARMa testDispersion" = testDispersion(res, plot = FALSE)$p.value)
  return(out)
}

out = runBenchmarks(returnStatistics, controlValues = seq(0, 1.5, 0.05), nRep = 500,
                    parallel = T)
plot(out, xlab = "Added dispersion sd", ylab = "Prop significant", main = "n = 200")
abline(h=0.05, col="blue")
abline(v=0.34, col="red")


out$summaries$propSignificant


returnStatistics <- function(control = 1){
  testData = createData(sampleSize = 200, family = poisson(),
                        overdispersion = control, fixedEffects = 0,
                        randomEffectVariance = 0)
  fittedModel <- glm(observedResponse ~ Environment1, data = testData, family = poisson())
  x = summary(fittedModel)
  res <- simulateResiduals(fittedModel = fittedModel, n = 250)
  out <- c("Type I error GLM slope" = x$coefficients[2,4],
           "Pearson testDispersion" = testDispersion(res, type = "PearsonChisq", 
                                                    plot = FALSE)$p.value)
  return(out)
}

out2 = runBenchmarks(returnStatistics, controlValues = seq(0, 1.5, 0.05), nRep = 500,
                    parallel = T)
plot(out2, xlab = "Added dispersion sd", ylab = "Prop significant", main = "n = 200")
abline(h=0.05, col="blue")


### gamlss

dado = createData(randomEffectVariance = 0)

library(gamlss)
mod <- gamlss(observedResponse ~ Environment1, data=dado, family=PO())

plot(mod)
wp(mod)
wp(mod, xvar=dado$Environment1, n.inter=2)
dtop(mod)
rqres.plot(mod)
rqres.plot(mod, type="QQ")

mod2 <- glm(observedResponse ~ Environment1, data = dado, family = poisson())
wp(resid = residuals(mod2, type="pearson"))




















































