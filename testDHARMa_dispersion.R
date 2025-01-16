## testing a differnt approach for the DHARMa default dispersion test

library(DHARMa)

## FIRST examinating testDispersion ####

#### expectedVar ####

#expectedVar = sd(simulationOutput$simulatedResponse)^2

# is the variance in the whole simulation matrix 
# it can't use var because var function won't treat the matrix as a vector
# it is a way to normalization - to ensure it has similar scale to Pearson residuals (but I didn't fully understand it)



####  function spread() ####

# line 432 of the `tests.R` file for the testDispersion function 
# function spread calculates the variance of the raw residuals in a different way than what we would calculate through a linear model
# instead of calculating the sum((observed - fitted)^2)/(n-1)
# it does sum(((observed - fitted) - mean(observed-fitted))^2)/(n-1)
# because it uses the var(obs-fitted)

# recreating the testDispersion changing the spread function
 
testDispersion2 <- function(simulationOutput, alternative = c("two.sided", "greater", "less"), plot = T, type = c("DHARMa", "PearsonChisq"), ...){
  
  alternative <- match.arg(alternative)
  type <- match.arg(type)
  
  out = list()
  out$data.name = deparse(substitute(simulationOutput))
  
  if(type == "DHARMa"){
    
    simulationOutput = ensureDHARMa(simulationOutput, convert = "Model")
    
    # if(class(simulationOutput$fittedModel) %in% c("glmerMod", "lmerMod"){
    #   if(!"re.form" %in% names(simulationOutput$additionalParameters) & is.null(simulationOutput$additionalParameters$re.form)) message("recommended to run conditional simulations for dispersion test, see help")
    #}
    
    if(simulationOutput$refit == F){
      
      expectedVar = sd(simulationOutput$simulatedResponse)^2
      spread <- function(x) {
        sum((x - simulationOutput$fittedPredictedResponse)^2)/(length(x)-1)
        }
      out = testGeneric(simulationOutput, summary = spread, alternative = alternative, methodName = "DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated", plot = plot, ...)
      names(out$statistic) = "dispersion"
    } else {
      
      observed = tryCatch(sum(residuals(simulationOutput$fittedModel, type = "pearson")^2), error = function(e) {
        message(paste("DHARMa: the requested tests requires pearson residuals, but your model does not implement these calculations. Test will return NA. Error message:", e))
        return(NA)
      })
      if(is.na(observed)) return(NA)
      expected = apply(simulationOutput$refittedPearsonResiduals^2 , 2, sum)
      out$statistic = c(dispersion = observed / mean(expected))
      names(out$statistic) = "dispersion"
      out$method = "DHARMa nonparametric dispersion test via mean deviance residual fitted vs. simulated-refitted"
      
      p = getP(simulated = expected, observed = observed, alternative = alternative)
      
      out$alternative = alternative
      out$p.value = p
      class(out) = "htest"
      
      if(plot == T) {
        #plotTitle = gsub('(.{1,50})(\\s|$)', '\\1\n', out$method)
        xLabel = paste("Simulated values, red line = fitted model. p-value (",out$alternative, ") = ", out$p.value, sep ="")
        
        hist(expected, xlim = range(expected, observed, na.rm=T ), col = "lightgrey", main = "", xlab = xLabel, breaks = 20, cex.main = 1)
        abline(v = observed, lwd= 2, col = "red")
        
        main = ifelse(out$p.value <= 0.05,
                      "Dispersion test significant",
                      "Dispersion test n.s.")
        
        title(main = main, cex.main = 1,
              col.main = ifelse(out$p.value <= 0.05, "red", "black"))
      }
    }
    
  } else if(type == "PearsonChisq"){
    
    if("DHARMa" %in% class(simulationOutput)){
      model = simulationOutput$fittedModel
    }
    else model = simulationOutput
    
    if(! alternative == "greater") message("Note that the chi2 test on Pearson residuals is biased for MIXED models towards underdispersion. Tests with alternative = two.sided or less are therefore not reliable. If you have random effects in your model, I recommend to test only with alternative = 'greater', i.e. test for overdispersion, or else use the DHARMa default tests which are unbiased. See help for details.")
    
    rdf <- df.residual(model)
    rp <- getPearsonResiduals(model)
    Pearson.chisq <- sum(rp^2)
    prat <- Pearson.chisq/rdf
    if(alternative == "greater") pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
    else if (alternative == "less") pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=TRUE)
    else if (alternative == "two.sided") pval <- min(min(pchisq(Pearson.chisq, df=rdf, lower.tail=TRUE), pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)) * 2,1)
    
    out$statistic = prat
    names(out$statistic) = "dispersion"
    out$parameter = rdf
    names(out$parameter) = "df"
    out$method = "Parametric dispersion test via mean Pearson-chisq statistic"
    out$alternative = alternative
    out$p.value = pval
    class(out) = "htest"
    # c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
    return(out)
  }
  
  return(out)
}




##### testing it #####

# small test to see if the tests would differ if using the "corrected"version

calculateStatistics <- function(control = 0){
  # data
  testData <- DHARMa::createData(sampleSize = 100,
                                 overdispersion = control,
                                 intercept = 0,
                                 numGroups = 10,
                                 randomEffectVariance = 0,
                                 family = poisson())
  # model
  fittedModel <- stats::glm(observedResponse~ Environment1, data = testData, 
                            family = poisson()) 
  #results
  out <- list()
  
  # DHARMa default residuals
  res <- simulateResiduals(fittedModel)
  out$DHA.p     <- testDispersion(res, type = "DHARMa",plot = F)$p.value
  out$DHA.stat  <- testDispersion(res, type = "DHARMa",plot = F)$statistic
  
  # DHARMa "corrected" residuals
  out$DHAc.p    <- testDispersion2(res, plot = F, type = "DHARMa")$p.value
  out$DHAc.stat <- testDispersion2(res, plot = F, type = "DHARMa")$statistic
  
  return(unlist(out))
}

overdisp = c(0,0.5,1)
out <- runBenchmarks(calculateStatistics, controlValues = overdisp,
                     nRep=100, parallel = T, exportGlobal = T)

library(ggplot2)

## Dispersion statistics
out$simulations %>%
  ggplot(aes(x=DHA.stat.dispersion,y=DHAc.stat.dispersion)) +
  geom_point() +
  geom_smooth(method="lm", col="blue")+
  geom_abline(slope=1, intercept = 0, col="red")+
  facet_wrap(~controlValues, scales="free")

lm(DHAc.stat.dispersion ~ DHA.stat.dispersion, 
   data = out$simulations[out$simulations$controlValues==0,])
lm(DHAc.stat.dispersion ~ DHA.stat.dispersion, 
   data = out$simulations[out$simulations$controlValues==0.5,])
lm(DHAc.stat.dispersion ~ DHA.stat.dispersion, 
   data = out$simulations[out$simulations$controlValues==1,])


# pvalues
out$simulations %>%
  ggplot(aes(x=DHA.p,y=DHAc.p)) +
  geom_point() +
  facet_wrap(~controlValues, scales="free") +
  geom_smooth(method="lm", col="blue")+
  geom_abline(slope=1, intercept = 0, col="red")


lm(DHAc.p ~ DHA.p, 
   data = out$simulations[out$simulations$controlValues==0,])
lm(DHAc.p ~ DHA.p, 
   data = out$simulations[out$simulations$controlValues==0.5,])


## it seems it doesn't change the test and dispersion statistics


## SECOND: trying to include the sd of simul observations ####

#example to do by hand

testData <- createData(family=poisson(), 
                       intercept=-2,
                       randomEffectVariance = 0,
                       overdispersion = 0)
fittedModel <- glm(observedResponse  ~ Environment1,
                   data = testData, family = poisson())

simulationOutput <- simulateResiduals(fittedModel, n = 50)
testDispersion(simulationOutput)

simRes <- simulationOutput$simulatedResponse

# Observation-wise SD
simulationSD.obs <-  apply(simulationOutput$simulatedResponse,1,sd)

## numbero of unique simulated values per observation



look <- data.frame(obs = testData$observedResponse,
                   pred = simulationOutput$fittedPredictedResponse,
                   sdSim = simulationSD.obs,
                   nuniq = nuniq)




if(min(simulationSD.obs)==0){
  #some ideas
  # newSD <- min(simulationSD.obs[simulationSD.obs != 0])
 # newSD <- min(simulationSD.obs[simulationSD.obs != 0])/simulationOutput$nSim
  newSD <- newSD <- min(simulationSD.obs[simulationSD.obs != 0])/max(nuniq)
  
  simulationSD.obs[simulationSD.obs==0] <- newSD
}

## observed variance weighted by the simulatedSD for each observation


spread2 <- function(x){
  val <- (x - simulationOutput$fittedPredictedResponse)/simulationSD.obs
  (sum(val^2)/(length(x)-1))
}

(observed = spread2(simulationOutput$observedResponse))

simulated = apply(simulationOutput$simulatedResponse, 2, spread2)
mean(simulated)

hist(simulated); abline(v=observed, col="red")
(dispersion = observed/mean(simulated))

getP(simulated, observed, "two.sided")

testDispersion(simulationOutput, type="PearsonChisq")

### rodar alguns testes













