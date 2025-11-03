### Dispersion Tests Project
## 
# jan 25

library(DHARMa)
library(tidyverse)
library(here)
library(lme4)


# varying parameters
overdispersion <- seq(0,1,0.10)
intercept <- c(-3,-1.5,0,1.5,3) # reinclude -3
## sample size depend on ngroups (see it within the loop)
ngroups <- 50 # 10, 50 and 100
nRep = 1000


#########################
##### Poisson GLMM #####
#########################


out.pois <- list()
#load(here("data", "5_glmmPois_power.Rdata")) # loading out.pois

for(m in ngroups){
  
  if(m == 10) sampleSize <- c(50,100,200,500,1000)
  if(m == 50) sampleSize <- c(100,200,500,1000) 
  if(m == 100) sampleSize <- c(200,500,1000)
  
  # if(m ==10) intercept <- c(-1.5,0,1.5,3) # exclude -3 for m=10 (doesn't work)
  
  for (k in sampleSize){
    
    for (i in intercept){
      
      # function to varying sampleSize
      calculateStatistics <- function(control = 0){
        # data
        testData <- createData(overdispersion = control,
                               sampleSize = k,
                               intercept = i,
                               numGroups = m,
                               randomEffectVariance = 1,
                               family = poisson())
        # model
        fittedModel <- glmer(observedResponse  ~ 
                               Environment1 + (1|group),
                             data = testData, family = poisson()) 
        #results
        out <- list()
        
        # pearson residual
        out$Pear.p.val <- testDispersion(fittedModel, plot = F, 
                                         type="PearsonChisq")$p.value
        out$Pear.stat <- testDispersion(fittedModel, plot = F, 
                                        type="PearsonChisq")$statistic
        
        # DHARMa default unconditional
        res <- simulateResiduals(fittedModel, refit=F, re.form=NA)
        out$dhaUN.p.val <- testDispersion(res, type = "DHARMa",plot = F)$p.value
        out$dhaUN.stat  <- testDispersion(res, type = "DHARMa",plot = F)$statistic
        
        # DHARMa default conditional
        res <- simulateResiduals(fittedModel, refit=F, re.form=NULL)
        out$dhaCO.p.val <- testDispersion(res, type = "DHARMa",plot = F)$p.value
        out$dhaCO.stat  <- testDispersion(res, type = "DHARMa",plot = F)$statistic
        
        # DHARMa refit
        # res <- simulateResiduals(fittedModel, refit=T, re.form=NA)
        # out$refUN.p.val <- testDispersion(res, plot = F, type = "DHARMa")$p.value
        # out$refUN.stat <- testDispersion(res, plot = F, type = "DHARMa")$statistic
        
        # DHARMa refit conditional
        res <- simulateResiduals(fittedModel, refit=T, re.form=NULL)
        out$refCO.p.val <- testDispersion(res, plot = F, type = "DHARMa")$p.value
        out$refCO.stat <- testDispersion(res, plot = F, type = "DHARMa")$statistic
        
        
        return(unlist(out))
      }
      
      out <- runBenchmarks(calculateStatistics, controlValues = overdispersion,
                           nRep=nRep, parallel = T, exportGlobal = T)
      out.pois[[length(out.pois) + 1]] <- out
      names(out.pois)[length(out.pois)] <- paste(m, k, i, sep="_")
      # saving sim results
      save(out.pois, file=here("data","5_glmmPois_power_50.Rdata"))
    }
  }
}




