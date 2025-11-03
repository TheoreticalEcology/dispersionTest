### Dispersion Tests Project
## 
# Mar 25

library(DHARMa)
library(lme4)
library(here)

# testing the approximation alternatives for DF in GLMMs

# only for Poisson

sampleSize <- c(200, 500, 1000)
intercept <- c(-1.5,0,1.5)
overdispersion = seq(0,1,0.1)
ngroups <-  c(10, 50, 100)
nRep = 1000


# approximate DF for glmms
source(here("functions_others", "approximateDFpearson.R"))

out.pois <- list()

for(s in ngroups) {
  for (k in sampleSize){
    for (i in intercept){

    # function to varying sampleSize
    calculateStatistics <- function(control = 0){
      # data
      testData <- createData(overdispersion = control,
                             sampleSize = k,
                             intercept = i,
                             numGroups = s,
                             randomEffectVariance = 1,
                             family = poisson())
      # model
      fittedModel <-glmer(observedResponse ~ Environment1 + (1|group), 
                          data = testData, family = poisson())
      #results
      out <- list()

      # Pearson naive
      pearNaive <- approximateDFpearson(fittedModel, testData, type="naive", 
                                        alternative = "two.sided")
      out$PearNaive.p <- pearNaive$p.value
      out$PearNaive.rdf <- pearNaive$parameter
      out$PearNaive.statistic <- pearNaive$statistic
      
      # Pearson Sat
      pearSat <- approximateDFpearson(fittedModel, testData, type="sat", 
                                        alternative = "two.sided")
      out$PearSat.p <- pearSat$p.value
      out$PearSat.rdf <- pearSat$parameter
      out$PearSat.statistic <- pearSat$statistic
      
      # Pearson KR
      pearKR <- approximateDFpearson(fittedModel, testData, type="KR", 
                                        alternative = "two.sided")
      out$PearKR.p <- pearKR$p.value
      out$PearKR.rdf <- pearKR$parameter
      out$PearKR.statistic <- pearKR$statistic
      
      # Pearson KR2
      pearKR2 <- approximateDFpearson(fittedModel, testData, type="KR2", 
                                     alternative = "two.sided")
      out$PearKR2.p <- pearKR2$p.value
      out$PearKR2.rdf <- pearKR2$parameter
      out$PearKR2.statistic <- pearKR2$statistic

      return(unlist(out))
    }

    out <- runBenchmarks(calculateStatistics, controlValues = overdispersion,
                         nRep = nRep, parallel = T, exportGlobal = T)
    out.pois[[length(out.pois) + 1]] <- out
    names(out.pois)[length(out.pois)] <- paste(s, k, i, sep="_")
    # saving sim results
    save(out.pois, file=here("data","8_approximateDFpearson.Rdata"))
  }
 }
}
