### Dispersion Tests Project
## Melina Leite
# Dez 24

library(DHARMa)
library(tidyverse)
library(here)

#####################################
##### Instructions & Simulations #####
#####################################

# 1) Simulating 1000 binomial prop datasets with different sample sizes, intercepts and overdispersion. Fixing the number of trials to 10
# 2) fitting them to correct GLM models
# 5) calculating power for the dispersion tests 
#     - Pearson-chisq, 
#     - DHARMa default (simulated residuals)
#     - DHARMa refit (boostrapped Pearson residuals)



# varying parameters
overdispersion <- seq(0,1,0.10)
sampleSize = c(10,20,50,100,200,500,1000,10000)
intercept <- c(-3,-1.5,0,1.5,3)


out.out <- list()
for (k in sampleSize){
  for (i in intercept){
  
  # function to varying sampleSize
  calculateStatistics <- function(control = 0){
    # data
    testData <- createData(overdispersion = control,
                                   sampleSize = k,
                                   intercept = i,
                                   numGroups = 10,
                                   randomEffectVariance = 0,
                                   binomialTrials = 10,
                                   family = binomial())
    # model
    fittedModel <-glm(cbind(observedResponse1, observedResponse0)  ~ 
                            Environment1, data = testData, family = binomial()) 
    #results
    out <- list()
    
    # pearson residual
    out$Pear.p.val <- testDispersion(fittedModel, plot = F, 
                                     type="PearsonChisq")$p.value
    out$Pear.stat <- testDispersion(fittedModel, plot = F, 
                                    type="PearsonChisq")$statistic
    # DHARMa default residuals
    res <- simulateResiduals(fittedModel)
    out$DHA.p.val <- testDispersion(res, type = "DHARMa",plot = F)$p.value
    out$DHA.stat  <- testDispersion(res, type = "DHARMa",plot = F)$statistic
    
    # DHARMa refit residuals -> bootstrapped Pearson
    res <- simulateResiduals(fittedModel, refit=T)
    out$Ref.p.val <- testDispersion(res, plot = F, type = "DHARMa")$p.value
    out$Ref.stat <- testDispersion(res, plot = F, type = "DHARMa")$statistic
    
    return(unlist(out))
  }
  
  out <- runBenchmarks(calculateStatistics, controlValues = overdispersion,
                       nRep=1000, parallel = T, exportGlobal = T)
  out.out[[length(out.out) + 1]] <- out
}
}

names(out.out) <- as.vector(unite(expand.grid(intercept,sampleSize), "sim"))$sim

# saving sim results
save(out.out,sampleSize,intercept,overdispersion, file=here("data", 
                                             "3_glmBin_power.Rdata"))



#####################################
#####          RESULTS          #####
#####################################
# 


# load(here("data", "3_glmBin_power.Rdata"))
# 
# 
# simuls <- list()
# for (i in 1:length(out.out)) {
#   sim <- out.out[[i]]$simulations
#   params <- strsplit(names(out.out)[[i]], "_")[[1]]
#   sim$intercept <- params[1]
#   sim$sampleSize <- params[2]
#   simuls <- rbind(simuls,sim)
# }
# names(simuls)[names(simuls)=="controlValues"] <- "overdispersion"
# 
# simuls <- simuls %>% pivot_longer(1:3, names_to = "test", values_to = "p.val") %>%
#   group_by(sampleSize,intercept,overdispersion, test) %>%
#   summarise(p.sig = sum(p.val<0.05,na.rm=T),
#             nsim = length(p.val[!is.na(p.val)]) ) 
# simuls$prop.sig <- simuls$p.sig/simuls$nsim
# simuls$intercept <- fct_relevel(simuls$intercept, "-3", "-1", "0", "1", "3")
# simuls$sampleSize <- fct_relevel(simuls$sampleSize, "10", "50", "100", "500")
# 
# ggplot(simuls, aes(x=overdispersion, y=prop.sig, col=test))+
#   geom_point(alpha=0.7) + geom_line(alpha=0.7) +
#   scale_color_discrete(
#     labels=c("Quantile Residuals", "Pearson Chi-squared",
#                               "Pearson Param. Bootstrap."))+
#   facet_grid(sampleSize~intercept) +
#   geom_hline(yintercept = 0.5, linetype="dotted") +
#   ggtitle("Binomial", subtitle = "1000 sim; Ntrials=10") +
#   theme(panel.background = element_rect(color="black"),
#         legend.position = "bottom")
# ggsave(here("figures", "3_glmBin_power.jpeg"), width=10, height = 8)














