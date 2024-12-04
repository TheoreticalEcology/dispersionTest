### Dispersion Tests Project
## Melina Leite
# Nov 24

library(DHARMa)
library(tidyverse)
library(here)

#####################################
##### Instructions & Simulating #####
#####################################

## KS GOF test to compare Pearson Stats and Chi-squared distribution


# 1) Simulating 1000 Poisson datasets with different sample sizes and intercepts
# -sampleSize: c(10,50,100,500)
# -intercept:  c(-3,-1,0,2,4)
# 2) fitting them to correct GLM model
# 3) calculating pearson statistics of the pearson residuals (sum(resˆ2)) 
# 4) comparing distribution of these statistics with the Chi-squared distribution with the same DF.
# 5) Repeating these steps 100 times to get the proportion of significant results per combination of parameters.


# sampleSizes
sampleSize = c(10,50,100,500)
intercept <- c(-3,-1,0,1,3)


# KS TEST
ks.p <- function(x, df) ks.test(x,"pchisq", df=df)$p.value

final.res <- list()
final.sims <- list()

for (k in 1:100) { # MANY SIMULATIONS TO HAVE A PROP OF SIG RESULTS
  set.seed(k)  
  result <- list()
  sims <- list()
  
  for (i in intercept){
    
    calculateStatistics <- function(control = 10){
      testData <- DHARMa::createData(sampleSize = control, intercept = i,
                                     numGroups = 10,
                                     randomEffectVariance = 0,
                                     family = poisson())
      
      fittedModel <- stats::glm(observedResponse ~ Environment1, 
                                data = testData, family = poisson()) 
      
      pearson <- residuals(fittedModel, "pearson")
      
      out <- list()
      out$Pear.stat <- sum(residuals(fittedModel, "pearson")^2)
      out$rdf <- df.residual(fittedModel)
      return(unlist(out))
    }
    
    out <- runBenchmarks(calculateStatistics, controlValues = sampleSize,
                         nRep=1000, parallel = T, exportGlobal = T)
    #simulations
    stats <- as.data.frame(apply(as.data.frame(out$simulations), 2, unlist))
    stats$intercept <- i
    sims <- rbind(stats, sims)
    
    #results
    res <- stats %>% group_by(controlValues) %>%
      summarise(ks.p = ks.p(Pear.stat,rdf))
    res$intercept <- i
    result <- rbind(result,res)
  }
  final.res[[k]]  <- result
  final.sims[[k]] <- sims
}

# saving results
save(final.res, final.sims, file = here("data", "1_glmPois_pearsonChisq.Rdata"))
#load(here("data", "1_glmPois_pearsonChisq.Rdata"))



### FIGURES
#
# library(cowplot);
# theme_set(theme_cowplot())
# library(patchwork)
# 
# ## Summary
# 
# # manipulaing results
# names(final.res) <- 1:length(final.res)
# 
# final.pois <- bind_rows(final.res, .id="sim") %>% group_by(controlValues,intercept) %>%
#   summarise(ks.sig = sum(ks.p<0.05))
# 
# for (i in 1:nrow(final.pois)) {
#   confs <- binom.test(final.pois$ks.sig[i], 100)$conf.int
#   final.pois$conf.low[i] <-  round(confs[1],2)
#   final.pois$conf.up[i] <-  round(confs[2],2)
# }
# 
# ggplot(final.pois, aes(y=ks.sig/100, x=as.factor(controlValues),
#                   col=as.factor(intercept))) +
#   geom_point(position = position_dodge(width=0.4)) +
#   geom_errorbar(aes(ymin=conf.low, ymax=conf.up),width = 0.1,
#                 position = position_dodge(width=0.4)) +
#   geom_hline(yintercept = 0.05, linetype="dashed") +
#   scale_color_discrete("intercept")+
#   xlab("sampleSize") + ylab("Prop of significant KS test") +
#   labs(title="Poisson", tag = "A)") +
#   theme(panel.border  = element_rect(color = "black"))
# 
# ggsave(here("figures", "1_glmPois_pearsonChisq.jpeg"), width = 6, height = 8)
# 
# 
# # trick to include the binomial plot in the same figure / whe the other code is also loaded
# # +
# #  pbin + plot_layout(ncol=1, guides="collect")
# #ggsave(here("figures", "1_glmBOTH_pearsonChisq.jpeg"), width = 6, height = 8)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # shape distributions
# 
# sims <- bind_rows(final.sims, .id="sim") 
# 
# # all simulations
# 
# sims %>% filter(controlValues == 10) %>%
#   ggplot(aes(x=Pear.stat, group=sim)) +
#   geom_density(col=2, alpha=0.5) +
#   facet_grid(controlValues~intercept)+
#   stat_function(fun = dchisq, args = list(df = 8), col="black")+
#   xlab("") + ylab("Density") +
#   theme(panel.border  = element_rect(color = "black"),
#         legend.position = "none") +
#   
#   sims %>% filter(controlValues == 50) %>%
#   ggplot(aes(x=Pear.stat, group=sim)) +
#   geom_density(col=2, alpha=0.5) +
#   facet_grid(controlValues~intercept)+
#   stat_function(fun = dchisq, args = list(df = 48), col="black")+
#   xlab("") + ylab("Density") +
#   theme(panel.border  = element_rect(color = "black"),
#         legend.position = "none") +
#   
#   sims %>% filter(controlValues == 100) %>%
#   ggplot(aes(x=Pear.stat, group=sim)) +
#   geom_density(col=2, alpha=0.5) +
#   facet_grid(controlValues~intercept)+
#   stat_function(fun = dchisq, args = list(df = 98), col="black")+
#   xlab("") + ylab("Density") +
#   theme(panel.border  = element_rect(color = "black"),
#         legend.position = "none") +
#   
#   sims %>% filter(controlValues == 500) %>%
#   ggplot(aes(x=Pear.stat, group=sim)) +
#   geom_density(col=2, alpha=0.5) +
#   facet_grid(controlValues~intercept)+
#   stat_function(fun = dchisq, args = list(df = 498), col="black")+
#   xlab("Pearson Stats") + ylab("Density") +
#   theme(panel.border  = element_rect(color = "black"),
#         legend.position = "none") +
#   plot_layout(ncol=1) + plot_annotation(title="Poisson: Pearson Stats distribution")
# ggsave(here("figures", "1_glmPois_pearsonChisq_distrib.jpeg"), width=15,
#        height=15)
# 
# 
# # Mean distribution of 100 simulations
# 
# sims.mean <- sims %>% group_by(sim, controlValues,intercept) %>% 
#   arrange(sim, controlValues, intercept, Pear.stat) %>%
#   mutate(n = 1:n()) %>%
#   group_by(controlValues,intercept,n) %>%
#   summarise(mean.pearson = mean(Pear.stat))
# 
# sims.mean %>% filter(controlValues == 10) %>%
#   ggplot(aes(x=mean.pearson)) +
#   facet_grid(controlValues~intercept)+
#   geom_density(col=2) +
#   stat_function(fun = dchisq, args = list(df = 8), col="black") +
#   xlab("") + ylab("Density") +
#   theme(panel.border  = element_rect(color = "black"),
#         legend.position = "none") +
#   
#   
#   sims.mean %>% filter(controlValues == 50) %>%
#   ggplot(aes(x=mean.pearson)) +
#   facet_grid(controlValues~intercept)+
#   geom_density(col=2) +
#   stat_function(fun = dchisq, args = list(df = 48), col="black") +
#   xlab("") + ylab("Density") +
#   theme(panel.border  = element_rect(color = "black"),
#         legend.position = "none") +
#   
#   
#   sims.mean %>% filter(controlValues == 100) %>%
#   ggplot(aes(x=mean.pearson)) +
#   facet_grid(controlValues~intercept)+
#   geom_density(col=2) +
#   stat_function(fun = dchisq, args = list(df = 98), col="black") +
#   xlab("") + ylab("Density") +
#   theme(panel.border  = element_rect(color = "black"),
#         legend.position = "none") +
#   
#   sims.mean %>% filter(controlValues == 500) %>%
#   ggplot(aes(x=mean.pearson)) +
#   facet_grid(controlValues~intercept)+
#   geom_density(col=2) +
#   stat_function(fun = dchisq, args = list(df = 498), col="black") +
#   xlab("Pearson Stats") + ylab("Density") +
#   theme(panel.border  = element_rect(color = "black"),
#         legend.position = "none") +
#   
#   plot_layout(ncol=1) + plot_annotation(title="Poisson: Mean Pearson Stats distribution")
# 
# ggsave(here("figures", "1_glmPois_pearsonChisq_distrib_MEAN.jpeg"), width=15,
#        height=15)













