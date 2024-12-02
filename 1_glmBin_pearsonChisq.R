### Dispersion Tests Project
## Melina Leite
# Nov 24

library(DHARMa)
library(tidyverse)
library(here)

#####################################
##### Instructions & Simulations #####
#####################################

## KS GOF test to compare Pearson Stats and Chi-squared distribution


# 1) Simulating 1000 Binomial proportion datasets with different sample sizes and intercepts. Fixing number of trials in 20.
#     -sampleSize: c(10,50,100,500)
#     -intercept:  c(-3,-1,0,2,4)
# 2) fitting them to correct GLM model
# 3) calculating pearson statistics of the pearson residuals (sum(resˆ2)) 
# 4) comparing distribution of these statistics with the Chi-squared distribution with the same DF.
# 5) Repeating these steps 100 times to get the proportion of significant results per combination of parameters.


# function to varying N
calculateStatistics <- function(control = 10){
  testData <- DHARMa::createData(sampleSize = control,
                                 numGroups = 1,
                                 binomialTrials = 20,
                                 family = binomial())
  
  fittedModel <- stats::glm(cbind(observedResponse1,observedResponse0) ~
                              Environment1, data = testData, family = binomial()) 
  
  pearson <- residuals(fittedModel, "pearson")
  
  out <- list()
  out$Pear.stat <- sum(residuals(fittedModel, "pearson")^2)
  out$rdf <- df.residual(fittedModel)
  return(unlist(out))
}

# sampleSizes
sampleSize = c(10,50,100,500)

# intercept (in a loop)
intercept <- c(-3,-1,0,2,4)


# KS TEST
ks.p <- function(x, df) ks.test(x,"pchisq", df=df)$p.value

final.res <- list()
final.sims <- list()

for (k in 1:100) { # MANY SIMULATIONS TO HAVE A PROP OF SIG RESULTS
  set.seed(k)  
  result <- list()
  sims <- list()
  
  for (i in intercept){
    out <- runBenchmarks(calculateStatistics, controlValues = sampleSize,
                         intercept = i,
                         nRep=1000, parallel = T)
    #simulations
    stats <- as.data.frame(apply(as.data.frame(out$simulations), 2, unlist))
    stats$intercept <- i
    sims <- rbind(stats,sims)
    
    #results
    res <- stats %>% group_by(controlValues) %>%
      summarise(ks.p = ks.p(Pear.stat,rdf))
    res$intercept <- i
    result <- rbind(result,res)
  }
  final.res[[k]] <- result
  final.sims[[k]] <- sims
}

# saving results
save(final.res, final.sims, file = here("data", "1_glmBin_pearsonChisq.Rdata"))
#load(here("data", "1_glmBin_pearsonChisq.Rdata"))





### FIGURES
#
library(cowplot);
theme_set(theme_cowplot())
library(patchwork)

## Summary

# manipulaing results
names(final.res) <- 1:length(final.res)

final.bin <- bind_rows(final.res, .id="sim") %>% group_by(controlValues,intercept) %>%
  summarise(ks.sig = sum(ks.p<0.05))

for (i in 1:nrow(final.bin)) {
  confs <- binom.test(final.bin$ks.sig[i], 100)$conf.int
  final.bin$conf.low[i] <-  round(confs[1],2)
  final.bin$conf.up[i] <-  round(confs[2],2)
}

pbin <- ggplot(final.bin, aes(y=ks.sig/100, x=as.factor(controlValues),
                  col=as.factor(intercept))) +
  geom_point(position = position_dodge(width=0.4)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up),width = 0.1,
                position = position_dodge(width=0.4)) +
  geom_hline(yintercept = 0.05, linetype="dashed") +
  ylim(0,1)+
  scale_color_discrete("intercept")+
  xlab("sampleSize") + ylab("Prop of significant KS test") +
  labs(title="Binomial", tag = "B)") +
  theme(panel.border  = element_rect(color = "black"))
pbin
ggsave(here("figures", "1_glmBin_pearsonChisq.jpeg"), width = 6, height = 4)


# shape distributions

sims <- bind_rows(final.sims, .id="sim") 

# all simulations

sims %>% filter(controlValues == 10) %>%
  ggplot(aes(x=Pear.stat, group=sim)) +
  geom_density(col=2, alpha=0.5) +
  facet_grid(~intercept)+
  stat_function(fun = dchisq, args = list(df = 8), col="black")+
  ggtitle("N=10") +
  xlab("Pearson Stats") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none") +
  
sims %>% filter(controlValues == 50) %>%
  ggplot(aes(x=Pear.stat, group=sim)) +
  geom_density(col=2, alpha=0.5) +
  facet_grid(~intercept)+
  stat_function(fun = dchisq, args = list(df = 48), col="black")+
  ggtitle("N=50") +
  xlab("Pearson Stats") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none") +
  
sims %>% filter(controlValues == 100) %>%
  ggplot(aes(x=Pear.stat, group=sim)) +
  geom_density(col=2, alpha=0.5) +
  facet_grid(~intercept)+
  stat_function(fun = dchisq, args = list(df = 98), col="black")+
  ggtitle("N=100") +
  xlab("Pearson Stats") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none") +

sims %>% filter(controlValues == 500) %>%
  ggplot(aes(x=Pear.stat, group=sim)) +
  geom_density(col=2, alpha=0.5) +
  facet_grid(~intercept)+
  stat_function(fun = dchisq, args = list(df = 498), col="black")+
  ggtitle("N=500") +
  xlab("Pearson Stats") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none") +
plot_layout(ncol=1) + plot_annotation(title="Pearson Stats distribution")
ggsave(here("figures", "1_glmBin_pearsonChisq_distrib.jpeg"), width=15,
       height=15)


# Mean distribution of 100 simulations

sims.mean <- sims %>% group_by(sim, controlValues,intercept) %>% 
  arrange(sim, controlValues, intercept, Pear.stat) %>%
  mutate(n = 1:n()) %>%
  group_by(controlValues,intercept,n) %>%
  summarise(mean.pearson = mean(Pear.stat))

sims.mean %>% filter(controlValues == 10) %>%
  ggplot(aes(x=mean.pearson)) +
  facet_grid(~intercept)+
  geom_density(col=2) +
  stat_function(fun = dchisq, args = list(df = 8), col="black") +
  ggtitle("N=10") +
  xlab("Pearson Stats") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none") +


sims.mean %>% filter(controlValues == 50) %>%
  ggplot(aes(x=mean.pearson)) +
  facet_grid(~intercept)+
  geom_density(col=2) +
  stat_function(fun = dchisq, args = list(df = 48), col="black") +
  ggtitle("N=50") +
  xlab("Pearson Stats") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none") +


sims.mean %>% filter(controlValues == 100) %>%
  ggplot(aes(x=mean.pearson)) +
  facet_grid(~intercept)+
  geom_density(col=2) +
  stat_function(fun = dchisq, args = list(df = 98), col="black") +
  ggtitle("N=100") +
  xlab("Pearson Stats") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none") +

sims.mean %>% filter(controlValues == 500) %>%
  ggplot(aes(x=mean.pearson)) +
  facet_grid(~intercept)+
  geom_density(col=2) +
  stat_function(fun = dchisq, args = list(df = 498), col="black") +
  ggtitle("N=500") +
  xlab("Pearson Stats") + ylab("Density") +
  theme(panel.border  = element_rect(color = "black"),
        legend.position = "none") +
  
plot_layout(ncol=1) + plot_annotation(title="Mean Pearson Stats distribution")

ggsave(here("figures", "1_glmBin_pearsonChisq_distrib_MEAN.jpeg"), width=15,
       height=15)

















