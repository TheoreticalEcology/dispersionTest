### Dispersion Tests Project
## 
# Oct 25

set.seed(1)
library(DHARMa)
library(tidyverse)
library(here)
library(lme4)
library(parallel)
library(cowplot)
theme_set(theme_cowplot())

# calculating runtime for GLMM tests

# parametric pearson 
# nonparametric pearson
# simulation-based response variance - conditional

# parameters
overdispersion <- 0.5
intercept <- 0
sampleSize = 1000
ngroups <- 100 
nRep = 1000

# paralell
nCores <- 7

# função para rodar uma simulação
run_sim <- function(i) {
  # gerar dados
  testData <- createData(overdispersion = overdispersion,
                         sampleSize = sampleSize,
                         intercept = intercept,
                         numGroups = ngroups,
                         randomEffectVariance = 1,
                         family = poisson())
  
  # ajustar modelo
  fittedModel <- glmer(observedResponse ~ Environment1 + (1 | group),
                       data = testData, family = poisson())
  
  # medir tempos
  paramP <- system.time(testDispersion(simulateResiduals(fittedModel),
                                       plot = FALSE, type = "PearsonChisq"))[3]
  
  nonparamP <- system.time(testDispersion(simulateResiduals(fittedModel,
                                                            refit = TRUE,
                                                            re.form = NULL),
                                          plot = FALSE, type = "DHARMa"))[3]
  
  simbased <- system.time(testDispersion(simulateResiduals(fittedModel,
                                                           refit = FALSE,
                                                           re.form = NULL),
                                         plot = FALSE, type = "DHARMa"))[3]
  
  # retornar resultados
  return(c(paramP = paramP, nonparamP = nonparamP, simbased = simbased))
}

# configurar cluster
cl <- makeCluster(nCores)
clusterEvalQ(cl, {
  library(lme4)
  library(DHARMa)
})
clusterExport(cl, c("createData", "overdispersion", "intercept", "sampleSize", "ngroups"))

# rodar simulações em paralelo
results_list <- parLapply(cl, 1:nRep, run_sim)

# encerrar cluster
stopCluster(cl)

# combinar resultados em data.frame
times <- as.data.frame(do.call(rbind, results_list)) %>%
  pivot_longer(1:3, names_to = "test", values_to = "time_s") %>%
  mutate(test = fct_recode(test, `param. Pearson` = "paramP.elapsed", 
                           `nonparam. Pearson` ="nonparamP.elapsed", 
                           `sim-based variance`="simbased.elapsed")) %>%
  mutate(test = fct_relevel(test, "param. Pearson","sim-based variance",
                            "nonparam. Pearson"))
save(times, file=here("data", "5_glmmm_runtime_test.Rdata"))

# stats
times %>% group_by(test) %>% summarise(mean = mean(time_s)) 
2790/0.07

# plot runing time
times %>%
  ggplot(aes(x=test, y=time_s)) + 
  scale_y_log10()+
  geom_hline(yintercept = c(0.1,1,10), col="gray", linetype="dashed")+
  ggbeeswarm::geom_quasirandom() +
  stat_summary( geom = "crossbar", fun=median, col="red", ) +
  annotate("text", x=c(0.55, 1.55, 2.55), y=c(0.075, 0.08, 32),
           label=c(0.066, 0.072, 27.9), col="red")+
  xlab("Dispersion test") +
  ylab("Runtime (seconds)")+
  theme(panel.background = element_rect(color="black"))
ggsave(here("figures", "5_glmm_runtime_test.jpeg"))
