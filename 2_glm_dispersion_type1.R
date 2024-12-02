### Dispersion Tests Project
## Melina Leite
# Nov 24

library(DHARMa)
library(tidyverse); library(cowplot);
theme_set(theme_cowplot())
library(here)
library(patchwork)

#####################################
##### Instructions & Simulating #####
#####################################

# 1) Simulating 1000 Poisson datasets with different sample sizes and intercepts
#       - sampleSize: c(10,50,100,500)
#       - intercept:  c(-3,-1,0,2,4)
# 2) fitting them to correct GLM models
# 5) calculating type I error rate for 
#     - Pearson-chisq, 
#     - DHARMa default (simulated residuals)
#     - DHARMa refit (boostrapped Pearson residuals)


# function to varying sampleSize
calculateStatistics <- function(control = 10){
  # data
  testData <- DHARMa::createData(sampleSize = control,
                                 numGroups = 1,
                                 family = poisson())
  # model
  fittedModel <- stats::glm(observedResponse ~ Environment1, 
                            data = testData, family = poisson()) 
  #results
  out <- list()
  
  # pearson residual
  out$Pear.p.val <- testDispersion(fittedModel, plot = F, 
                                   type="PearsonChisq")$p.value
  # DHARMa default residuals
  res <- simulateResiduals(fittedModel)
  out$DHA.p.val<- testDispersion(res, type = "DHARMa",plot = F)$p.value
  
  # DHARMa refit residuals -> bootstrapped Pearson
  res <- simulateResiduals(fittedModel, refit=T)
  out$Ref.p.val <- testDispersion(res, plot = F, type = "DHARMa")$p.value
  return(unlist(out))
}

# varying intercept in a loop 
intercept <- c(-3,-1,0,1,4)
sampleSize = c(10,50,100,500)

out.out <- list()
for (i in intercept){

  out <- runBenchmarks(calculateStatistics, controlValues = sampleSize,
                       intercept = i,
                       nRep=10000, parallel = T)
  out.out[[length(out.out) + 1]] <- out
}

names(out.out) <- intercept

# saving sim results
save(out.out,sampleSize,intercept, file=here("data", 
                                                "2_glm_dispersion_type1.Rdata"))

#load(here("data", "2_glm_dispersion_type1.Rdata"))


# prep data
# using the summary table doesn't work because there are NAs in the DHARMa results

# prop.sig <- list()
# for (i in 1:length(out.out)) {
#   prop.sig <- rbind(prop.sig, out.out[[i]]$summaries$propSignificant)
# }
# prop.sig$sampleSize = prop.sig$controlValues
# prop.sig$intercept = rep(intercept, each=length(unique(sampleSize)))
# 
# prop.sig <- prop.sig %>% pivot_longer(c(Pear.p.val,DHA.p.val, Ref.p.val), 
#                                       names_to = "tests", values_to = "p.val")
# prop.sig$binom.test <- binom.test(prop.sig$p.val,n = 1000,p=0.05)



simuls <- list()
for (i in 1:length(out.out)) {
  sim <- out.out[[i]]$simulations
  sim$intercept <- names(out.out)[i]
  simuls <- rbind(simuls,sim)
}
names(simuls)[names(simuls)=="controlValues"] <- "sampleSize"

simuls <- simuls %>% pivot_longer(1:3, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize,intercept, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]) )
simuls$prop.sig <- simuls$p.sig/simuls$nsim
for (i in 1:nrow(simuls)) {
  btest <- binom.test(simuls$p.sig[i], n=simuls$nsim[i], p=0.05)
  simuls$p.bin0.05[i] <- btest$p.value 
  simuls$conf.low[i] <- btest$conf.int[1]
  simuls$conf.up[i] <- btest$conf.int[2]
  }
simuls$intercept <- as.factor(as.numeric(simuls$intercept))

#######################################
##### Figures and summary results #####
#######################################

#### Type I error rate for the dispersion tests #####

simuls$asterisk <- paste0(round(simuls$prop.sig,3),ifelse(simuls$p.bin0.05 < 0.05,"*",""))


ggplot(simuls, aes(y=rev(as.factor(sampleSize)), x=as.factor(intercept), 
                   fill=prop.sig)) +
  facet_grid(rows = vars(test)) +
  geom_tile() +
  scale_y_discrete("sampleSize", labels=c(500,100,50,20,10)) +
  scale_x_discrete("intercept", position = "top") +
  geom_text(aes(label=asterisk)) +
  scale_fill_gradient2(name="type I", low=4, high=2, midpoint=0.05,
                       mid="white") +
  ggtitle("Type I error rates") 





ggplot(simuls, aes(y = prop.sig, x=intercept, col=test)) +
  facet_wrap(~sampleSize) +
  geom_point(position = position_dodge(width=0.8)) +
  geom_errorbar(position = position_dodge(width=0.8),
    aes(ymin=conf.low, ymax=conf.up), width = 0.1)+
  geom_line(aes(x=as.numeric(intercept)),
            position = position_dodge(width=0.8))+
  geom_hline(yintercept = 0.05, linetype="dotted") 



ggplot(simuls, aes(y = prop.sig, x=as.factor(sampleSize), col=intercept)) +
  facet_wrap(~test, labeller = as_labeller(c(`DHA.p.val`="DHARMa default" ,
                                             `Pear.p.val`="Pearson-Chisq" ,
                                             `Ref.p.val`="DHARMa refit"))) +
  geom_point(position = position_dodge(width=0.8)) +
  geom_errorbar(position = position_dodge(width=0.8),
                aes(ymin=conf.low, ymax=conf.up), width = 0.1)+
  geom_hline(yintercept = 0.05, linetype="dotted")+
  ggtitle("Dispersion tests: Type I error rates for GLM Poisson",
          sub = "95% CIs with exact Binomial tests; 10000 simulations") +
  xlab("Sample size") +
  ylab("Type I error") +
  theme(panel.background  = element_rect(color = "black"))
ggsave(here("figures", "2_glm_dispersion_type1.jpeg"), width=10)






















