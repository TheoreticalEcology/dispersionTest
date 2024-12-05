### Dispersion Tests Project
## Melina Leite
# Nov 24

library(DHARMa)
library(tidyverse); library(cowplot);
theme_set(theme_cowplot())
library(here)
library(patchwork)


####################
##### Binomial #####
####################

#load(here("data","2_glmBin_type1.Rdata"))

# prep data
# 
# simuls <- list()
# for (i in 1:length(out.out)) {
#   sim <- out.out[[i]]$simulations
#   sim$intercept <- names(out.out)[i]
#   simuls <- rbind(simuls,sim)
# }
# names(simuls)[names(simuls)=="controlValues"] <- "sampleSize"
# 
# simuls <- simuls %>% pivot_longer(1:3, names_to = "test", values_to = "p.val") %>%
#   group_by(sampleSize,intercept, test) %>%
#   summarise(p.sig = sum(p.val<0.05,na.rm=T),
#             nsim = length(p.val[!is.na(p.val)]) )
# simuls$prop.sig <- simuls$p.sig/simuls$nsim
# for (i in 1:nrow(simuls)) {
#   btest <- binom.test(simuls$p.sig[i], n=simuls$nsim[i], p=0.05)
#   simuls$p.bin0.05[i] <- btest$p.value
#   simuls$conf.low[i] <- btest$conf.int[1]
#   simuls$conf.up[i] <- btest$conf.int[2]
# }
# simuls$intercept <- as.factor(as.numeric(simuls$intercept))
# 
# simuls$test <- factor(simuls$test, levels = c("Pear.p.val", "Ref.p.val","DHA.p.val"))
# 
# 
# 
# ggplot(simuls, aes(y = prop.sig, x=as.factor(sampleSize), col=intercept)) +
#   facet_wrap(~test, labeller = as_labeller(c(`DHA.p.val`="Quantile residuals" ,
#                                              `Pear.p.val`="Pearson-Chisq" ,
#                                              `Ref.p.val`="Pearson Param. Bootstrapping"))) +
#   geom_point(position = position_dodge(width=0.8)) +
#   geom_errorbar(position = position_dodge(width=0.8),
#                 aes(ymin=conf.low, ymax=conf.up), width = 0.1)+
#   geom_hline(yintercept = 0.05, linetype="dotted")+
#   geom_line(aes(x=as.numeric(as.factor(sampleSize))),
#             position = position_dodge(width=0.8))+
#   ggtitle("GLM Binomial",
#           sub = "95% CIs with exact Binomial tests; 10000 simulations; Ntrials=10") +
#   xlab("Sample size") +
#   ylab("Type I error") +
#   theme(panel.background  = element_rect(color = "black"))
# ggsave(here("figures", "2_glmBin_type1.jpeg"), width=9, height = 5)

###################
##### Poisson #####
###################


#load(here("data", "2_glmPois_type1.Rdata"))

# simuls <- list()
# for (i in 1:length(out.out)) {
#   sim <- out.out[[i]]$simulations
#   sim$intercept <- names(out.out)[i]
#   simuls <- rbind(simuls,sim)
# }
# names(simuls)[names(simuls)=="controlValues"] <- "sampleSize"
# 
# simuls <- simuls %>% pivot_longer(1:3, names_to = "test", values_to = "p.val") %>%
#   group_by(sampleSize,intercept, test) %>%
#   summarise(p.sig = sum(p.val<0.05,na.rm=T),
#             nsim = length(p.val[!is.na(p.val)]) )
# simuls$prop.sig <- simuls$p.sig/simuls$nsim
# for (i in 1:nrow(simuls)) {
#   btest <- binom.test(simuls$p.sig[i], n=simuls$nsim[i], p=0.05)
#   simuls$p.bin0.05[i] <- btest$p.value
#   simuls$conf.low[i] <- btest$conf.int[1]
#   simuls$conf.up[i] <- btest$conf.int[2]
#   }
# simuls$intercept <- as.factor(as.numeric(simuls$intercept))
# 
# simuls$test <- factor(simuls$test, levels = c("Pear.p.val", "Ref.p.val","DHA.p.val"))
# 
# 
# #######################################
# ##### Figures and summary results #####
# #######################################
# 
# #### Type I error rate for the dispersion tests #####
# 
# 
# ggplot(simuls, aes(y = prop.sig, x=as.factor(sampleSize), col=intercept)) +
#   facet_wrap(~test, labeller = as_labeller(c(`DHA.p.val`="Quantile residuals" ,
#                                              `Pear.p.val`="Pearson-Chisq" ,
#                                              `Ref.p.val`="Pearson Param. Bootstrapping"))) +
#   geom_point(position = position_dodge(width=0.8)) +
#   geom_errorbar(position = position_dodge(width=0.8),
#                 aes(ymin=conf.low, ymax=conf.up), width = 0.1)+
#   geom_hline(yintercept = 0.05, linetype="dotted")+
#   geom_line(aes(x=as.numeric(as.factor(sampleSize))),
#             position = position_dodge(width=0.8))+
#   ggtitle("GLM Poisson",
#           sub = "95% CIs with exact Binomial tests; 10,000 simulations") +
#   xlab("Sample size") +
#   ylab("Type I error") +
#   scale_y_sqrt(breaks = c(0,0.01,0.05,0.2,0.4,0.6))+
#   theme(panel.background  = element_rect(color = "black"))
# ggsave(here("figures", "2_glmPois_type1.jpeg"), width=9, height = 5)


