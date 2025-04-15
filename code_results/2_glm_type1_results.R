### Dispersion Tests Project
## Melina Leite
# Nov 24

library(DHARMa)
library(tidyverse); library(cowplot);
theme_set(theme_cowplot())
library(here)
library(patchwork)

# plot Colors
source(here("code_results", "plotColors.R"))


##############-###
#### Binomial ####
##############-###

load(here("data","2_glmBin_type1.Rdata"))

# prep data

simuls.bin <- list()
for (i in 1:length(out.bin)) {
  sim <- out.bin[[i]]$simulations
  sim$intercept <- names(out.bin)[i]
  simuls.bin <- rbind(simuls.bin,sim)
}
names(simuls.bin)[names(simuls.bin)=="controlValues"] <- "sampleSize"



##### p values #####
p.bin <- simuls.bin %>% dplyr::select(Pear.p.val,DHA.p.val,Ref.p.val, replicate,
                                intercept, sampleSize) %>%
  pivot_longer(1:3, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize,intercept, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]) )
p.bin$prop.sig <- p.bin$p.sig/p.bin$nsim
for (i in 1:nrow(p.bin)) {
  btest <- binom.test(p.bin$p.sig[i], n=p.bin$nsim[i], p=0.05)
  p.bin$p.bin0.05[i] <- btest$p.value
  p.bin$conf.low[i] <- btest$conf.int[1]
  p.bin$conf.up[i] <- btest$conf.int[2]
}
p.bin$intercept <- as.factor(as.numeric(p.bin$intercept))
p.bin$test <- factor(p.bin$test, levels = c("Pear.p.val", "Ref.p.val","DHA.p.val"))



##### Dispersion statistics ####
stats.bin <- simuls.bin %>% dplyr::select(Pear.stat.dispersion,
                                          DHA.stat.dispersion,
                                          Ref.stat.dispersion, replicate,
                                         intercept, sampleSize) %>%
  pivot_longer(1:3, names_to = "test", values_to = "Dispersion") %>%
  mutate(intercept = as.factor(as.numeric(intercept)))



##### distribution p.values ####
pvals.bin <- simuls.bin %>% dplyr::select(Pear.p.val,DHA.p.val,Ref.p.val, replicate,
                                          intercept, sampleSize) %>%
  pivot_longer(1:3, names_to = "test", values_to = "p.val")  %>%
  mutate(intercept = as.factor(as.numeric(intercept)))



##### find the callibrated alpha for 0.05% p-value: ####
alpha.bin <- pvals.bin %>% group_by(sampleSize, intercept, test) %>%
  summarise(alpha = quantile(ecdf(p.val), 0.05))
# saving it together with alpha.pois to use in script 3_glm_power_results.R
# to callibrate the power



#############-###
#### Poisson ####
#############-###

load(here("data", "2_glmPois_type1.Rdata"))

simuls.pois <- list()
for (i in 1:length(out.pois)) {
  sim <- out.pois[[i]]$simulations
  sim$intercept <- names(out.pois)[i]
  simuls.pois <- rbind(simuls.pois,sim)
}
names(simuls.pois)[names(simuls.pois)=="controlValues"] <- "sampleSize"

p.pois <- simuls.pois %>% dplyr::select(Pear.p.val,DHA.p.val,Ref.p.val, replicate,
                                        intercept, sampleSize) %>%
  pivot_longer(1:3, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize,intercept, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]) )
p.pois$prop.sig <- p.pois$p.sig/p.pois$nsim
for (i in 1:nrow(p.pois)) {
  btest <- binom.test(p.pois$p.sig[i], n=p.pois$nsim[i], p=0.05)
  p.pois$p.bin0.05[i] <- btest$p.value
  p.pois$conf.low[i] <- btest$conf.int[1]
  p.pois$conf.up[i] <- btest$conf.int[2]
}
p.pois$intercept <- as.factor(as.numeric(p.pois$intercept))

p.pois$test <- factor(p.pois$test, levels = c("Pear.p.val", "Ref.p.val","DHA.p.val"))

##### Dispersion statistics ####
stats.pois <- simuls.pois %>% dplyr::select(Pear.stat.dispersion,
                                          DHA.stat.dispersion,
                                          Ref.stat.dispersion, replicate,
                                          intercept, sampleSize) %>%
  pivot_longer(1:3, names_to = "test", values_to = "Dispersion") %>%
  mutate(intercept = as.factor(as.numeric(intercept)))


##### distribution p.values ####

pvals.pois <- simuls.pois %>% dplyr::select(Pear.p.val,DHA.p.val,Ref.p.val, replicate,
                                            intercept, sampleSize) %>%
  pivot_longer(1:3, names_to = "test", values_to = "p.val")  %>%
  mutate(intercept = as.factor(as.numeric(intercept)))


##### find the callibrated alpha for 0.05 p-value: ####

alpha.pois <- pvals.pois %>% group_by(sampleSize, intercept, test) %>%
  summarise(alpha = quantile(ecdf(p.val), 0.05))
alpha.pois
save(alpha.pois, alpha.bin, file=here("data", "2_callibrated_alphaLevels.Rdata"))



##################################-###
#### Figures and summary results #####
##################################-###


##### Type I error rate for the dispersion tests #####

f.bin <- ggplot(p.bin, aes(y = prop.sig, x=as.factor(sampleSize), 
                           col=intercept)) +
  facet_wrap(~test, labeller = 
               as_labeller(c(`DHA.p.val` = "2) Sim-based residuals" ,
                            `Pear.p.val` = "1a) Pearson Chi-squared" ,
                           `Ref.p.val` = "1b) Pearson param. bootstrap"))) +
  geom_point(position = position_dodge(width=0.8)) +
  geom_errorbar(position = position_dodge(width=0.8), col = "black",
                aes(ymin=conf.low, ymax=conf.up, group=intercept), width = 0.1)+
  geom_hline(yintercept = 0.05, linetype="dotted")+
  geom_line(aes(x=as.numeric(as.factor(sampleSize))),
            position = position_dodge(width=0.8))+
  scale_color_manual(values = col.intercept)+
  ggtitle("GLM Binomial") +
  xlab("Sample size") +
  ylab("Type I error") +
  theme(panel.background  = element_rect(color = "black"),
        legend.position = c(0.9,0.75),
        axis.text.x = element_text(angle=45, hjust=1))+
  labs(tag="B)")
f.bin
#ggsave(here("figures", "2_glmBin_type1.jpeg"), width=10, height = 5)

f.pois <- ggplot(p.pois, aes(y = prop.sig, x=as.factor(sampleSize), 
                             col=intercept)) +
  facet_wrap(~test, labeller = as_labeller(c(`DHA.p.val` = "2) Sim-based residuals" ,
                                             `Pear.p.val` = "1a) Pearson Chi-squared" ,
                                    `Ref.p.val` = "1b) Pearson param. bootstrap"))) +
  geom_point(position = position_dodge(width=0.8)) +
  geom_errorbar(position = position_dodge(width=0.8), col="black",
                aes(ymin=conf.low, ymax=conf.up, group=intercept), width = 0.1)+
  geom_hline(yintercept = 0.05, linetype="dotted")+
  geom_line(aes(x=as.numeric(as.factor(sampleSize))),
            position = position_dodge(width=0.8))+
  scale_color_manual(values = col.intercept)+
  ggtitle("GLM Poisson") +
  xlab("Sample size") +
  ylab("Type I error") +
  scale_y_sqrt(breaks = c(0,0.01,0.05,0.2,0.4,0.6))+
  theme(panel.background  = element_rect(color = "black"),
        legend.position = c(0.9,0.75),
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(tag="A)")
f.pois
#ggsave(here("figures", "2_glmPois_type1.jpeg"), width=10, height = 5)


## both distributions

dats <- bind_rows(list(Poisson = p.pois, Binomial = p.bin), .id="model") %>%
  mutate(model = fct_relevel(model, "Poisson", "Binomial")) 

dats %>%
  ggplot(aes(y = prop.sig, x=as.factor(sampleSize), col=intercept,
             group=intercept)) +
  facet_grid(model~test, scales="free",
             labeller = as_labeller(c(`DHA.p.val` = "2) Sim-based dispersion" ,
                                       `Pear.p.val` = "1a) Pearson Chi-squared",
                                `Ref.p.val` = "1b) Pearson param. bootstrap",
                                         `Binomial` = "Binomial",
                                      `Poisson` = "Poisson"))) +
  scale_y_sqrt(breaks = c(0,0.01,0.05,0.2,0.4,0.6))+
  geom_point(position = position_dodge(width=0.8), col="white") +
  geom_hline(yintercept = 0.05, linetype="dotted")+
  geom_hline(yintercept = 0, col="gray")+
  geom_line(aes(x=as.numeric(as.factor(sampleSize))),
            position = position_dodge(width=0.8))+
  geom_errorbar(position = position_dodge(width=0.8), #col="black",
                aes(ymin=conf.low, ymax=conf.up, group=intercept), width = 0.1)+
  scale_color_manual(values = col.intercept)+
  xlab("Sample size") +
  ylab("Type I error") +
  theme(panel.background  = element_rect(color = "black"),
        legend.position = c(0.9,0.86),
        legend.background  = element_rect(fill="#F0F0F0"),
        axis.text.x = element_text(angle=45, hjust=1))

ggsave(here("figures", "2_glm_type1.jpeg"), width=10, height = 7)



##### Dispersion statistics #### 

# summary stats

d.bin <- stats.bin %>% group_by(test, sampleSize, intercept) %>%
  summarise(mean = mean(Dispersion,na.rm=T),
            median = median(Dispersion, na.rm=T),
            sd = sd(Dispersion, na.rm=T)) %>%
 # pivot_longer(c(mean,median), names_to = "stat", values_to = "Dispersion") %>%
  ggplot(aes(x=sampleSize, y=median, col=test))+
  geom_point()+ geom_line()+
  scale_x_log10()+
  facet_grid(~intercept) +
  scale_color_manual(values=col.tests[c(4,1,2)],
    labels=c("Sim-based dispersion", "Pearson Chi-squared",
             "Pearson Param. Bootstrap."))+
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Binomial") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom",
        axis.text.x = element_text(angle=45, hjust=1))+
  ylab("Dispersion")+
  labs(tag="B)")
d.bin

d.pois <- stats.pois %>% group_by(test,sampleSize, intercept) %>%
  summarise(mean = mean(Dispersion,na.rm=T),
            median = median(Dispersion, na.rm=T),
            sd = sd(Dispersion, na.rm=T)) %>%
  #pivot_longer(c(mean,median), names_to = "stat", values_to = "Dispersion") %>%
  ggplot(aes(x=sampleSize, y=median, col=test))+
  geom_point()+ geom_line()+
  #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))+
  scale_x_log10()+
  facet_grid(~intercept) +
  scale_color_manual(values=col.tests[c(4,1,2)],
                     labels=c("Sim-based dispersion", "Pearson Chi-squared",
                              "Pearson Param. Bootstrap."))+
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Poisson") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1)) +
  ylab("Dispersion")+
  labs(tag="A)")
d.pois

d.pois + d.bin +
  plot_layout(ncol=1)

ggsave(here("figures", "2_glm_dispersionStats.jpeg"), width=10, height = 8)


##### distribution p values ####

ggplot(pvals.bin, aes(x=p.val, col=test))+
  geom_density()+
  geom_hline(yintercept=1, linetype="dotted")+
  facet_grid(sampleSize ~ intercept, scales="free") +
  scale_color_manual(values=col.tests[c(4,1,2)],
                     labels=c("Sim-based dispersion", "Pearson Chi-squared",
                              "Pearson Param. Bootstrap."))+
  theme(panel.background = element_rect(color="black"),
        axis.text.x = element_text(angle=45,  hjust=1),
        legend.position = "bottom") +
  ggtitle("Binomial: distribution P-values")
ggsave(here("figures", "2_glmBin_distribPvals.jpeg"), width=10, height = 15)


ggplot(pvals.pois, aes(x=p.val, col=test))+
  geom_density()+
  geom_hline(yintercept=1, linetype="dotted")+
  facet_grid(sampleSize ~ intercept, scales="free") +
  scale_color_manual(values=col.tests[c(4,1,2)],
                     labels=c("Sim-based dispersion", "Pearson Chi-squared",
                              "Pearson Param. Bootstrap."))+
  theme(panel.background = element_rect(color="black"),
        axis.text.x = element_text(angle=45,  hjust=1),
        legend.position = "bottom") +
  ggtitle("Poisson: distribution P-values")
ggsave(here("figures", "2_glmPois_distribPvals.jpeg"), width=10, height = 15)


##### ecdf p values ####

ggplot(pvals.bin, aes(x=p.val, col=test))+
  stat_ecdf()+
  #xlim(0,0.25) + ylim(0,0.2)+
  facet_grid(sampleSize ~ intercept) +
  geom_hline(yintercept=0.05, linetype="dotted")+
  theme(panel.background = element_rect(color="black"),
        axis.text.x = element_text(angle=45,  hjust=1),
        legend.position = "bottom")+
  ggtitle("Binomial: ecdf p-values")
#ggsave(here("figures", "2_glmBin_ecdfPvals.jpeg"), width=10, height = 15)

ggplot(pvals.pois, aes(x=p.val, col=test))+
  stat_ecdf()+
  #xlim(0,0.25) + ylim(0,0.2)+
  facet_grid(sampleSize ~ intercept) +
  geom_hline(yintercept=0.05, linetype="dotted")+
  theme(panel.background = element_rect(color="black"),
        axis.text.x = element_text(angle=45,  hjust=1),
        legend.position = "bottom")+
  ggtitle("Poisson: ecdf p-values")
#ggsave(here("figures", "2_glmPois_ecdfPvals.jpeg"), width=10, height = 15)

