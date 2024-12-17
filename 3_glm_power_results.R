### Dispersion Tests Project
## Melina Leite
# Dez 24

library(DHARMa)
library(tidyverse)
library(here)
library(cowplot);
theme_set(theme_cowplot())
library(patchwork)

load(here("data", "2_callibrated_alphaLevels.Rdata")) # callibrated alpha level
# created in script 2_glm_type1_results.R

##################-###
###### Binomial ######
##################-###

load(here("data", "3_glmBin_power.Rdata")) # simulated data


simuls.bin <- list()
for (i in 1:length(out.bin)) {
  sim <- out.bin[[i]]$simulations
  params <- strsplit(names(out.bin)[[i]], "_")[[1]]
  sim$intercept <- params[1]
  sim$sampleSize <- params[2]
  simuls.bin <- rbind(simuls.bin,sim)
}
names(simuls.bin)[names(simuls.bin)=="controlValues"] <- "overdispersion"

# "RAW" Power ######

p.bin <- simuls.bin %>% dplyr::select(Pear.p.val,DHA.p.val,Ref.p.val, replicate,
                                      overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:3, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize,intercept,overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]) )
p.bin$prop.sig <- p.bin$p.sig/p.bin$nsim
p.bin$intercept <- fct_relevel(p.bin$intercept, "-3", "-1.5", "0", "1.5", "3")
p.bin$sampleSize <- as.factor(as.numeric(p.bin$sampleSize))


ggplot(p.bin, aes(x=overdispersion, y=prop.sig, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Quantile Residuals", "Pearson Chi-squared",
             "Pearson Param. Bootstrap."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Binomial", subtitle = "1000 sim; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom")
ggsave(here("figures", "3_glmBin_power.jpeg"), width=10, height = 15)




## Callibrated Power #####

cp.bin <- simuls.bin %>% dplyr::select(Pear.p.val,DHA.p.val,Ref.p.val, replicate,
                             overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:3, names_to = "test", values_to = "p.val") %>%
  mutate(sampleSize = as.numeric(sampleSize)) %>%
  left_join(alpha.bin, by=c("sampleSize", "intercept", "test")) %>%
  mutate(significance = p.val <alpha) %>%
  group_by(sampleSize,intercept,overdispersion, test) %>%
  summarise(p.sig = sum(significance, na.rm = T),
            nsim = length(p.val[!is.na(p.val)]) )
cp.bin$prop.sig <- cp.bin$p.sig/cp.bin$nsim
cp.bin$intercept <- fct_relevel(cp.bin$intercept, "-3", "-1.5", "0", "1.5", "3")
cp.bin$sampleSize <- as.factor(as.numeric(cp.bin$sampleSize))


ggplot(cp.bin, aes(x=overdispersion, y=prop.sig, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Quantile Residuals", "Pearson Chi-squared",
             "Pearson Param. Bootstrap."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Binomial: callibrated power", subtitle = "1000 sim; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom")
ggsave(here("figures", "3_glmBin_powerCalibrated.jpeg"), width=10, height = 15)




# figure statistics #####

st.bin <- simuls.bin %>% dplyr::select(Pear.stat.dispersion,DHA.stat.dispersion,
                                       Ref.stat.dispersion, replicate,
                                            overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:3, names_to = "test", values_to = "Disp.stats")  %>%
  group_by(sampleSize,intercept,overdispersion, test) %>%
    summarise(mean.stat = mean(Disp.stats, na.rm=T))
st.bin$intercept <- fct_relevel(st.bin$intercept, "-3", "-1.5", "0", "1.5", "3")
st.bin$sampleSize <- as.factor(as.numeric(st.bin$sampleSize))


ggplot(st.bin, aes(x=overdispersion, y=mean.stat, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Quantile Residuals", "Pearson Chi-squared",
             "Pearson Param. Bootstrap."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Binomial: dispersion statistics", subtitle = "1000 sim; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom")
ggsave(here("figures", "3_glmBin_dispersionStats.jpeg"), width=10, height = 15)


#################-##
##### Poisson  #####
################-###


load(here("data", "3_glmPois_power.Rdata"))

simuls.pois <- list()
for (i in 1:length(out.pois)) {
  sim <- out.pois[[i]]$simulations
  params <- strsplit(names(out.pois)[[i]], "_")[[1]]
  sim$intercept <- params[1]
  sim$sampleSize <- params[2]
  simuls.pois <- rbind(simuls.pois,sim)
}
names(simuls.pois)[names(simuls.pois)=="controlValues"] <- "overdispersion"


# "RAW" Power #####
p.pois <- simuls.pois %>% dplyr::select(Pear.p.val,DHA.p.val,Ref.p.val, replicate,
                                        overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:3, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize,intercept,overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]) )
p.pois$prop.sig <- p.pois$p.sig/p.pois$nsim
p.pois$intercept <- fct_relevel(p.pois$intercept, "-3", "-1.5", "0", "1.5", "3")
p.pois$sampleSize <- as.factor(as.numeric(p.pois$sampleSize))

ggplot(p.pois, aes(x=overdispersion, y=prop.sig, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Quantile Residuals", "Pearson Chi-squared",
             "Pearson Param. Bootstrap."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Poisson", subtitle = "1000 sim; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom")
ggsave(here("figures", "3_glmPois_power.jpeg"), width=10, height = 15)


# Calibrated power #####

cp.pois <- simuls.pois %>% dplyr::select(Pear.p.val,DHA.p.val,Ref.p.val, replicate,
                                       overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:3, names_to = "test", values_to = "p.val") %>%
  mutate(sampleSize = as.numeric(sampleSize)) %>%
  left_join(alpha.pois, by=c("sampleSize", "intercept", "test")) %>%
  mutate(significance = p.val <alpha) %>%
  group_by(sampleSize,intercept,overdispersion, test) %>%
  summarise(p.sig = sum(significance, na.rm = T),
            nsim = length(p.val[!is.na(p.val)]) )
cp.pois$prop.sig <- cp.pois$p.sig/cp.pois$nsim
cp.pois$intercept <- fct_relevel(cp.pois$intercept, "-3", "-1.5", "0", "1.5", "3")
cp.pois$sampleSize <- as.factor(as.numeric(cp.pois$sampleSize))



ggplot(cp.pois, aes(x=overdispersion, y=prop.sig, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_discrete(
    labels=c("Quantile Residuals", "Pearson Chi-squared",
             "Pearson Param. Bootstrap."))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Poisson: callibrated power", subtitle = "1000 sim; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom")
ggsave(here("figures", "3_glmPois_powerCalibrated.jpeg"), width=10, height = 15)



### figure statistics #####

st.pois <- simuls.pois %>% dplyr::select(Pear.stat.dispersion,DHA.stat.dispersion,
                                         Ref.stat.dispersion, replicate,
                                         overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:3, names_to = "test", values_to = "Disp.stats")  %>%
  group_by(sampleSize,intercept,overdispersion, test) %>%
  summarise(mean.stat = mean(Disp.stats, na.rm=T))
st.pois$intercept <- fct_relevel(st.pois$intercept, "-3", "-1.5", "0", "1.5", "3")
st.pois$sampleSize <- as.factor(as.numeric(st.pois$sampleSize))


ggplot(st.pois, aes(x=overdispersion, y=mean.stat, col=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_y_log10()+
  scale_color_discrete(
    labels=c("Quantile Residuals", "Pearson Chi-squared",
             "Pearson Param. Bootstrap."))+
  facet_grid(sampleSize~intercept, scales="free_y") +
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Poisson: dispersion statistics", subtitle = "1000 sim; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom")
ggsave(here("figures", "3_glmPois_dispersionStats.jpeg"), width=10, height = 15)


# FIGURE DISPERSION STATISTICS together ####

dispersion <- bind_rows(list(Poisson = st.pois, Binomial = st.bin), .id= "model") %>%
  filter(intercept == 0, sampleSize %in% c(10,20,50,100)) %>%
  pivot_wider(names_from = test, values_from = mean.stat) %>%
  mutate(reldif_DHA_Pear = (DHA.stat.dispersion - Pear.stat.dispersion)/DHA.stat.dispersion,
         reldif_DHA_REf = (DHA.stat.dispersion - Ref.stat.dispersion)/DHA.stat.dispersion)


ggplot(dispersion, aes(x=overdispersion, y=reldif_DHA_Pear, col=sampleSize)) +
  geom_point(size=2) + geom_line()+
  facet_grid(~model)+
  ylab("Relative diff. Dispersion stats \n (Sim-based x Pearson)")+
  ylim(-0.4,0.2)+
  geom_hline(yintercept = 0, linetype="dotted") 
ggsave(here("figures", "3_glm_DISPERSION_diff.jpeg"), height=5, width=10)


