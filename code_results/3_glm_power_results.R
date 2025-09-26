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

# plot Colors
source(here("code_results", "plotColors.R"))


##############-###
#### Binomial ####
##############-###

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

###### "RAW" Power ######

p.bin <- simuls.bin %>% dplyr::select(Pear.p.val,DHA.p.val,Ref.p.val, replicate,
                                      overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:3, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize,intercept,overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]) )
p.bin$prop.sig <- p.bin$p.sig/p.bin$nsim
p.bin$intercept <- fct_relevel(p.bin$intercept, "-3", "-1.5", "0", "1.5", "3")
p.bin$sampleSize <- as.factor(as.numeric(p.bin$sampleSize))

##### Callibrated Power #####

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


##### figure ####
bind_rows(list(uncalibrated=p.bin, calibrated= cp.bin), .id="model") %>%
  ggplot(aes(x=overdispersion, y=prop.sig, col=test, linetype = model))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_manual( values= col.tests[c(4,1,2)],
                      labels=c("Sim-based response variance", 
                               "param. Pearson residuals",
                               "nonparam. Pearson residuals"))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Binomial power", subtitle = "1000 sim; Ntrials = 10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom")+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave(here("figures", "3_glmBin_power.jpeg"), width=10, height = 15)



##### figure statistics #####

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
  scale_color_manual( values= col.tests[c(4,1,2)],
                      labels=c("Sim-based response variance", 
                               "param. Pearson residuals",
                               "nonparam. Pearson residuals"))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Binomial: dispersion statistics", subtitle = "1000 sim; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom")
ggsave(here("figures", "3_glmBin_dispersionStats.jpeg"), width=10, height = 15)


###############-##
#### Poisson  ####
##############-###


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


##### "RAW" Power #####
p.pois <- simuls.pois %>% dplyr::select(Pear.p.val,DHA.p.val,Ref.p.val, replicate,
                                        overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:3, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize,intercept,overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]) )
p.pois$prop.sig <- p.pois$p.sig/p.pois$nsim
p.pois$intercept <- fct_relevel(p.pois$intercept, "-3", "-1.5", "0", "1.5", "3")
p.pois$sampleSize <- as.factor(as.numeric(p.pois$sampleSize))

##### Calibrated power #####

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


##### figure ####
bind_rows(list(uncalibrated=p.pois, calibrated= cp.pois), .id="model") %>%
  ggplot(aes(x=overdispersion, y=prop.sig, col=test, linetype = model))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_color_manual( values= col.tests[c(4,1,2)],
                      labels=c("Sim-based response variance", 
                               "param. Pearson residuals",
                               "nonparam. Pearson residuals"))+
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Poisson power", subtitle = "1000 sim") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom")+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave(here("figures", "3_glmPois_power.jpeg"), width=10, height = 15)


##### figure statistics #####

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
  scale_color_manual( values= col.tests[c(4,1,2)],
                      labels=c("Sim-based response variance", 
                               "param. Pearson residuals",
                               "nonparam. Pearson residuals"))+
  facet_grid(sampleSize~intercept, scales="free_y") +
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Poisson: dispersion statistics", subtitle = "1000 simulations") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom")
ggsave(here("figures", "3_glmPois_dispersionStats.jpeg"), width=10, height = 15)





## FIGURE POWER together ####

# raw and callibrated
pow <- bind_rows(list(Poisson_uncalibrated = p.pois, 
                      Binomial_uncalibrated = p.bin,
                      Poisson_calibrated = cp.pois, 
                      Binomial_calibrated = cp.bin), .id= "model") %>%
  separate(model, c("model", "calibration")) %>%
  mutate(model = fct_relevel(model, "Poisson", "Binomial"))

# intercept = 0, some sample sizes
sub.pow <- pow %>% filter(intercept == 0, sampleSize %in% c(10,100,1000))

ggplot(sub.pow, aes(x=overdispersion, y=prop.sig, col=test, linetype = calibration)) +
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  ylab("Power") + xlab("Overdispersion") +
  scale_color_manual(values = col.tests[c(4,1,2)], 
                     labels = c("Sim-based response variance", 
                                "param. Pearson residuals",
                                 "nonparam. Pearson residuals")) +
  facet_grid(model~sampleSize, 
            labeller = as_labeller(c("Binomial"= "Binomial",
                                     "Poisson" = "Poisson",
                                     "10" = "n = 10",
                                     "100" = "n = 100",
                                     "1000" = "n = 1000"))) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  guides(
    color = guide_legend(position="inside"),
    linetype   = guide_legend(position ="inside")) +
  theme(panel.background = element_rect(color="black"),
        legend.spacing.y = unit(0, "cm"),
        legend.text = element_text(size=9),
        legend.key.spacing.y = unit(-0.1, "lines"),
        legend.title = element_blank(),
        legend.box.background = element_rect(fill = "gray94", color="gray94"),
        legend.position = c(0.01,0.87))

ggsave(here("figures", "3_glm_both_power.jpeg"), width=10, height = 6)



## FIGURE DISPERSION STATISTICS together ####

dispersion <- bind_rows(list(Poisson = st.pois, Binomial = st.bin), .id= "model") %>%
  pivot_wider(names_from = test, values_from = mean.stat) %>%
  mutate(reldif_DHA_Pear = (DHA.stat.dispersion - Pear.stat.dispersion)/DHA.stat.dispersion,
         reldif_DHA_Ref = (DHA.stat.dispersion - Ref.stat.dispersion)/DHA.stat.dispersion)


##### DHARMa stats X Pearson Chi-sq ####
small.text <- data.frame(label = c("Sim-based higher than Pearson", 
                                       "Sim-based lower than Pearson", "", ""),
                         model= rep(c("Binomial", "Poisson"), each=2),
                         x = rep(0.7,4), y=c(0.03,-0.025, 0,0 ))
#all results
ggplot(dispersion, aes(x=overdispersion, y=reldif_DHA_Pear, col=sampleSize)) +
  geom_point(size=2) + geom_line()+
  facet_grid(intercept~model, scales="free")+
  ylab("Relative diff. Dispersion stats \n (Sim-based x param. Pearson)")+
  geom_hline(yintercept = 0, linetype="dotted") +
  theme(panel.background = element_rect(color="black"))

###### figure to present intercept == 0 ####
dispersion %>% filter(intercept == 0) %>%
ggplot(aes(x=overdispersion, y=reldif_DHA_Pear, col=sampleSize)) +
  geom_point(size=2) + geom_line()+
  facet_grid(~model)+
  ylab("Relative diff. Dispersion stats \n (Sim-based x param. Pearson)")+
  geom_hline(yintercept = 0, linetype="dotted")  +
  ylim(-0.4,0.4) +
  geom_text(data=small.text, aes(x=x, y=y, label = label, group=model, 
                                 colour=""),
           colour="black") +
  theme(panel.background = element_rect(color="black"))
ggsave(here("figures", "3_glm_DISP_diff_DHA-Pear.jpeg"), height=4, width=9)



##### DHARMa stats X Pearson Bootstrapping ####
dispersion %>% filter(intercept == 0) %>%
ggplot(aes(x=overdispersion, y=reldif_DHA_Ref, col=sampleSize)) +
  geom_point(size=2) + geom_line()+
  facet_grid(~model)+
  ylab("Relative diff. Dispersion stats \n (Sim-based x nonparam. Pearson)")+
  ylim(-0.4,0.2)+
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_text(data=small.text, aes(x=x, y=y, label = label, group=model, colour=""),
            colour="black")+
  theme(panel.background = element_rect(color="black"))
ggsave(here("figures", "3_glm_DISP_diff_DHA-Ref.jpeg"), height=4, width=9)
