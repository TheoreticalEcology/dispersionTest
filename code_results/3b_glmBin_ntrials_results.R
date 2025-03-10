### Dispersion Tests Project
## Melina Leite
# Mar 25

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

load(here("data","3b_glmBin_ntrials.Rdata"))

#simulations
bin  <- map_dfr(trial.bin, "simulations", .id="ntrial") %>%  
  rename("overdispersion" = "controlValues")

# add ntrials 10
load(here("data", "3_glmBin_power.Rdata"))
bin10  <- out.bin$`0_500`$simulations %>%
  mutate(ntrial = "10") %>%
  rename("overdispersion" = "controlValues")

bins <- bind_rows(bin, bin10) %>%
  mutate(ntrial = fct_relevel(ntrial, "5", "10", "20"))


##### p values #####
p.bin <- bins %>% dplyr::select(Pear.p.val,DHA.p.val,Ref.p.val, replicate,
                               overdispersion, ntrial) %>%
  pivot_longer(1:3, names_to = "test", values_to = "p.val") %>%
  group_by(ntrial, overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]) )
p.bin$prop.sig <- p.bin$p.sig/p.bin$nsim
for (i in 1:nrow(p.bin)) {
  btest <- binom.test(p.bin$p.sig[i], n=p.bin$nsim[i], p=0.05)
  p.bin$p.bin0.05[i] <- btest$p.value
  p.bin$conf.low[i] <- btest$conf.int[1]
  p.bin$conf.up[i] <- btest$conf.int[2]
}
p.bin$test <- factor(p.bin$test, levels = c("Pear.p.val", "Ref.p.val","DHA.p.val"))

#power
power <- p.bin %>%
  ggplot(aes(x=overdispersion, y=prop.sig, col=test)) +
  geom_point() + geom_line()+
  scale_color_manual(values = col.tests[c(1,2,4)],
                     labels = c("Pearson Chi-squared",
                                "Pearson Param. Bootstrap.",
                                "Sim-based dispersion")) +
  facet_grid(~ntrial, labeller = as_labeller(c("10" = "10 trials",
                                               "5" = "5 trials",
                                               "20" = "20 trials"))) +
  ylab("Power")+
  theme(panel.background = element_rect(color="black"),
        legend.position = "none") +
  geom_hline(yintercept = 0.05, linetype="dotted")
power



#type 1
type1 <- p.bin %>% filter(overdispersion == 0) %>%
  ggplot(aes(x=ntrial, y=prop.sig, col=test))+
  geom_point(position = position_dodge(width = 0.2)) + 
  geom_line(aes(x=as.numeric(ntrial)),
            position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0.05, linetype="dashed") +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up), width=0.05,
                position = position_dodge(width = 0.2)) +
 scale_color_manual(values = col.tests[c(1,2,4)],
                   labels = c("Pearson Chi-squared",
                              "Pearson Param. Bootstrap.",
                              "Sim-based dispersion"))  +
  ylab("Type I error") + xlab("Number of trials")+
  theme(panel.background = element_rect(color="black"),
        legend.position = "right")
type1


##### Dispersion statistics ####
stats.bin <- bins %>% dplyr::select(Pear.stat.dispersion,
                                   DHA.stat.dispersion,
                                   Ref.stat.dispersion, 
                                   overdispersion, ntrial) %>%
  pivot_longer(1:3, names_to = "test", values_to = "Dispersion") %>%
  group_by(ntrial, overdispersion,test) %>% 
  summarise(mean.stat = mean(Dispersion, na.rm=T))
stats.bin$test <- factor(stats.bin$test, levels = c("Pear.stat.dispersion", 
                                                    "Ref.stat.dispersion",
                                                    "DHA.stat.dispersion"))


disp <- stats.bin %>%
  ggplot(aes(x=overdispersion, y=mean.stat, col=test)) +
  geom_point() + geom_line()+
  scale_color_manual(values = col.tests[c(1,2,4)],
                     labels = c("Pearson Chi-squared",
                                "Pearson Param. Bootstrap.",
                                "Sim-based dispersion")) +
  facet_grid(~ntrial, labeller = as_labeller(c("10" = "10 trials",
                                               "5" = "5 trials",
                                               "20" = "20 trials"))) +
  ylab("Dispersion statistic")+
  theme(panel.background = element_rect(color="black"),
        legend.position = "none")
disp

# Figures


(power /disp) + 
  (type1 + plot_spacer()) +
  plot_layout(ncol = 1)
ggsave(here("figures", "3b_glmBin_ntrials.jpeg"))


