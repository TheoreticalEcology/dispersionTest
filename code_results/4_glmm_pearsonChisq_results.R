### Dispersion Tests Project
## Melina Leite
# Sep 26

library(DHARMa)
library(tidyverse)
library(here)
library(cowplot);
theme_set(theme_cowplot())
library(patchwork)


# plot Colors
source(here("functions_others", "plotColors.R"))
source(here("functions_others", "mcse.R"))



## Binomial and Poisson ####

load(here("data", "4_glmmBin_pearsonChisq.Rdata")) 
load(here("data", "4_glmmPois_pearsonChisq.Rdata"))

## Run Time ####
# in minutes
(tspent.bin <- map_dbl(out.bin, "time"))
(tspent.pois <- map_dbl(out.pois, "time"))
mean(c(tspent.bin,tspent.pois)) 


## simulations ####
simuls.bin <- map_dfr(out.bin, "simulations", .id="ngroups")  %>%
  rename("overdispersion" = "controlValues")

simuls.pois <- map_dfr(out.pois, "simulations", .id="ngroups")  %>%
  rename("overdispersion" = "controlValues",
         "Pear2.p.val" = "Pear2.p.var",
         "PearG.p.val" = "PearG.p.var")

## power ####
p.bin <- simuls.bin %>% dplyr::select(Pear2.p.val, PearG.p.val, replicate,
                                      ngroups,
                                      overdispersion) %>%
  pivot_longer(1:2, names_to = "test", values_to = "p.val") %>%
  mutate(ngroups = fct_relevel(ngroups, "10", "20", "50", "100")) %>%
  group_by(ngroups, overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]))
p.bin$prop.sig <- p.bin$p.sig/p.bin$nsim

for (i in 1:nrow(p.bin)) {
  btest <- binom.test(p.bin$p.sig[i], n=p.bin$nsim[i], p=0.05)
  p.bin$p.bin0.05[i] <- btest$p.value
  p.bin$conf.low[i] <- btest$conf.int[1]
  p.bin$conf.up[i] <- btest$conf.int[2]
}

p.pois <- simuls.pois %>% dplyr::select(Pear2.p.val, PearG.p.val, replicate,
                                      ngroups,
                                      overdispersion) %>%
  pivot_longer(1:2, names_to = "test", values_to = "p.val") %>%
  mutate(ngroups = fct_relevel(ngroups, "10", "20", "50", "100")) %>%
  group_by(ngroups, overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]))
p.pois$prop.sig <- p.pois$p.sig/p.pois$nsim

for (i in 1:nrow(p.pois)) {
  btest <- binom.test(p.pois$p.sig[i], n=p.pois$nsim[i], p=0.05)
  p.pois$p.bin0.05[i] <- btest$p.value
  p.pois$conf.low[i] <- btest$conf.int[1]
  p.pois$conf.up[i] <- btest$conf.int[2]
}

power <- bind_rows(list(Binomial = p.bin, Poisson=p.pois), .id="model") %>%
  mutate(model = fct_relevel(model, "Poisson", "Binomial")) %>%
  mutate(test = fct_recode(test, `two-way`="Pear2.p.val",
                           `greater`= "PearG.p.val"))



## dispersion statistics ####

st.bin <- simuls.bin %>% 
  mutate(ngroups = fct_relevel(ngroups, "10", "20", "50", "100")) %>%
  group_by(ngroups, overdispersion) %>%
  summarise(mean.stat = mean(Pear.stat.dispersion, na.rm=T),
            sd.stat = sd(Pear.stat.dispersion, na.rm=T),
            mcse = mcse_mean(Pear.stat.dispersion))

st.pois <- simuls.pois %>% 
  mutate(ngroups = fct_relevel(ngroups, "10", "20", "50", "100")) %>%
  group_by(ngroups, overdispersion) %>%
  summarise(mean.stat = mean(Pear.stat.dispersion, na.rm=T),
            sd.stat = sd(Pear.stat.dispersion, na.rm=T),
            mcse = mcse_mean(Pear.stat.dispersion))

disper <- bind_rows(list(Binomial = st.bin, Poisson=st.pois), .id="model")%>%
  mutate(model = fct_relevel(model, "Poisson", "Binomial"))



###### figure power #####


fig.power <- 
  power %>% filter(ngroups %in% c(10,50,100)) %>%
  ggplot( aes(x=overdispersion, y=prop.sig, col=model,
                               linetype=test))+
  geom_point(alpha=0.7, shape=1) + geom_line(alpha=0.7) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up),
                width = 0,
                linetype = "solid",
                show.legend = FALSE) +
  annotate("rect", xmin = -0.05, xmax = 0.05, ymin = 0, ymax = 1,
           alpha = .1,fill = "blue")+
  scale_linetype_discrete(name= "Test")+ scale_color_discrete("Model")+
  facet_wrap(~ngroups, labeller = as_labeller(c(`10`= "m = 10 groups",
                                                `50`= "m = 50 groups",
                                                `100`= "m = 100 groups",
                                                `Binomial` = "Binomial",
                                                `Poisson` = "Poisson"))) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  theme(panel.background = element_rect(color="black"),
    legend.position ="inside",
    legend.position.inside = c(0.23,0.5),
    legend.box.background = element_rect(color="gray94", fill="gray94")) +
  labs(tag="A)") +
  ylab("Power") + ylim(0,1)
fig.power




###### figure dispersion stat ####


fig.disp <- disper %>% filter(ngroups %in% c(10,50,100)) %>%
  ggplot(aes(x=overdispersion, y=mean.stat, col=model))+
  geom_point(alpha=0.7, shape=1) + geom_line( alpha=0.7) +
  geom_errorbar(aes(ymin=mean.stat-1.96*mcse, ymax=mean.stat+1.96*mcse),
                width = 0,
                linetype = "solid",
                show.legend = FALSE)+
  facet_grid(~ngroups, labeller = as_labeller(c(`10`= "m = 10 groups",
                                                `50`= "m = 50 groups",
                                                `100`= "m = 100 groups"))) +
  geom_hline(yintercept = 1, linetype="dotted", col="gray") +
  scale_color_discrete("Model")+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
           alpha = 0.1,fill = "red")+
  theme(panel.background = element_rect(color="black"),
        legend.position ="inside",
        legend.position.inside = c(0.018,0.83),
        legend.box.background = element_rect(color="gray94", fill="gray94"))+
        
  labs(tag="B)") +
  ylab("Dispersion statistics")+
  scale_y_log10()
fig.disp
fig.power + fig.disp + plot_layout(ncol=1)  +
  plot_annotation(title="Chi-squared Pearson dispersion tests for GLMMs",
                  theme = theme(plot.title = element_text(hjust=0.5)))


ggsave(here("figures", "4_glmm_pearsonChisq.pdf"), width=12, height = 8)
