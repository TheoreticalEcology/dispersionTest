### Dispersion Tests Project
## Melina Leite
# Dez 24

library(DHARMa)
library(tidyverse)
library(here)
library(cowplot);
theme_set(theme_cowplot())
library(patchwork)



## TESTS WITH TYPE I ERROR FROM THE POWER SIMULATIONS

load(here("data", "5_glmmBin_power.Rdata")) # simulated data
#simulations
simuls.bin <- map_dfr(out.bin, "simulations", .id="model") %>%
  separate(model, c("ngroups", "sampleSize", "intercept"), sep="_") %>%
  rename("overdispersion" = "controlValues")

# power
p.bin <- simuls.bin %>% filter(overdispersion == 0) %>%
  dplyr::select(Pear.p.val, dhaUN.p.val, dhaCO.p.val,
                                      refUN.p.val, refCO.p.val, replicate,
                                      ngroups,
                                      overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:5, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize, ngroups, intercept, overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]))
p.bin$prop.sig <- p.bin$p.sig/p.bin$nsim
p.bin$intercept <- fct_relevel(p.bin$intercept, "-3", "-1.5", "0", "1.5", "3")
p.bin$sampleSize <- as.factor(as.numeric(p.bin$sampleSize))
for (i in 1:nrow(p.bin)) {
  btest <- binom.test(p.bin$p.sig[i], n=p.bin$nsim[i], p=0.05)
  p.bin$p.bin0.05[i] <- btest$p.value
  p.bin$conf.low[i] <- btest$conf.int[1]
  p.bin$conf.up[i] <- btest$conf.int[2]
}
p.bin$test <- fct_relevel(p.bin$test, "Pear.p.val", "refUN.p.val", "refCO.p.val",
                          "dhaUN.p.val", "dhaCO.p.val")

ggplot(p.bin, aes(y = prop.sig, x=as.factor(sampleSize), col=intercept)) +
  facet_grid(~test, labeller = as_labeller(c(`dhaUN.p.val`="Sim-based unconditional",
                                             `dhaCO.p.val`="Sim-based conditional",
                                             `Pear.p.val`="Pearson-Chisq" ,
                                             `refUN.p.val`= "Pearson ParBoot. unconditional",
                                             `refCO.p.val`= "Pearson ParBoot. conditional"))) +
  geom_point(position = position_dodge(width=0.8)) +
  geom_errorbar(position = position_dodge(width=0.8),
                aes(ymin=conf.low, ymax=conf.up), width = 0.1)+
  geom_hline(yintercept = 0.05, linetype="dotted")+
  geom_line(aes(x=as.numeric(as.factor(sampleSize))),
            position = position_dodge(width=0.8))+
  ggtitle("GLM Binomial", subtitle="100 sim only") +
  xlab("Sample size") +
  ylab("Type I error") +
  theme(panel.background  = element_rect(color = "black"),
        legend.position = "bottom",
        axis.text.x = element_text(angle=45, hjust=1))+
  labs(tag="B)")
