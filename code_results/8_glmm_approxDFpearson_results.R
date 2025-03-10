### Dispersion Tests Project
## Melina Leite
# Mar 25

library(DHARMa)
library(tidyverse)
library(here)
library(cowplot);
theme_set(theme_cowplot())
library(patchwork)

# plot Colors
source(here("code_results", "plotColors.R"))


load(here("data", "8_approximateDFpearson.Rdata")) # simulated data

(tspent <- map_dbl(out.pois, "time")) # seconds


#simulations
simuls.pois <- map_dfr(out.pois, "simulations", .id="model") %>%
  separate(model, c("ngroups", "sampleSize", "intercept"), sep="_") %>%
  rename("overdispersion" = "controlValues")

# power
p.pois <- simuls.pois %>% dplyr::select(ends_with(".p"), replicate,
                                        ngroups,
                                        overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:4, names_to = "test", values_to = "p.val") %>%
  group_by(sampleSize, ngroups, intercept, overdispersion, test) %>%
  summarise(p.sig = sum(p.val<0.05,na.rm=T),
            nsim = length(p.val[!is.na(p.val)]))
p.pois$prop.sig <- p.pois$p.sig/p.pois$nsim
p.pois$intercept <- fct_relevel(p.pois$intercept, "-1.5", "0", "1.5")
p.pois$ngroups <- fct_relevel(p.pois$ngroups, "10", "50", "100")
p.pois$sampleSize <- as.factor(as.numeric(p.pois$sampleSize))

###### figure power ####

p.pois %>%
  ggplot(aes(x=overdispersion, y=prop.sig, col=test, linetype=ngroups, shape=ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  ggtitle("Poisson", subtitle = "100 sim") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=2, byrow=TRUE))
#ggsave(here("figures", "5_glmmPois_powerALL.jpeg"), width=12, height = 15)
