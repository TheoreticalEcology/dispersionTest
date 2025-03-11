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

(tspent <- map_dbl(out.pois, "time")) # minutes


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
  geom_hline(yintercept = c(0.05,0.5), linetype="dotted") +
  ylim(0,1)+
  ggtitle("Poisson", subtitle = "1000 sim") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=2, byrow=TRUE))
#ggsave(here("figures", "5_glmmPois_powerALL.jpeg"), width=12, height = 15)


p.pois %>% filter(ngroups == "100") %>%
  ggplot(aes(x=overdispersion, y=prop.sig, col=test, linetype=ngroups, shape=ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = c(0.05,0.5), linetype="dotted") +
  ylim(0,1)+
  ggtitle("Poisson", subtitle = "1000 sim") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=2, byrow=TRUE))



###### Type 1 error #####

p.pois %>% filter(overdispersion == 0) %>% ungroup() %>%
  mutate(ngroups = fct_relevel(ngroups, "10", "50", "100")) %>%
  ggplot(aes(x=sampleSize, y=prop.sig, col=intercept))+
  geom_point( position = position_dodge(width = 0.3))+
  geom_hline(yintercept = 0.05, linetype="dotted")+
  geom_line(aes(x=as.numeric(sampleSize)),
            position = position_dodge(width = 0.3))+
  facet_grid(ngroups~test)+
  theme(panel.background = element_rect(color="black"),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position = c(0.01,0.9)) 
#ggsave(here("figures", "5_glmmPois_type1.jpeg"), width=12, height = 8)






###### dispersion stat ####
d.pois <- simuls.pois %>% dplyr::select(ends_with("dispersion"), replicate, ngroups,
                                        overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:4, names_to = "test", values_to = "dispersion") %>%
  group_by(sampleSize, ngroups, intercept,overdispersion, test) %>%
  summarise(mean.stat = mean(dispersion, na.rm=T))
d.pois$intercept <- fct_relevel(d.pois$intercept, "-1.5", "0", "1.5")
d.pois$ngroups <- fct_relevel(d.pois$ngroups, "10", "50", "100")
d.pois$sampleSize <- as.factor(as.numeric(d.pois$sampleSize))


# all groups
d.pois %>% filter(test != "Pear.stat.dispersion",
                  mean.stat <1000) %>%
  ggplot( aes(x=overdispersion, y=mean.stat, col=test, linetype=ngroups, shape=ngroups))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  scale_y_log10()+
  facet_grid(sampleSize~intercept, scales="free") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0.5, ymax = 1,
           alpha = .1,fill = "red")+
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Poisson: dispersion statistics", subtitle = "100 sim") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom") + 
  guides(color=guide_legend(nrow=2, byrow=TRUE))
#ggsave(here("figures", "5_glmmPois_dispersionStatsALL.jpeg"), width=12, height = 15)





