### Dispersion Tests Project
## Melina Leite
# Dez 24

library(DHARMa)
library(tidyverse)
library(here)
library(cowplot);
theme_set(theme_cowplot())
library(patchwork)

####################
##### Binomial #####
####################

load(here("data", "3_glmBin_power.Rdata"))


simuls.bin <- list()
for (i in 1:length(out.bin)) {
  sim <- out.bin[[i]]$simulations
  params <- strsplit(names(out.bin)[[i]], "_")[[1]]
  sim$intercept <- params[1]
  sim$sampleSize <- params[2]
  simuls.bin <- rbind(simuls.bin,sim)
}
names(simuls.bin)[names(simuls.bin)=="controlValues"] <- "overdispersion"


# figure Power 
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

# figure statistics

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
  ggtitle("Binomial", subtitle = "1000 sim; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom")
ggsave(here("figures", "3_glmBin_dispersionStats.jpeg"), width=10, height = 15)


# Dispersion statistics for 0 overdispersion

st.1.bin <- simuls.bin %>% 
  filter(overdispersion == 0) %>%
  dplyr::select(Pear.stat.dispersion,DHA.stat.dispersion,
                Ref.stat.dispersion, replicate,
                overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:3, names_to = "test", values_to = "Disp.stats") %>%
  group_by(sampleSize,intercept,overdispersion, test) %>%
  summarise(median.stat = median(Disp.stats, na.rm=T))
st.1.bin$intercept <- fct_relevel(st.1.bin$intercept, "-3", "-1.5", "0", "1.5", "3")
st.1.bin$sampleSize <- as.numeric(st.1.bin$sampleSize)

ggplot(st.1.bin, aes(x=sampleSize, y=median.stat, col=test))+
  geom_point()+ geom_line()+
  facet_grid(~intercept) +
  scale_x_log10()+
  geom_hline(yintercept = 1, linetype="dashed")


st.1.bin %>% filter(sampleSize %in% c(10,20,50,100)) %>%
  ggplot(aes(x=test, y=Disp.stats, col=test))+
  geom_violin()+
  facet_grid(sampleSize~intercept)+
  geom_hline(yintercept = 1, linetype="dashed") +
  scale_x_discrete(
    labels=c("Quantile Residuals", "Pearson Chi-squared",
             "Pearson Param. Bootstrap.")) +
  theme(axis.text.x = element_text(angle=45,hjust=1),
        legend.position = "none")




####################
##### Poisson  #####
####################


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


# figure Power 
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

# figure statistics

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
  facet_grid(sampleSize~intercept) +
  geom_hline(yintercept = 1, linetype="dotted", col="gray")+
  ggtitle("Poisson", subtitle = "1000 sim; Ntrials=10") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom")
ggsave(here("figures", "3_glmPois_dispersionStats.jpeg"), width=10, height = 15)


# Dispersion statistics for 0 overdispersion

st.1.pois <- simuls.pois %>% 
  filter(overdispersion == 0) %>%
  dplyr::select(Pear.stat.dispersion,DHA.stat.dispersion,
                Ref.stat.dispersion, replicate,
                overdispersion, intercept, sampleSize) %>%
  pivot_longer(1:3, names_to = "test", values_to = "Disp.stats")
st.1.pois$intercept <- fct_relevel(st.1.pois$intercept, "-3", "-1.5", "0", "1.5", "3")
st.1.pois$sampleSize <- as.factor(as.numeric(st.1.pois$sampleSize))


st.1.pois %>% filter(sampleSize %in% c(10,20,50,100)) %>%
  ggplot(aes(x=test, y=Disp.stats, col=test))+
  geom_violin()+
  facet_grid(sampleSize~intercept)+
  geom_hline(yintercept = 1, linetype="dashed") +
  scale_x_discrete(
    labels=c("Quantile Residuals", "Pearson Chi-squared",
             "Pearson Param. Bootstrap.")) +
  theme(axis.text.x = element_text(angle=45,hjust=1),
        legend.position = "none")
