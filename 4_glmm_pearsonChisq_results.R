### Dispersion Tests Project
## Melina Leite
# Dez 24

library(DHARMa)
library(tidyverse)
library(here)
library(cowplot);
theme_set(theme_cowplot())
library(patchwork)


# colors
cols = c("#556B2F", "#A2CD5A", "#1874CD", "#104E8B", "#7EC0EE",
                               "#CD5555", "#8B3A3A", "#FF8C69")


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
  btest <- binom.test(p.bin$p.sig[i], n=p.bin$nsim[i])
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
  btest <- binom.test(p.pois$p.sig[i], n=p.pois$nsim[i])
  p.pois$p.pois0.05[i] <- btest$p.value
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
            sd.stat = sd(Pear.stat.dispersion, na.rm=T))

st.pois <- simuls.pois %>% 
  mutate(ngroups = fct_relevel(ngroups, "10", "20", "50", "100")) %>%
  group_by(ngroups, overdispersion) %>%
  summarise(mean.stat = mean(Pear.stat.dispersion, na.rm=T),
            sd.stat = sd(Pear.stat.dispersion, na.rm=T))

disper <- bind_rows(list(Binomial = st.bin, Poisson=st.pois), .id="model")%>%
  mutate(model = fct_relevel(model, "Poisson", "Binomial"))



###### figure power #####


fig.power <- 
  power %>% filter(ngroups %in% c(10,50,100)) %>%
  ggplot( aes(x=overdispersion, y=prop.sig, col=model,
                               linetype=test))+
  geom_point(alpha=0.7) + geom_line(alpha=0.7) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up), width = 0.01) +
  annotate("rect", xmin = -0.05, xmax = 0.05, ymin = 0, ymax = 1,
           alpha = .1,fill = "blue")+
  scale_linetype_discrete(name= "Pearson test")+
  facet_wrap(~ngroups, labeller = as_labeller(c(`10`= "m = 10 groups",
                                                `50`= "m = 50 groups",
                                                `100`= "m = 100 groups",
                                                `Binomial` = "Binomial",
                                                `Poisson` = "Poisson"))) +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  theme(panel.background = element_rect(color="black"),
    legend.position ="right") +
  labs(tag="A)") +
  ylab("Power") + ylim(0,1)
fig.power




###### figure dispersion stat ####


fig.disp <- disper %>% filter(ngroups %in% c(10,50,100)) %>%
  ggplot(aes(x=overdispersion, y=mean.stat, col=model))+
  geom_point(alpha=0.7) + geom_line( alpha=0.7) +
  geom_errorbar(aes(ymin=mean.stat-sd.stat, ymax=mean.stat+sd.stat), width=0.02)+
  facet_grid(~ngroups, labeller = as_labeller(c(`10`= "m = 10 groups",
                                                `50`= "m = 50 groups",
                                                `100`= "m = 100 groups"))) +
  geom_hline(yintercept = 1, linetype="dotted", col="gray") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
           alpha = 0.1,fill = "red")+
  theme(panel.background = element_rect(color="black"),
        legend.position = "none")+
  labs(tag="B)") +
  ylab("Dispersion statistics")+
  scale_y_log10()
fig.disp
fig.power + fig.disp + plot_layout(ncol=1, guides="collect")  +
  plot_annotation(title="Pearson Chi-squared dispersion test for GLMMs",
                  theme = theme(plot.title = element_text(hjust=0.5)))


ggsave(here("figures", "4_glmm_pearsonChisq.jpeg"), width=12, height = 8)



## testing insets for type 1 error
# from: https://clarewest.github.io/blog-posts/ggplotInset.html
# library(purrr)
# 
# df <- power %>% filter(overdispersion == 0, ngroups %in% c(10,50,100)) 
# 
# ## A function to plot the inset
# get_inset2 <- function(df){
#   p <- ggplot(data = df %>% ungroup() %>%
#                 group_by(ngroups,model,test, overdispersion) %>%
#                 slice(1),
#               aes(x=test,y=prop.sig, col=model,group=overdispersion)) +
#     geom_point(position = position_dodge(width=0.5)) +
#     geom_errorbar(aes(ymin=conf.low,ymax=conf.up ), width=0.1)+
#      ylab("Type 1") + xlab("")+
#     theme(legend.position="none",
#           axis.text.x = element_text(angle=45, hjust=1),
#           panel.background = element_rect(color="black"))
#   return(p)
# }
# 
# 
# insets <- df %>%
#   split(f = .$ngroups) %>%
#   purrr::map(~annotation_custom2(
#     grob = ggplotGrob(get_inset2(.)),
#     data = expand.grid(ngroups=unique(.$ngroups),
#                       overdispersion=0,
#                       model=unique(.$model),
#                       test = unique(.$test),
#                       prop.sig = unique(.$prop.sig)),
#     ymin = 0, ymax=0.5, xmin=0.35, xmax=1))
# 
# fig.power + insets
