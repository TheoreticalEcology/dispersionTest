### Dispersion Tests Project
## Results: power of dispersion tests under UNDERDISPERSION (GLMs)
## Simulations in code_simulations/3_glm_power_underdisp.R
##   Poisson : Conway-Maxwell-Poisson (mean parametrisation), nu = 1/(1-u)
##   Binomial: Conway-Maxwell-binomial (mean-matched), nu = 1/(1-u), Ntrials = 10
## Melina Leite
# Sep 26

library(DHARMa)
library(tidyverse)
library(here)
library(cowplot);
theme_set(theme_cowplot())
library(patchwork)

load(here("data", "2_callibrated_alphaLevels.Rdata")) # callibrated alpha level
# created in script 2_glm_type1_results.R. Valid here because u = 0 is exactly
# the Poisson / binomial null (CMP and CMB with nu = 1).

# plot Colors
source(here("functions_others", "plotColors.R"))
source(here("functions_others", "mcse.R"))

test.labels <- c("Sim-based residual variance",
                 "Chi-squared Pearson",
                 "Param. bootstrap Pearson")


##############-###
#### helpers  ####
##############-###

# DHARMaBenchmark list -> one long data.frame
bindSims <- function(outList){
  simuls <- list()
  for (i in 1:length(outList)) {
    sim <- outList[[i]]$simulations
    params <- strsplit(names(outList)[[i]], "_")[[1]]
    sim$intercept <- params[1]
    sim$sampleSize <- params[2]
    simuls <- rbind(simuls, sim)
  }
  names(simuls)[names(simuls)=="controlValues"] <- "underdispersion"
  simuls
}

# proportion of significant tests + exact binomial 95% CI (MC error)
# pcols: columns with p-values; alphaTab: NULL (alpha = 0.05) or calibrated alphas
powerTable <- function(simuls, pcols = c("Pear.p.val","DHA.p.val","Ref.p.val"),
                       alphaTab = NULL){
  p <- simuls %>% dplyr::select(all_of(pcols), replicate,
                                underdispersion, intercept, sampleSize) %>%
    pivot_longer(all_of(pcols), names_to = "test", values_to = "p.val") %>%
    mutate(sampleSize = as.numeric(sampleSize))
  if (is.null(alphaTab)) {
    p$alpha <- 0.05
  } else {
    p <- p %>% left_join(alphaTab %>% ungroup() %>%
                           mutate(intercept = as.character(intercept)),
                         by = c("sampleSize", "intercept", "test"))
  }
  p <- p %>% mutate(significance = p.val < alpha) %>%
    group_by(sampleSize, intercept, underdispersion, test) %>%
    summarise(p.sig = sum(significance, na.rm = T),
              nsim = length(p.val[!is.na(p.val)]), .groups = "drop")
  p$prop.sig <- p$p.sig/p$nsim
  p$p.bin0.05 <- p$conf.low <- p$conf.up <- NA_real_
  for (i in 1:nrow(p)) {
    btest <- binom.test(p$p.sig[i], n = p$nsim[i], p = 0.05)
    p$p.bin0.05[i] <- btest$p.value
    p$conf.low[i] <- btest$conf.int[1]
    p$conf.up[i] <- btest$conf.int[2]
  }
  p$intercept <- fct_relevel(p$intercept, "-3", "-1.5", "0", "1.5", "3")
  p$sampleSize <- as.factor(p$sampleSize)
  p
}

# mean dispersion statistic + MCSE
statTable <- function(simuls){
  st <- simuls %>% dplyr::select(Pear.stat.dispersion, DHA.stat.dispersion,
                                 Ref.stat.dispersion, replicate,
                                 underdispersion, intercept, sampleSize) %>%
    pivot_longer(1:3, names_to = "test", values_to = "Disp.stats")  %>%
    group_by(sampleSize, intercept, underdispersion, test) %>%
    summarise(mean.stat = mean(Disp.stats, na.rm=T),
              mcse = mcse_mean(Disp.stats), .groups = "drop")
  st$intercept <- fct_relevel(st$intercept, "-3", "-1.5", "0", "1.5", "3")
  st$sampleSize <- as.factor(as.numeric(st$sampleSize))
  st
}

# full power figure (sampleSize x intercept)
powerFig <- function(dat, title, subtitle){
  ggplot(dat, aes(x=underdispersion, y=prop.sig, col=test, linetype = model))+
    geom_point(alpha=0.7, shape=1) + geom_line(alpha=0.7) +
    geom_errorbar(aes(ymin=conf.low, ymax=conf.up),
                  alpha=0.7,
                  width = 0,
                  linetype = "solid",
                  show.legend = FALSE) +
    scale_color_manual( values= col.tests[c(4,1,2)], labels = test.labels)+
    facet_grid(sampleSize~intercept) +
    geom_hline(yintercept = 0.5, linetype="dotted") +
    xlab("Underdispersion (u; nu = 1/(1-u))") + ylab("Power") +
    ggtitle(title, subtitle = subtitle) +
    theme(panel.background = element_rect(color="black"),
          legend.position = "bottom")+
    guides(color=guide_legend(nrow=2, byrow=TRUE))
}

# dispersion statistics figure
statFig <- function(st, title, subtitle){
  ggplot(st, aes(x=underdispersion, y=mean.stat, col=test))+
    geom_point(alpha=0.7, shape=1) + geom_line(alpha=0.7) +
    geom_errorbar(aes(ymin=mean.stat-1.96*mcse, ymax=mean.stat+1.96*mcse),
                  alpha=0.7,
                  width = 0,
                  linetype = "solid",
                  show.legend = FALSE) +
    # expected dispersion ratio under the DGP (1-u), ignoring small-mean limits
    # (annotate, not geom_abline: abline's "intercept" column clashes with the facet)
    annotate("segment", x = 0, xend = 0.9, y = 1, yend = 0.1,
             linetype = "dashed", col = "gray50") +
    scale_color_manual( values= col.tests[c(4,1,2)], labels = test.labels)+
    facet_grid(sampleSize~intercept) +
    geom_hline(yintercept = 1, linetype="dotted", col="gray")+
    xlab("Underdispersion (u)") + ylab("Mean dispersion statistic") +
    ggtitle(title, subtitle = subtitle) +
    theme(panel.background = element_rect(color="black"),
          legend.position = "bottom")
}


##############-###
#### Binomial ####
##############-###

load(here("data", "3_glmBin_power_underdisp.Rdata")) # simulated data
simuls.bin <- bindSims(out.bin)

###### "RAW" Power (two-sided, alpha = 0.05) ######
p.bin <- powerTable(simuls.bin)

##### Callibrated Power (two-sided) #####
cp.bin <- powerTable(simuls.bin, alphaTab = alpha.bin)

##### One-sided power (alternative = "less", alpha = 0.05) #####
lp.bin <- powerTable(simuls.bin, pcols = c("Pear.p.less","DHA.p.less","Ref.p.less")) %>%
  mutate(test = str_replace(test, "p.less", "p.val"))  # same labels/colours

##### figure ####
bind_rows(list(uncalibrated=p.bin, calibrated= cp.bin, `one-sided (less)` = lp.bin),
          .id="model") %>%
  powerFig("Binomial power: underdispersion (CMB)", "1000 sim; Ntrials = 10")
ggsave(here("figures", "3_glmBin_power_underdisp.pdf"), width=10, height = 15)


##### figure statistics #####
st.bin <- statTable(simuls.bin)
statFig(st.bin, "Binomial: dispersion statistics, underdispersion (CMB)",
        "1000 sim; Ntrials=10; dashed = 1-u")
ggsave(here("figures", "3_glmBin_dispersionStats_underdisp.pdf"), width=10, height = 15)


###############-##
#### Poisson  ####
##############-###

load(here("data", "3_glmPois_power_underdisp.Rdata"))
simuls.pois <- bindSims(out.pois)

##### "RAW" Power #####
p.pois <- powerTable(simuls.pois)

##### Calibrated power #####
cp.pois <- powerTable(simuls.pois, alphaTab = alpha.pois)

##### One-sided power (alternative = "less") #####
lp.pois <- powerTable(simuls.pois, pcols = c("Pear.p.less","DHA.p.less","Ref.p.less")) %>%
  mutate(test = str_replace(test, "p.less", "p.val"))

##### figure ####
bind_rows(list(uncalibrated=p.pois, calibrated= cp.pois, `one-sided (less)` = lp.pois),
          .id="model") %>%
  powerFig("Poisson power: underdispersion (CMP)", "1000 sim")
ggsave(here("figures", "3_glmPois_power_underdisp.pdf"), width=10, height = 15)


##### figure statistics #####
st.pois <- statTable(simuls.pois)
statFig(st.pois, "Poisson: dispersion statistics, underdispersion (CMP)",
        "1000 simulations; dashed = 1-u")
ggsave(here("figures", "3_glmPois_dispersionStats_underdisp.pdf"), width=10, height = 15)



## INDUCED DISPERSION (Rev. 2, Comment 13) ####
# Mean Pearson dispersion (glm, correct model) = population-level dispersion
# actually induced by each simulation setting
induced <- bind_rows(list(Poisson = st.pois, Binomial = st.bin), .id = "model") %>%
  filter(test == "Pear.stat.dispersion") %>%
  dplyr::select(model, sampleSize, intercept, underdispersion, mean.stat, mcse)
write_csv(induced, here("data", "3_glm_induced_dispersion_underdisp.csv"))

# large n only (n = 1000, 10000): the achieved dispersion by intercept
induced %>% filter(sampleSize %in% c(1000, 10000)) %>%
  ggplot(aes(x = underdispersion, y = mean.stat, col = intercept)) +
  geom_point(shape = 1) + geom_line() +
  annotate("segment", x = 0, xend = 0.9, y = 1, yend = 0.1,
           linetype = "dashed", col = "gray50") +
  scale_color_manual(values = col.intercept) +
  facet_grid(sampleSize ~ model,
             labeller = labeller(sampleSize = function(x) paste("n =", x))) +
  ylab("Achieved Pearson dispersion") + xlab("Underdispersion (u)") +
  theme(panel.background = element_rect(color="black"))
ggsave(here("figures", "3_glm_induced_dispersion_underdisp.pdf"), width=8, height = 6)



## FIGURE POWER together ####

pow <- bind_rows(list(Poisson_uncalibrated = p.pois,
                      Binomial_uncalibrated = p.bin,
                      Poisson_calibrated = cp.pois,
                      Binomial_calibrated = cp.bin), .id= "model") %>%
  separate(model, c("model", "calibration")) %>%
  mutate(model = fct_relevel(model, "Poisson", "Binomial"))

# intercept = 0, some sample sizes
sub.pow <- pow %>% filter(intercept == 0, sampleSize %in% c(10,100,1000))

ggplot(sub.pow, aes(x=underdispersion, y=prop.sig, col=test, linetype = calibration)) +
  geom_point(alpha=0.7, shape=1) + geom_line(alpha=0.7) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up),
                alpha=0.7,
                width = 0,
                linetype = "solid",
                show.legend = FALSE) +
  ylab("Power") + xlab("Underdispersion") +
  scale_color_manual(values = col.tests[c(4,1,2)], labels = test.labels) +
  facet_grid(model~sampleSize,
             labeller = as_labeller(c("Binomial"= "Binomial (CMB)",
                                      "Poisson" = "Poisson (CMP)",
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

ggsave(here("figures", "3_glm_both_power_underdisp.pdf"), width=10, height = 6)



## FIGURE POWER over- AND underdispersion together ####
# x-axis: underdispersion plotted as -u (left), overdispersion as sd of the
# added noise (right). NB: the two halves are on different parameter scales.

load(here("data", "3_glmBin_power.Rdata"))   # overwrites out.bin
load(here("data", "3_glmPois_power.Rdata"))  # overwrites out.pois
# (helpers name the control column "underdispersion"; here it holds the
#  overdispersion sd, and is renamed to the signed axis "disp" below)
toSigned <- function(p, sign) p %>% rename(disp = underdispersion) %>%
  mutate(disp = sign * disp)

pow.over <- bind_rows(list(Poisson  = powerTable(bindSims(out.pois)),
                           Binomial = powerTable(bindSims(out.bin))),
                      .id = "model") %>% toSigned(1) %>% mutate(direction = "over")
pow.under <- bind_rows(list(Poisson = p.pois, Binomial = p.bin), .id = "model") %>%
  toSigned(-1) %>% mutate(direction = "under")

bind_rows(pow.over, pow.under) %>%
  filter(intercept == 0, sampleSize %in% c(10,100,1000)) %>%
  mutate(model = fct_relevel(model, "Poisson", "Binomial")) %>%
  ggplot(aes(x = disp, y = prop.sig, col = test,
             group = interaction(test, direction))) +
  geom_point(alpha=0.7, shape=1) + geom_line(alpha=0.7) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up), alpha=0.7, width = 0,
                show.legend = FALSE) +
  geom_vline(xintercept = 0, col = "gray60") +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  scale_color_manual(values = col.tests[c(4,1,2)], labels = test.labels) +
  facet_grid(model ~ sampleSize,
             labeller = as_labeller(c("Binomial"= "Binomial",
                                      "Poisson" = "Poisson",
                                      "10" = "n = 10",
                                      "100" = "n = 100",
                                      "1000" = "n = 1000"))) +
  xlab("← underdispersion (u)          overdispersion (σ) →") +
  ylab("Power (two-sided, α = 0.05)") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom", legend.title = element_blank())
ggsave(here("figures", "3_glm_both_power_over_under.pdf"), width=10, height = 6)
