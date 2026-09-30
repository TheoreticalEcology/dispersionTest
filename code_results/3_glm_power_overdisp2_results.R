### Dispersion Tests Project
## Results: power of dispersion tests under alternative OVERDISPERSION mechanisms
## negative binomial (Poisson) and beta-binomial (binomial), constant dispersion
## factor phi (Rev. 2, Comments 2 and 13). Simulations: 3_glm_power_overdisp2.R
## Melina Leite
# Sep 26

library(DHARMa)
library(tidyverse)
library(here)
library(cowplot);
theme_set(theme_cowplot())
library(patchwork)

load(here("data", "2_callibrated_alphaLevels.Rdata")) # callibrated alpha level
# created in script 2_glm_type1_results.R. Valid here because phi = 1 is exactly
# the Poisson / binomial null (NB and beta-binomial with phi = 1).

# plot Colors
source(here("functions_others", "plotColors.R"))
source(here("functions_others", "mcse.R"))

test.labels <- c("Sim-based residual variance",
                 "Chi-squared Pearson",
                 "Param. bootstrap Pearson")
names(col.intercept) <- c("-3", "-1.5", "0", "1.5", "3") # same colour per intercept in all figures

# phi is on a multiplicative scale -> log2 x-axis
phi.breaks <- c(1, 1.5, 2, 3, 5)
scale_x_phi <- function(...) scale_x_continuous(transform = "log2", breaks = phi.breaks, ...)
xlab.phi <- xlab(expression("Dispersion factor" ~ phi ~ "(log scale)"))


##############-###
#### helpers  ####
##############-###

# load one object from an .Rdata file without overwriting the global environment
loadObj <- function(file, obj){
  e <- new.env()
  load(file, envir = e)
  e[[obj]]
}

# intercept as factor with only the levels present, in numeric order
fctIntercept <- function(x){
  lev <- c("-3", "-1.5", "0", "1.5", "3")
  factor(as.character(x), levels = lev[lev %in% as.character(x)])
}

# DHARMaBenchmark list -> one long data.frame
# ctrl: name for the control variable (phi, overdispersion, underdispersion)
bindSims <- function(outList, ctrl = "phi"){
  simuls <- lapply(names(outList), function(nm){
    params <- strsplit(nm, "_")[[1]]
    outList[[nm]]$simulations %>%
      mutate(intercept = params[1], sampleSize = params[2])
  }) %>% bind_rows()
  names(simuls)[names(simuls) == "controlValues"] <- ctrl
  simuls
}

# proportion of significant tests + exact binomial 95% CI (MC error) + MCSE
# pcols: columns with p-values; alphaTab: NULL (alpha = 0.05) or calibrated alphas
powerTable <- function(simuls, pcols = c("Pear.p.val","DHA.p.val","Ref.p.val"),
                       alphaTab = NULL, ctrl = "phi"){
  p <- simuls %>% dplyr::select(all_of(pcols), replicate,
                                all_of(ctrl), intercept, sampleSize) %>%
    pivot_longer(all_of(pcols), names_to = "test", values_to = "p.val") %>%
    mutate(sampleSize = as.numeric(sampleSize))
  if (is.null(alphaTab)) {
    p$alpha <- 0.05
  } else {
    p <- p %>% left_join(alphaTab %>% ungroup() %>%
                           mutate(intercept = as.character(intercept),
                                  sampleSize = as.numeric(as.character(sampleSize))) %>%
                           dplyr::select(sampleSize, intercept, test, alpha),
                         by = c("sampleSize", "intercept", "test"))
    if (any(is.na(p$alpha))) stop("calibrated alpha missing for some settings")
  }
  p <- p %>% mutate(significance = p.val < alpha) %>%
    group_by(across(all_of(c("sampleSize", "intercept", ctrl, "test")))) %>%
    summarise(p.sig = sum(significance, na.rm = T),
              nsim = sum(!is.na(p.val)), .groups = "drop")
  p$prop.sig <- p$p.sig/p$nsim
  p$mcse <- sqrt(p$prop.sig * (1 - p$prop.sig) / p$nsim)
  p$p.bin0.05 <- p$conf.low <- p$conf.up <- NA_real_
  for (i in 1:nrow(p)) {
    btest <- binom.test(p$p.sig[i], n = p$nsim[i], p = 0.05)
    p$p.bin0.05[i] <- btest$p.value
    p$conf.low[i] <- btest$conf.int[1]
    p$conf.up[i] <- btest$conf.int[2]
  }
  p$intercept <- fctIntercept(p$intercept)
  p$sampleSize <- factor(p$sampleSize, levels = sort(unique(p$sampleSize)))
  p
}

# mean dispersion statistic + MCSE
statTable <- function(simuls, ctrl = "phi"){
  st <- simuls %>% dplyr::select(Pear.stat.dispersion, DHA.stat.dispersion,
                                 Ref.stat.dispersion, replicate,
                                 all_of(ctrl), intercept, sampleSize) %>%
    pivot_longer(1:3, names_to = "test", values_to = "Disp.stats")  %>%
    group_by(across(all_of(c("sampleSize", "intercept", ctrl, "test")))) %>%
    summarise(mean.stat = mean(Disp.stats, na.rm=T),
              mcse = mcse_mean(Disp.stats), .groups = "drop")
  st$intercept <- fctIntercept(st$intercept)
  st$sampleSize <- factor(as.numeric(st$sampleSize),
                          levels = sort(unique(as.numeric(st$sampleSize))))
  st
}

# full power figure (sampleSize x intercept)
powerFig <- function(dat, title, subtitle){
  ggplot(dat, aes(x=phi, y=prop.sig, col=test, linetype = model))+
    geom_point(alpha=0.7, shape=1) + geom_line(alpha=0.7) +
    geom_errorbar(aes(ymin=conf.low, ymax=conf.up),
                  alpha=0.7,
                  width = 0,
                  linetype = "solid",
                  show.legend = FALSE) +
    scale_color_manual( values= col.tests[c(4,1,2)], labels = test.labels)+
    scale_x_phi() +
    facet_grid(sampleSize~intercept) +
    geom_hline(yintercept = 0.5, linetype="dotted") +
    xlab.phi + ylab("Power") +
    ggtitle(title, subtitle = subtitle) +
    theme(panel.background = element_rect(color="black"),
          legend.position = "bottom")+
    guides(color=guide_legend(nrow=2, byrow=TRUE))
}

# dispersion statistics figure
statFig <- function(st, title, subtitle){
  ggplot(st, aes(x=phi, y=mean.stat, col=test))+
    geom_point(alpha=0.7, shape=1) + geom_line(alpha=0.7) +
    geom_errorbar(aes(ymin=mean.stat-1.96*mcse, ymax=mean.stat+1.96*mcse),
                  alpha=0.7,
                  width = 0,
                  linetype = "solid",
                  show.legend = FALSE) +
    # expected dispersion under the DGP: y = phi (exact for NB and beta-binomial)
    # (annotate, not geom_abline: abline's "intercept" column clashes with the facet)
    annotate("segment", x = 1, xend = max(phi.breaks), y = 1, yend = max(phi.breaks),
             linetype = "dashed", col = "gray50") +
    scale_color_manual( values= col.tests[c(4,1,2)], labels = test.labels)+
    scale_x_phi() +
    scale_y_continuous(transform = "log2") +
    facet_grid(sampleSize~intercept) +
    geom_hline(yintercept = 1, linetype="dotted", col="gray")+
    xlab.phi + ylab("Mean dispersion statistic (log scale)") +
    ggtitle(title, subtitle = subtitle) +
    theme(panel.background = element_rect(color="black"),
          legend.position = "bottom")
}


##############-###
#### Binomial ####
##############-###

out.bin <- loadObj(here("data", "3_glmBin_power_overdisp2.Rdata"), "out.bin")
simuls.bin <- bindSims(out.bin)

###### "RAW" Power (two-sided, alpha = 0.05) ######
p.bin <- powerTable(simuls.bin)

##### Callibrated Power (two-sided) #####
cp.bin <- powerTable(simuls.bin, alphaTab = alpha.bin)

##### One-sided power (alternative = "greater", alpha = 0.05) #####
gp.bin <- powerTable(simuls.bin, pcols = c("Pear.p.greater","DHA.p.greater","Ref.p.greater")) %>%
  mutate(test = str_replace(test, "p.greater", "p.val"))  # same labels/colours

##### figure ####
bind_rows(list(uncalibrated=p.bin, calibrated= cp.bin, `one-sided (greater)` = gp.bin),
          .id="model") %>%
  powerFig("Binomial power: overdispersion (beta-binomial)", "1000 sim; Ntrials = 10")
ggsave(here("figures", "3_glmBin_power_overdisp2.pdf"), width=10, height = 15)


##### figure statistics #####
st.bin <- statTable(simuls.bin)
statFig(st.bin, "Binomial: dispersion statistics, overdispersion (beta-binomial)",
        "1000 sim; Ntrials = 10; dashed = phi")
ggsave(here("figures", "3_glmBin_dispersionStats_overdisp2.pdf"), width=10, height = 15)


###############-##
#### Poisson  ####
##############-###

out.pois <- loadObj(here("data", "3_glmPois_power_overdisp2.Rdata"), "out.pois")
simuls.pois <- bindSims(out.pois)

##### "RAW" Power #####
p.pois <- powerTable(simuls.pois)

##### Calibrated power #####
cp.pois <- powerTable(simuls.pois, alphaTab = alpha.pois)

##### One-sided power (alternative = "greater") #####
gp.pois <- powerTable(simuls.pois, pcols = c("Pear.p.greater","DHA.p.greater","Ref.p.greater")) %>%
  mutate(test = str_replace(test, "p.greater", "p.val"))

##### figure ####
bind_rows(list(uncalibrated=p.pois, calibrated= cp.pois, `one-sided (greater)` = gp.pois),
          .id="model") %>%
  powerFig("Poisson power: overdispersion (negative binomial)", "1000 sim")
ggsave(here("figures", "3_glmPois_power_overdisp2.pdf"), width=10, height = 15)


##### figure statistics #####
st.pois <- statTable(simuls.pois)
statFig(st.pois, "Poisson: dispersion statistics, overdispersion (negative binomial)",
        "1000 simulations; dashed = phi")
ggsave(here("figures", "3_glmPois_dispersionStats_overdisp2.pdf"), width=10, height = 15)


## TABLE power + MCSE (Rev. 1, Comment 2: MC errors in the repository) ####
bind_rows(list(Poisson_uncalibrated = p.pois, Binomial_uncalibrated = p.bin,
               Poisson_calibrated = cp.pois, Binomial_calibrated = cp.bin,
               Poisson_greater = gp.pois, Binomial_greater = gp.bin), .id = "model") %>%
  separate(model, c("model", "calibration"), sep = "_") %>%
  mutate(dgp = ifelse(model == "Poisson", "negative binomial", "beta-binomial")) %>%
  dplyr::select(model, dgp, calibration, sampleSize, intercept, phi, test,
                nsim, p.sig, prop.sig, mcse, conf.low, conf.up) %>%
  write_csv(here("data", "3_glm_power_overdisp2_MCSE.csv"))



## FIGURE POWER together (as figures/3_glm_both_power.pdf) ####

pow <- bind_rows(list(Poisson_uncalibrated = p.pois,
                      Binomial_uncalibrated = p.bin,
                      Poisson_calibrated = cp.pois,
                      Binomial_calibrated = cp.bin), .id= "model") %>%
  separate(model, c("model", "calibration")) %>%
  mutate(model = fct_relevel(model, "Poisson", "Binomial"))

# intercept = 0, some sample sizes
sub.pow <- pow %>% filter(intercept == 0, sampleSize %in% c(10,100,1000))

ggplot(sub.pow, aes(x=phi, y=prop.sig, col=test, linetype = calibration)) +
  geom_point(alpha=0.7, shape=1) + geom_line(alpha=0.7) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up),
                alpha=0.7,
                width = 0,
                linetype = "solid",
                show.legend = FALSE) +
  ylab("Power") + xlab.phi +
  scale_x_phi() +
  scale_color_manual(values = col.tests[c(4,1,2)], labels = test.labels) +
  facet_grid(model~sampleSize,
             labeller = as_labeller(c("Binomial"= "Binomial (beta-binomial)",
                                      "Poisson" = "Poisson (neg. binomial)",
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

ggsave(here("figures", "3_glm_both_power_overdisp2.pdf"), width=10, height = 6)



## INDUCED DISPERSION (Rev. 2, Comment 13) ####
# Mean Pearson dispersion (glm, correct mean model) = population-level
# dispersion actually induced by each simulation setting.
# NB / beta-binomial: expected to be = phi for every intercept (by construction).
induced <- bind_rows(list(Poisson = st.pois, Binomial = st.bin), .id = "model") %>%
  filter(test == "Pear.stat.dispersion") %>%
  dplyr::select(model, sampleSize, intercept, phi, mean.stat, mcse)
write_csv(induced, here("data", "3_glm_induced_dispersion_overdisp2.csv"))

# same for the default DHARMa mechanism (Gaussian noise on the linear predictor,
# 3_glm_power.R), where the induced dispersion depends on the intercept
simuls.gauss <- list(
  Poisson  = bindSims(loadObj(here("data", "3_glmPois_power.Rdata"), "out.pois"),
                      ctrl = "overdispersion"),
  Binomial = bindSims(loadObj(here("data", "3_glmBin_power.Rdata"), "out.bin"),
                      ctrl = "overdispersion"))
st.gauss <- bind_rows(lapply(simuls.gauss, statTable, ctrl = "overdispersion"),
                      .id = "model")

induced.gauss <- st.gauss %>%
  filter(test == "Pear.stat.dispersion", sampleSize == 1000) %>%
  mutate(model = fct_relevel(model, "Poisson", "Binomial"))

g.ind1 <- ggplot(induced.gauss, aes(x = overdispersion, y = mean.stat, col = intercept)) +
  geom_point(shape = 1) + geom_line() +
  geom_hline(yintercept = 1, linetype = "dotted", col = "gray") +
  scale_color_manual(values = col.intercept, limits = names(col.intercept)) + # same legend in A and B
  scale_y_continuous(transform = "log2") +
  facet_grid(~ model) +
  ylab("Achieved Pearson dispersion") +
  xlab(expression("Gaussian noise on linear predictor," ~ sigma)) +
  ggtitle("A) Default DHARMa mechanism (3_glm_power.R)") +
  theme(panel.background = element_rect(color="black"))

g.ind2 <- induced %>% filter(sampleSize == 1000) %>%
  mutate(model = fct_relevel(model, "Poisson", "Binomial")) %>%
  ggplot(aes(x = phi, y = mean.stat, col = intercept)) +
  geom_point(shape = 1) + geom_line() +
  annotate("segment", x = 1, xend = 5, y = 1, yend = 5,
           linetype = "dashed", col = "gray50") +
  scale_color_manual(values = col.intercept, limits = names(col.intercept)) + # same legend in A and B
  scale_x_phi() + scale_y_continuous(transform = "log2") +
  facet_grid(~ model, labeller = as_labeller(c("Binomial"= "Binomial (beta-binomial)",
                                               "Poisson" = "Poisson (neg. binomial)"))) +
  ylab("Achieved Pearson dispersion") + xlab.phi +
  ggtitle("B) Negative binomial / beta-binomial (dashed = phi)") +
  theme(panel.background = element_rect(color="black"))

g.ind1 / g.ind2 + plot_layout(guides = "collect") +
  plot_annotation(caption = "Mean Pearson dispersion of the correct-mean GLM, n = 1000, 1000 simulations")
ggsave(here("figures", "3_glm_induced_dispersion_overdisp2.pdf"), width=8, height = 8)



## FIGURE POWER x DATA-GENERATING MECHANISM (Rev. 2, Comment 2) ####
# Power of the three tests against the INDUCED dispersion (mean Pearson
# dispersion at n = 1000 for the same intercept and DGP level), which puts all
# data-generating mechanisms on one common scale:
#   - Gaussian noise on the linear predictor (3_glm_power.R)
#   - negative binomial / beta-binomial (this script)
#   - Conway-Maxwell-Poisson / -binomial underdispersion (3_glm_power_underdisp.R)

# power tables with a common control column "level"
powGauss <- bind_rows(lapply(simuls.gauss, powerTable, ctrl = "overdispersion"),
                      .id = "model") %>%
  rename(level = overdispersion) %>% mutate(dgp = "Gaussian noise")

powPhi <- bind_rows(list(Poisson = p.pois, Binomial = p.bin), .id = "model") %>%
  rename(level = phi) %>% mutate(dgp = "NB / beta-binomial")

# induced dispersion for each (model, dgp, intercept, level)
indGauss <- induced.gauss %>% transmute(model = as.character(model), intercept,
                                        level = overdispersion, disp = mean.stat,
                                        dgp = "Gaussian noise")
indPhi <- induced %>% filter(sampleSize == 1000) %>%
  transmute(model, intercept, level = phi, disp = mean.stat, dgp = "NB / beta-binomial")

powList <- list(powGauss, powPhi)
indList <- list(indGauss, indPhi)

# underdispersion (if simulated): left side of the dispersion axis
fUnder <- c(Poisson = here("data", "3_glmPois_power_underdisp.Rdata"),
            Binomial = here("data", "3_glmBin_power_underdisp.Rdata"))
if (all(file.exists(fUnder))) {
  simU <- list(Poisson  = bindSims(loadObj(fUnder["Poisson"], "out.pois"), ctrl = "level"),
               Binomial = bindSims(loadObj(fUnder["Binomial"], "out.bin"), ctrl = "level"))
  powList[[3]] <- bind_rows(lapply(simU, powerTable, ctrl = "level"), .id = "model") %>%
    mutate(dgp = "CMP / CMB (underdispersion)")
  indList[[3]] <- bind_rows(lapply(simU, statTable, ctrl = "level"), .id = "model") %>%
    filter(test == "Pear.stat.dispersion", sampleSize == 1000) %>%
    transmute(model, intercept, level, disp = mean.stat, dgp = "CMP / CMB (underdispersion)")
}

powDGP <- bind_rows(powList) %>%
  left_join(bind_rows(indList), by = c("model", "dgp", "intercept", "level")) %>%
  mutate(model = fct_relevel(model, "Poisson", "Binomial"),
         dgp = fct_relevel(dgp, "Gaussian noise", "NB / beta-binomial"))

# induced dispersion of all DGPs (n = 1000) for the repository
bind_rows(indList) %>% write_csv(here("data", "3_glm_induced_dispersion_allDGPs.csv"))

powDGP %>%
  filter(intercept == 0, sampleSize %in% c(10, 100, 1000)) %>%
  ggplot(aes(x = disp, y = prop.sig, col = test, linetype = dgp, shape = dgp,
             group = interaction(test, dgp))) +
  geom_point(alpha = 0.7) + geom_line(alpha = 0.7) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.up), alpha = 0.7, width = 0,
                linetype = "solid", show.legend = FALSE) +
  geom_vline(xintercept = 1, col = "gray60") +
  geom_hline(yintercept = 0.05, linetype = "dashed", col = "gray40") +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  scale_color_manual(values = col.tests[c(4,1,2)], labels = test.labels) +
  scale_shape_manual(values = c(1, 2, 0)) +
  scale_x_continuous(transform = "log2", breaks = c(0.25, 0.5, 1, 2, 4, 8)) +
  facet_grid(model ~ sampleSize,
             labeller = as_labeller(c("Binomial"= "Binomial",
                                      "Poisson" = "Poisson",
                                      "10" = "n = 10",
                                      "100" = "n = 100",
                                      "1000" = "n = 1000"))) +
  xlab("Induced dispersion (mean Pearson dispersion at n = 1000; log scale)") +
  ylab("Power (two-sided, alpha = 0.05)") +
  theme(panel.background = element_rect(color="black"),
        legend.position = "bottom", legend.box = "vertical",
        legend.title = element_blank())
ggsave(here("figures", "3_glm_power_DGPcomparison.pdf"), width=10, height = 7)



## FIGURE DISPERSION STATISTICS together ####
# NB: for k/n binomial models, DHARMa 0.5.0 computes the sim-based statistic from
# observed/simulated successes minus the fitted PROPORTION (getObservedResponse vs
# getFitted), so between-observation variance of the mean enters the statistic and
# pulls it towards 1 when the slope is not 0 (Poisson not affected). This shows
# up as "Sim-based lower than Pearson" in the Binomial panels below.

dispersion <- bind_rows(list(Poisson = st.pois, Binomial = st.bin), .id= "model") %>%
  dplyr::select(-mcse) %>%
  pivot_wider(names_from = test, values_from = mean.stat) %>%
  mutate(reldif_DHA_Pear = (DHA.stat.dispersion - Pear.stat.dispersion)/DHA.stat.dispersion,
         reldif_DHA_Ref = (DHA.stat.dispersion - Ref.stat.dispersion)/DHA.stat.dispersion)

# text labels (Binomial panel only, as in 3_glm_power_results.R)
small.text <- data.frame(label = c("Sim-based higher than Pearson",
                                   "Sim-based lower than Pearson"),
                         model = "Binomial", x = 3, y = c(0.02, -0.02))

##### DHARMa stats X Pearson Chi-sq ####
#all results
ggplot(dispersion, aes(x=phi, y=reldif_DHA_Pear, col=sampleSize)) +
  geom_point(size=2, shape=1) + geom_line()+
  scale_x_phi() +
  facet_grid(intercept~model, scales="free")+
  ylab("Relative diff. Dispersion stats \n (Sim-based x Chi-squared Pearson)")+
  xlab.phi +
  geom_hline(yintercept = 0, linetype="dotted") +
  theme(panel.background = element_rect(color="black"))

###### figure to present intercept == 0 ####
dispersion %>% filter(intercept == 0) %>%
  ggplot(aes(x=phi, y=reldif_DHA_Pear, col=sampleSize)) +
  geom_point(size=2, shape=1) + geom_line()+
  scale_x_phi() +
  facet_grid(~model)+
  ylab("Relative diff. Dispersion stats \n (Sim-based x Chi-squared Pearson)")+
  xlab.phi +
  geom_hline(yintercept = 0, linetype="dotted")  +
  geom_text(data = small.text, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3.5, colour = "black") +
  theme(panel.background = element_rect(color="black"),
        text = element_text(size=12),
        axis.text = element_text(size=10))
ggsave(here("figures", "3_glm_DISP_diff_DHA-Pear_overdisp2.pdf"), height=4, width=9)


##### DHARMa stats X Pearson Bootstrapping ####
dispersion %>% filter(intercept == 0) %>%
  ggplot(aes(x=phi, y=reldif_DHA_Ref, col=sampleSize)) +
  geom_point(size=2, shape=1) + geom_line()+
  scale_x_phi() +
  facet_grid(~model)+
  ylab("Relative diff. Dispersion stats \n (Sim-based x Param. bootstrap Pearson)")+
  xlab.phi +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_text(data = small.text, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3.5, colour = "black") +
  theme(panel.background = element_rect(color="black"),
        text = element_text(size=12),
        axis.text = element_text(size=10))
ggsave(here("figures", "3_glm_DISP_diff_DHA-Ref_overdisp2.pdf"), height=4, width=9)



## SUMMARY for the text (intercept = 0, uncalibrated two-sided power) ####
p.sum <- bind_rows(list(Poisson = p.pois, Binomial = p.bin), .id = "model") %>%
  filter(intercept == 0, phi %in% c(1, 1.5, 2, 3)) %>%
  mutate(test = recode(test, DHA.p.val = "sim-based", Pear.p.val = "chisq Pearson",
                       Ref.p.val = "pboot Pearson"),
         power = sprintf("%.2f", prop.sig)) %>%
  dplyr::select(model, sampleSize, phi, test, power) %>%
  pivot_wider(names_from = test, values_from = power)
print(p.sum, n = Inf)
