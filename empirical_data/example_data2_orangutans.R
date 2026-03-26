library(glmmTMB)
library(DHARMa)
library(tidyverse)


#### model 2c ----

dat2c <- readxl::read_excel("empirical_data/data2_laumer_etal_2025_datasets.xlsx",
                           sheet = "model 2c")
dat2c$z.log.AgeAtThatTime2 <- dat2c$z.log.AgeAtThatTime^2

mod2c <- glmmTMB(NumberOfBodyParts ~ 
                   (z.log.AgeAtThatTime + I(z.log.AgeAtThatTime^2)) * Wild_Zoo + 
                   (1 + z.log.AgeAtThatTime + z.log.AgeAtThatTime2||FocalName) + 
                   (1|follow.in.subj), 
                 family=truncated_poisson,
                 data=dat2c)
summary(mod2c)

pear <- residuals(mod2c, type="pearson")

res2c <- simulateResiduals(mod2c)
plot(res2c)

(t2cSim <- testDispersion(res2c))
(t2cPear <- testDispersion(res2c, type = "PearsonChisq", alternative = "greater"))

res2cRefit <- simulateResiduals(mod2c, refit=T, n=100) 
(t2cNPear <- testDispersion(res2cRefit))



(tab2c <- data.frame(Model = rep("Zero-Truncated Poisson", 3),
                    Test = c("Sim-based residual variance", 
                             "Parametric Pearson residuals", 
                             "Nonparametric Pearson Residuals"),
                    Dispersion.parameter = round(c(t2cSim$statistic, 
                                                   t2cPear$statistic,
                                                   t2cNPear$statistic),3),
                    P.value = round(c(t2cSim$p.value, t2cPear$p.value, 
                                      t2cNPear$p.value),2)))



##### mod 2c with compois -----

mod2cCMP <- glmmTMB(NumberOfBodyParts ~ 
                   (z.log.AgeAtThatTime + I(z.log.AgeAtThatTime^2)) * Wild_Zoo + 
                   (1 + z.log.AgeAtThatTime + z.log.AgeAtThatTime2||FocalName) + 
                   (1|follow.in.subj), 
                 family=truncated_compois,
                 data=dat2c)
summary(mod2cCMP)

pear <- residuals(mod2cCMP, type="pearson")


res2cCMP <- simulateResiduals(mod2cCMP)
plot(res2cCMP)

(t2cCMPSim <- testDispersion(res2cCMP))
#pearson residuals is not implemented in truncated_compois models
#t2cCMPPear <- testDispersion(res2cCMP, type = "PearsonChisq", alternative = "greater")

#res2cCMPRefit <- simulateResiduals(mod2cCMP, refit=T, n=100) 
#t2cCMPNPear <- testDispersion(res2cCMPRefit)



tab2cCMP <- data.frame(Model = rep("Zero-Truncated CMP", 3),
                    Test = c("Sim-based residual variance", 
                             "Parametric Pearson residuals", 
                             "Nonparametric Pearson Residuals"),
                    Dispersion.parameter = round(c(t2cCMPSim$statistic, 
                                                   NA,
                                                   NA),3),
                    P.value = round(c(t2cCMPSim$p.value, NA, 
                                      NA),2))

##### mod 2c with NB -----

mod2cNB <- glmmTMB(NumberOfBodyParts ~ 
                      (z.log.AgeAtThatTime + I(z.log.AgeAtThatTime^2)) * Wild_Zoo + 
                      (1 + z.log.AgeAtThatTime + z.log.AgeAtThatTime2||FocalName) + 
                      (1|follow.in.subj), 
                    family=truncated_nbinom2(),
                    data=dat2c)
summary(mod2cNB)

res2cNB <- simulateResiduals(mod2cNB)
plot(res2cNB)

(t2cNBSim <- testDispersion(res2cNB))
(t2cNBPear <- testDispersion(res2cNB, type = "PearsonChisq", alternative = "greater"))

res2cNBRefit <- simulateResiduals(mod2cNB, refit=T, n=100) 
(t2cNBNPear <- testDispersion(res2cNBRefit))



(tab2cNB <- data.frame(Model = rep("Zero-Truncated NB", 3),
                       Test = c("Sim-based residual variance", 
                                "Parametric Pearson residuals", 
                                "Nonparametric Pearson Residuals"),
                       Dispersion.parameter = round(c(t2cNBSim$statistic, 
                                                      t2cNBPear$statistic,
                                                      t2cNBNPear$statistic),3),
                       P.value = round(c(t2cNBSim$p.value, t2cNBPear$p.value, 
                                         t2cNBNPear$p.value),2)))



# all results ----
bind_rows(tab1a, tab2a, tab2b, tab2c, tab2cCMP)


bind_rows(list(m2a= tab2a, m2b=tab2b,m2bCMP=tab2bCMP, m2c=tab2c, m2cCMP=tab2cCMP),.id="model")

bind_rows(list(m2c=tab2c, m2cCMP=tab2cCMP),.id="model")


AIC(mod2cCMP) - AIC(mod2c)

## mod 2c predictions ----


vals <-seq(min(dat2c$z.log.AgeAtThatTime), max(dat2c$z.log.AgeAtThatTime),length.out=30)
breaks <-seq(min(dat2c$z.log.AgeAtThatTime), max(dat2c$z.log.AgeAtThatTime),length.out=6)
values <- round(exp(mean(log(dat2c$AgeAtThatTime)) +  breaks*sd(log(dat2c$AgeAtThatTime))),2)


library(ggeffects)
pred2c <- ggemmeans(mod2c, terms=c("z.log.AgeAtThatTime [vals]", "Wild_Zoo")) |> 
  as.data.frame()
pred2cCMP <- ggemmeans(mod2cCMP, terms=c("z.log.AgeAtThatTime [vals]", "Wild_Zoo"))|> 
  as.data.frame()

preds <- bind_rows(list(trunc_Pois = pred2c, trunc_CMP = pred2cCMP), .id="model") %>%
  filter(group == "wild" | x<=3.4041 & x>=-2.1613)
  
  
ggplot(preds, aes(x=x, y=exp(predicted))) +
  geom_line(aes( col=group)) +
  geom_ribbon(aes(ymin=exp(conf.low), ymax=exp(conf.high),fill=group), alpha=0.2)+
  facet_grid(~model) +
  geom_point(data=dat2c, aes(x=z.log.AgeAtThatTime, y=NumberOfBodyParts,
                             col=Wild_Zoo), alpha=0.1) +
  scale_y_continuous(breaks = c(1,2,3,4,5,10,15)) +
  scale_x_continuous(breaks=breaks, labels=values)+
  scale_color_manual(values= c( "purple4","darkorange"))+
  scale_fill_manual(values= c( "purple4","darkorange")) +
  xlab("Age in years (log scale)") +
  ylab("Number of body parts")


broom.mixed::tidy(mod2c)
broom.mixed::tidy(mod2cCMP)










#### model 1a ----

dat1 <- readxl::read_excel("empirical_data/data2_laumer_etal_2025_datasets.xlsx",
                           sheet = "model 1a")
dat1$z.log.AgeAtThatTime_CS2 <- dat1$z.log.AgeAtThatTime_CS^2
dat1$time = log(dat1$VisibleObsTime/6)


##### NB2 (as in the paper)-----

mod1a <- glmmTMB(NrExplorationEvents_ALL_CS_NEWEST ~ 
                   (z.log.AgeAtThatTime_CS + z.log.AgeAtThatTime_CS2) * Wild_Zoo +
                   offset(time) + 
                   (z.log.AgeAtThatTime_CS + z.log.AgeAtThatTime_CS2||FocalName), 
                 ziformula= ~ z.log.AgeAtThatTime_CS + z.log.AgeAtThatTime_CS2+
                   Wild_Zoo, 
                 family=nbinom2,
                 data=dat1,
)
summary(mod1a)

res1a <- simulateResiduals(mod1a)
plot(res1a)
plotResiduals(res1a, form=~.)


(t1aSim <- testDispersion(res1a))
(t1aPear <- testDispersion(res1a, type = "PearsonChisq", alternative = "greater"))


# however nonparametric pearson residuals doesn't work
#res1aRefit <- simulateResiduals(mod1a, refit=T)
# why? offset!!
newData <- dat1
newData$NrExplorationEvents_ALL_CS_NEWEST = getFitted(mod1a)
refittedModel = update(mod1a, data = newData)


# result
(tab1a <- data.frame(Model = rep("Zero-inflated NB", 3),
                     Test = c("Sim-based residual variance", 
                              "Parametric Pearson residuals", 
                              "Nonparametric Pearson Residuals"),
                     Dispersion.parameter = round(c(t1aSim$statistic, 
                                                    t1aPear$statistic,
                                                    NA),2),
                     P.value = round(c(t1aSim$p.value, t1aPear$p.value, NA),3)))




##### NB1 (new)----

mod1a1 <- glmmTMB(NrExplorationEvents_ALL_CS_NEWEST ~ 
                    (z.log.AgeAtThatTime_CS + z.log.AgeAtThatTime_CS2) * Wild_Zoo +
                    offset(time) + 
                    (z.log.AgeAtThatTime_CS + z.log.AgeAtThatTime_CS2||FocalName), 
                  ziformula= ~ z.log.AgeAtThatTime_CS + z.log.AgeAtThatTime_CS2+
                    Wild_Zoo, 
                  family=nbinom1,
                  data=dat1,
)
summary(mod1a1)

res1a1 <- simulateResiduals(mod1a1)
plot(res1a1)
plotResiduals(res1a1, form=~.)


(t1a1Sim <- testDispersion(res1a1))
(t1a1Pear <- testDispersion(res1a1, type = "PearsonChisq", alternative = "greater"))


##### NB12 (new)----

mod1a12 <- glmmTMB(NrExplorationEvents_ALL_CS_NEWEST ~ 
                     (z.log.AgeAtThatTime_CS + z.log.AgeAtThatTime_CS2) * Wild_Zoo +
                     offset(time) + 
                     (z.log.AgeAtThatTime_CS + z.log.AgeAtThatTime_CS2||FocalName), 
                   ziformula= ~ z.log.AgeAtThatTime_CS + z.log.AgeAtThatTime_CS2+
                     Wild_Zoo, 
                   family=nbinom12,
                   data=dat1,
)
summary(mod1a12)

res1a12 <- simulateResiduals(mod1a12)
plot(res1a12)
plotResiduals(res1a12, form=~.)


(t1a12Sim <- testDispersion(res1a12))
#nbinom12 doesn't have Pearson residuals
#(t1a12Pear <- testDispersion(res1a12, type = "PearsonChisq", alternative = "greater"))



AIC(mod1a, mod1a1, mod1a12)





#### model 2a ----

dat2a <- readxl::read_excel("data2_laumer_etal_2025_datasets.xlsx",
                            sheet = "model 2a")
dat2a$z.log.AgeAtThatTime2 <- dat2a$z.log.AgeAtThatTime^2

mod2a <- glmmTMB(NumberOfManipulations ~ 
                   (z.log.AgeAtThatTime + z.log.AgeAtThatTime2) * Wild_Zoo +
                   (1+z.log.AgeAtThatTime + z.log.AgeAtThatTime2||FocalName)+
                   (1|follow.in.subj), 
                 family=truncated_nbinom2,
                 data=dat2a)
summary(mod2a)

res2a <- simulateResiduals(mod2a)
plot(res2a)
plotResiduals(res2a, form=~.)

(t2aSim <- testDispersion(res2a))
(t2aPear <- testDispersion(res2a, type = "PearsonChisq", alternative = "greater"))

res2aRefit <- simulateResiduals(mod2a, refit=T, n=100) 
(t2aNPear <- testDispersion(res2aRefit))



(tab2a <- data.frame(Model = rep("Zero-Truncated NB", 3),
                     Test = c("Sim-based residual variance", 
                              "Parametric Pearson residuals", 
                              "Nonparametric Pearson Residuals"),
                     Dispersion.parameter = round(c(t2aSim$statistic, 
                                                    t2aPear$statistic,
                                                    t2aNPear$statistic),3),
                     P.value = round(c(t2aSim$p.value, t2aPear$p.value, 
                                       t2aNPear$p.value),2)))




#### model 2b ----

dat2b <- readxl::read_excel("data2_laumer_etal_2025_datasets.xlsx",
                            sheet = "model 2b")
dat2b$z.log.AgeAtThatTime2 <- dat2b$z.log.AgeAtThatTime^2

mod2b <- glmmTMB(NumberOfManipulations ~ 
                   (z.log.AgeAtThatTime+ z.log.AgeAtThatTime2) * Wild_Zoo +
                   (1+z.log.AgeAtThatTime+ z.log.AgeAtThatTime2||FocalName)+
                   (1|follow.in.subj), 
                 family=truncated_poisson,
                 data=dat2b)
summary(mod2b)

res2b <- simulateResiduals(mod2b)
plot(res2b)

(t2bSim <- testDispersion(res2b))
(t2bPear <- testDispersion(res2b, type = "PearsonChisq", alternative = "greater"))

res2bRefit <- simulateResiduals(mod2b, refit=T, n=100) 
(t2bNPear <- testDispersion(res2bRefit))



(tab2b <- data.frame(Model = rep("Zero-Truncated Poisson", 3),
                     Test = c("Sim-based residual variance", 
                              "Parametric Pearson residuals", 
                              "Nonparametric Pearson Residuals"),
                     Dispersion.parameter = round(c(t2bSim$statistic, 
                                                    t2bPear$statistic,
                                                    t2bNPear$statistic),3),
                     P.value = round(c(t2bSim$p.value, t2bPear$p.value, 
                                       t2bNPear$p.value),2)))


##### model 2b - CMP

mod2bCMP <- glmmTMB(NumberOfManipulations ~ 
                      (z.log.AgeAtThatTime+ z.log.AgeAtThatTime2) * Wild_Zoo +
                      (1+z.log.AgeAtThatTime+ z.log.AgeAtThatTime2||FocalName)+
                      (1|follow.in.subj), 
                    family=truncated_compois,
                    data=dat2b)
summary(mod2bCMP)

res2bCMP <- simulateResiduals(mod2bCMP)
plot(res2bCMP)

(t2bCMPSim <- testDispersion(res2bCMP))


(tab2bCMP <- data.frame(Model = rep("Zero-Truncated CMP", 3),
                        Test = c("Sim-based residual variance", 
                                 "Parametric Pearson residuals", 
                                 "Nonparametric Pearson Residuals"),
                        Dispersion.parameter = round(c(t2bCMPSim$statistic, 
                                                       NA,
                                                       NA),3),
                        P.value = round(c(t2bCMPSim$p.value, NA, 
                                          NA),2)))



