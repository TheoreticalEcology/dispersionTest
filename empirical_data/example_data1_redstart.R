library(blmeco)
library(DHARMa)
library(MASS)
library(glmmTMB)
library(ggplot2)


### 
data("redstart")
redstart$elevation.z <- as.numeric(scale(redstart$elevation))
redstart$elevation.z2 <- redstart$elevation.z^2


#### Poisson model ----

mod5 <- glm(counts ~ elevation.z + elevation.z2 + forests,
            data= redstart,
            family=poisson())
summary(mod5)
anova(mod5)

# model presents significant effexts of elevation and forest (p=0.02)

res5 <- simulateResiduals(mod5)

testDispersion(res5) # simulation-based residual variance test
testDispersion(res5, type = "PearsonChisq")
res5Refit <- simulateResiduals(mod5, refit = T) 
testDispersion(res5Refit) # non-parametric Pearson

# all tests presented very similar results of overdispersion

plot(res5)
testSpatialAutocorrelation(res5, x=~x, y=~y)


#### Spatial model ----

redstart$x.s <- scale(redstart$x)
redstart$y.s <- scale(redstart$y)

redstart$pos <- numFactor(redstart$x.s, redstart$y.s)
redstart$group <- factor(rep(1, nrow(redstart)))

modS <- glmmTMB(counts ~ elevation.z + elevation.z2 + forests + exp(pos + 0 | group),
            data= redstart,
            family=poisson())

summary(modS)

resS <- simulateResiduals(modS)

testDispersion(resS) # simulation-based residual variance test
testDispersion(resS, type = "PearsonChisq", alternative = "greater")
resSRefit <- simulateResiduals(modS, refit = T, n=100) 
testDispersion(resSRefit) # non-parametric Pearson


plot(resS)
resSs <- simulateResiduals(modS, rotation = "estimated")
testSpatialAutocorrelation(resSs, x=~x.s, y=~y.s)







#### Negative binomial ----

mod6 <- glm.nb(counts ~ elevation.z + elevation.z2 + forests,
            data= redstart)
summary(mod6)
anova(mod6)
# no more significant effect of forest (p=0.10), only elevation.

res6 <- simulateResiduals(mod6)

testDispersion(res6) 
testDispersion(res6, type = "PearsonChisq")
res6Refit <- simulateResiduals(mod6, refit = T) 
testDispersion(res6Refit) # non-parametric Pearson

# all tests presented very similar results, no overdispersion (1.02 to 1.04)

plot(res6)
res6s <- simulateResiduals(mod6, rotation = "estimated")
testSpatialAutocorrelation(res6, x=~x, y=~y)



library(ggeffects)
pred5 <- ggpredict(mod5, terms="forests", back_transform = T) |> as.data.frame()
pred6 <- ggpredict(mod6, terms="forests", back_transform = T) |> as.data.frame()

ggplot(redstart, aes(x=forests, y=counts)) + geom_point(alpha=0.2) +
  geom_line(data=pred5, aes(x=x, y=predicted), col="darkorange") +
  geom_ribbon(data=pred5, aes(x=x,y=predicted, ymin=conf.low, 
                              ymax=conf.high), alpha=0.2, fill="darkorange") +
  
  geom_line(data=pred6, aes(x=x, y=predicted), col="purple4") +
  geom_ribbon(data=pred6, aes(x=x,y=predicted, ymin=conf.low, 
                             ymax=conf.high), alpha=0.2, fill="purple4") +
  ylim(0,3) +
  ylab("Number of breeding pairs") + xlab("Forest cover")

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
coefs <- bind_rows(list(Poisson=broom::tidy(mod5),
               `negative binomial` = broom::tidy(mod6)), .id="model") %>%
  filter(term !="(Intercept)")
ggplot(coefs, aes(x=estimate, y=term, col=model, group=model)) + 
  geom_point(position = position_dodge(0.7)) +
  geom_linerange(aes(xmin=estimate - 1.96*std.error,xmax=estimate+1.96*std.error),position = position_dodge(0.7)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme(panel.background = element_rect(color="black"))
 





pred5 <- ggpredict(mod5, terms="forests", back_transform = T) |> as.data.frame()
pred6 <- ggpredict(mod6, terms="forests", back_transform = T) |> as.data.frame()
predS <- ggpredict(modS, terms="forests", back_transform = T) |> as.data.frame()

ggplot(redstart, aes(x=forests, y=counts)) + geom_point(alpha=0.2) +
  geom_line(data=pred5, aes(x=x, y=predicted), col="darkorange") +
  geom_ribbon(data=pred5, aes(x=x,y=predicted, ymin=conf.low, 
                              ymax=conf.high), alpha=0.2, fill="darkorange") +
  
  geom_line(data=pred6, aes(x=x, y=predicted), col="purple4") +
  geom_ribbon(data=pred6, aes(x=x,y=predicted, ymin=conf.low, 
                              ymax=conf.high), alpha=0.2, fill="purple4") +
  geom_line(data=predS, aes(x=x, y=predicted), col="darkgreen") +
  geom_ribbon(data=predS, aes(x=x,y=predicted, ymin=conf.low, 
                              ymax=conf.high), alpha=0.2, fill="darkgreen") +
  ylim(0,3) +
  ylab("Number of breeding pairs") + xlab("Forest cover")

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
coefs <- bind_rows(list(Poisson=broom::tidy(mod5),
                        `negative binomial` = broom::tidy(mod6), spatialglm=broom.mixed::tidy(modS)), .id="model") %>%
  filter(term !="(Intercept)")
ggplot(coefs, aes(x=estimate, y=term, col=model, group=model)) + 
  geom_point(position = position_dodge(0.7)) +
  geom_linerange(aes(xmin=estimate - 1.96*std.error,xmax=estimate+1.96*std.error),position = position_dodge(0.7)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme(panel.background = element_rect(color="black"))


















