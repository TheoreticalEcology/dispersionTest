# Dispersion testes paper ##
## Code for generating a figure as example of over and underdispersion for figure 1


library(glmmTMB)
library(DHARMa)
library(MASS)
library(ggeffects)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot() + 
            theme(panel.background = element_rect(color = "black")))
library(patchwork)


#### UNDERDISPERSION EXAMPLE ####

set.seed(123)
n <- 100
x <- runif(n, 0, 1)

# Define true model parameters
beta0 <- 1.7
beta1 <- 0.18
mu <- exp(beta0 + beta1 * x)

n_trials <- ceiling(max(mu))+1
y_under <- rbinom(n, size = n_trials, prob = mu / n_trials)

dataU <- data.frame(y_under = y_under, x = x)


# Fit Poisson GLM (ignoring that the true variance is lower than mu)
model_under <- glm(y_under ~ x, family = poisson, data=dataU)
summary(model_under)

plot(ggpredict(model_under, terms="x[all]", interval="prediction",
               ci_level = 0.95), show_data = T, ci_style="dot")

pred_under <- as.data.frame(ggpredict(model_under, terms="x[all]", 
                                      interval="prediction"))

# using Conway-Maxwell-Poisson
munder <- glmmTMB(y_under ~ x, data=dataU, family = compois())
summary(munder)

plot(ggpredict(munder, terms="x[all]", interval="prediction", 
               ci_level = 0.95), show_data = T, ci_style="dot") 

predUok <- as.data.frame(ggpredict(munder, terms="x[all]",
                                   interval="prediction"))
predUok$conf.low2 <- round(predUok$conf.low,0)
predUok$conf.high2 <- floor(predUok$conf.high)

# FIGURE UNDERDISP ----


ggplot(dataU, aes(x=x, y=y_under)) + geom_point() +
  stat_quantile(formula=y ~ x, quantiles=c(0.025, 0.975), linetype="dashed",
                col="black") +
  geom_line(data=pred_under, aes(x=x, y=predicted), size=1.2, col="darkred")+
  geom_ribbon(data=pred_under, aes(x=x, y=predicted, ymin=conf.low, 
                                   ymax=conf.high), alpha=0.1, col="red", 
              linetype= "dotted", fill="red")+
  ylab("Y") + xlab("X") + 
  ylim(0,18)+
  annotate("segment", x=0, xend=0.06, y=17, yend=17, linetype="dashed", 
           col="black") +
  annotate("text", x=0.08, y=17, label="data dispersion (central 95%)", hjust=0, size=4)+
  annotate("rect", xmin=0, xmax=0.06, ymin=13.5, ymax=15, 
           linetype="dotted",
           col="red", alpha=0.1, fill="red") +
  annotate("text", x=0.08, y=14.5, 
           label="assumed model dispersion \n (95% prediction interval)", 
           hjust=0, size=4)

ggsave("figures/00_fig1_underdisp_example.png",device="png", height = 10, width = 12,units="cm")




# OVERDISPERSION EXAMPLE -----

set.seed(123)
beta0 <- 1.5
beta1 <- 0
mu <- exp(beta0 + beta1 * x)

size_param <- 0.6  # small value to induce overdispersion
y_over <- rnbinom(n, mu = mu, size = size_param)


dataO <- data.frame(y_over = y_over, x= x)

# Fit Poisson GLM (ignoring the extra variance)
model_over <- glm(y_over ~ x, family = poisson, data=dataO)
summary(model_over)

plot(ggpredict(model_over, terms="x[all]", interval = "prediction"), show_data = T)

pred_over <- as.data.frame(ggpredict(model_over, terms="x[all]", interval = "prediction"))

# using negative binomial
mover <- glm.nb(y_over ~ x, data=dataO)
summary(mover)
plot(ggpredict(mover, terms="x[all]", interval = "prediction"), show_data = T)

predOok <- as.data.frame(ggpredict(mover, terms="x[all]",interval = "prediction"))
predOok$conf.low2 <- round(predOok$conf.low,0)
predOok$conf.high2 <- floor(predOok$conf.high)-9

# FIGURE OVERDISPERSION ----

ggplot(dataO, aes(x=x, y=y_over)) + geom_point() +
  stat_quantile(formula=y~exp(x), quantiles=c(0.025, 0.975), linetype="dashed",
                col="black") +
  geom_line(data=pred_over, aes(x=x, y=predicted), size=1.2, col="darkred")+
  geom_ribbon(data=pred_over, aes(x=x, y=predicted, ymin=conf.low, 
                                   ymax=conf.high), alpha=0.1, col="red", 
              linetype= "dotted", fill="red")+
  ylab("Y") + xlab("X")  + 
  ylim(0,50)+
  annotate("rect", xmin=0, xmax=0.06, ymin=40, ymax=44, 
           linetype="dotted", col="red", alpha=0.1,fill="red")+
  annotate("text", x=0.08, y=42, hjust=0, size=4,
           label="assumed model dispersion \n (95% prediction interval)") +
  annotate("text", x=0.08, y=48, label="data dispersion (central 95%)", 
           hjust=0, size=4) +
  annotate("segment", xend=0, x=0.06, y=48, yend=48, linetype="dashed",
          linewidth=0.9,  col="black")


ggsave("figures/00_fig1_overdisp_example.png",device="png", height = 10, width = 12,units="cm")

