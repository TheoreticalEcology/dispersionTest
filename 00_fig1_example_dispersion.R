library(glmmTMB)
library(DHARMa)
library(MASS)
library(ggeffects)
library(tidyverse);
library(cowplot)
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

# FIGURE OVERDISPERSION v2 ----

ggplot(dataO, aes(x=x, y=y_over)) + geom_point() +
  geom_ribbon(data=pred_over, aes(x=x, y=predicted, ymin=conf.low, 
                                  ymax=conf.high), alpha=0.1, col="red", 
              linetype= "dotted", fill="red")+
  geom_ribbon(data=predOok, aes(x=x, y=predicted, ymin=conf.low2, 
                                ymax=conf.high2), alpha=0.1, col="blue", 
              linetype= "dotted",fill="blue") +
  ylab("Y") + xlab("X")  +
  
  annotate("rect", xmin=0.08, xmax=0.12, ymin=80, ymax=85, 
           col="red", alpha=0.4,fill="red")+
  annotate("text", x=0.05, y=90, label="Prediction Interval (95%)", hjust=0,
           size=5)+
  annotate("text", x=0.13, y=82.5, label="Overdispersed model", hjust=0)+
  
  annotate("rect", xmin=0.08, xmax=0.12, ymin=74, ymax=79, 
           col="blue", alpha=0.4,fill="blue")+
  annotate("text", x=0.13, y=76.5, label="Correct model", hjust=0)+
  theme_cowplot() 

ggsave("figures/00_fig1_overdisp_example2.png",device="png", height = 10, width = 12,units="cm")











#### OLD ---------------------------------------------------------------------

summary(model_over)

res_over <- simulateResiduals(model_over)
plot(res_over)
testDispersion(res_over, type = "PearsonChisq")
testDispersion(res_over, type = "DHARMa")






po<-ggplot(dataO, aes(y=y_over, x=x))+ #geom_point() +
  geom_line(data=predOok, aes(x=x, y=predicted), col="blue") +
  geom_ribbon(data=predOok, aes(x=x, y=predicted, 
                                ymin=conf.low, ymax=conf.high),
              col="blue", fill="blue",alpha=0.2) +
  geom_line(data=pred_over, aes(x=x, y=predicted), col="red") +
  geom_ribbon(data=pred_over, aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high),
              fill="red",alpha=0.4) +
  ylab("Y") + xlab("X") + ylim(0,10.2) +
  
  annotate("rect", xmin=0.08, xmax=0.12, ymin=9.8, ymax=10.2, 
           col="red", alpha=0.4,fill="red")+
  
  annotate("text", x=0.13, y=10, label="Overdispersed model", hjust=0)+
  
  annotate("rect", xmin=0.08, xmax=0.12, ymin=8.8, ymax=9.2, 
           col="blue", alpha=0.4,fill="blue")+
  
  annotate("text", x=0.13, y=9, label="Correct model", hjust=0)

po
ggsave("figures/00_fig1_overdisp_example.png", device="png", height = 7.5, width = 7.5,units="cm")



# FIGURE RESIDUALS ----


plot(res_over)
plot(res_ok)

pear_over <- data.frame(model = rep(c("overdisp", "correct"), each=100),
                         residual = c(residuals(model_over, type="pearson"),
                                      residuals(mover, type="pearson")),
                         predicted = c(predict(model_over),
                                       predict(mover)))
ggplot(pear_over, aes(x=predicted, y=residual, col=model)) +
  geom_point(alpha=0.7) +
  geom_quantile(quantiles = c(0.25,  0.75), size = 2) +
  geom_hline(yintercept = 0, linetype="dashed")+
  theme_cowplot()+
  theme(legend.position = "top")
ggsave("figures/Pearson_residuals_overdisp.jpeg", device="png", height = 10,
       width = 12,units="cm")


quant_over <- data.frame(model = rep(c("overdisp", "correct"), each=100),
                        residual = c(res_over$scaledResiduals,
                                     res_ok$scaledResiduals),
                        predicted = c(res_over$fittedPredictedResponse,
                                      res_ok$fittedPredictedResponse))
ggplot(quant_over, aes(x=predicted, y=residual, col=model)) +
  geom_point(alpha=0.7) +
  geom_quantile(quantiles = c(0.25,  0.75), size = 2) +
  geom_hline(yintercept = 0.5, linetype="dashed")+
  theme_cowplot()+
  theme(legend.position = "top")
ggsave("figures/Quantile_residuals_overdisp.jpeg", device="png", height = 10,
       width = 12,units="cm")

# resdata <- data.frame(res_over = ,
#                       pred_over = ,
#                       res_ok = ,
#                       pred_ok = )
# ggplot(resdata, aes(x=pred_over, y=res_over))+
#   geom_point() +
#   geom_hline(yintercept = 0.5, linetype="dotted") +
#   theme_cowplot() +
#   ggplot(resdata, aes(x=pred_ok, y=res_ok))+
#   geom_point(col="red") +
#   geom_hline(yintercept = 0.5, linetype="dotted") +
#   theme_cowplot()


## figures with insets for slopes ICs

slopes = data.frame(model = c("Underdispersed \n Poisson", "Corrected \n Conway-Maxwell-Poisson", 
                              "Overdispersed \n Poisson", "Corrected \n Negative Binomial"),
                    estimate = c(coef(model_under)[2],
                                 fixef(munder)$cond[2],
                                 coef(model_over)[2],
                                 coef(mover)[2]),
                    conf.low = c(confint(model_under)[2,1]-0.05,
                                 confint(munder)[2,1],
                                 confint(model_over)[2,1],
                                 confint(mover)[2,1]),
                    conf.up = c(confint(model_under)[2,2]+0.05,
                                 confint(munder)[2,2],
                                 confint(model_over)[2,2],
                                 confint(mover)[2,2]),
                    ypos = c(1.05,0.95,0.75,0.65),
                    col= c( "red", "blue","red", "blue")
                    ) %>%
  mutate(model = fct_relevel(model,  "Corrected \n Negative Binomial", "Overdispersed \n Poisson",
                             "Corrected \n Conway-Maxwell-Poisson", "Underdispersed \n Poisson"))

slopes %>% 
  ggplot(aes(x=estimate, y=model, col=col)) + geom_point() +
  scale_color_manual(values = c("red", "blue"))+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_vline(xintercept = 0.6, linetype="dotted")+
  geom_linerange(aes(xmin=conf.low, xmax=conf.up)) +
  xlab("Slope estimate")+ ylab("")+
  theme(legend.position = "none")


insU <- slopes %>% slice(1,2)  %>%
  ggplot(aes(x=estimate, y=model, col=col)) + geom_point() +
  scale_color_manual(values = c( "blue","red"))+
  scale_x_continuous(breaks=c(0,0.4,0.8))+
  scale_y_discrete( labels= c("Correc.", "Underdis."))+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_linerange(aes(xmin=conf.low, xmax=conf.up)) +
  xlab("")+ ylab("")+ ggtitle("Slope")+
  theme(legend.position = "none",
        axis.text = element_text(size=10),
        title = element_text(size=10, hjust=0.5),
        panel.background = element_rect(color="black", fill="white")) +
  annotate("segment", y=2, x=slopes[1,3], yend=1, xend=slopes[2,3], linetype="dotted")+
  annotate("segment", y=2, x=slopes[1,4], yend=1, xend=slopes[2,4], linetype="dotted")
#ggsave("coefs_under_example.jpeg", heigh=4, width=5)


pu + ylim(0,13) + inset_element(insU, left=0.1,bottom = 0.5,right = 0.8,top=1) 
#ggsave("figures/00_fig1_underdisp_example.png",device="png", height = 5, width = 5)

insO <- slopes %>% slice(3:4)  %>%
  ggplot(aes(x=estimate, y=model, col=col)) + geom_point() +
  scale_color_manual(values = c( "blue","red"))+
  scale_x_continuous(breaks=c(0,1,2), limits = c(0,2))+
  scale_y_discrete( labels= c("Correc.", "Overdis."))+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_linerange(aes(xmin=conf.low, xmax=conf.up)) +
  xlab("")+ ylab("")+ ggtitle("Slope")+
  theme(legend.position = "none",
        axis.text = element_text(size=10),
        title = element_text(size=10, hjust=0.5),
        panel.background = element_rect(color="black", fill="white")) +
  annotate("segment", y=2, x=slopes[3,3], yend=1, xend=slopes[4,3], linetype="dotted")+
  annotate("segment", y=2, x=slopes[3,4], yend=1, xend=slopes[4,4], linetype="dotted")
#ggsave("coefs_over_example.jpeg", heigh=4, width=5)

po + ylim(0,15) + inset_element(insO, left=0.1,bottom = 0.5,right = 0.8,top=1) 
#ggsave("figures/00_fig1_overdisp_example.png", device="png", height = 20, width = 20)



# Getting dispersion results ----

library(broom);library(broom.mixed)
tidy(model_under)
tidy(munder) %>% dplyr::select(-effect,-component)

tidy(model_over)
tidy(mover)

testDispersion(res_under, type = "PearsonChisq")
testDispersion(simulateResiduals(munder), type = "PearsonChisq")

testDispersion(res_over, type = "PearsonChisq")
testDispersion(simulateResiduals(mover), type = "PearsonChisq")
performance::check_overdispersion(mover)

########munder#########################
# Compare Dispersion Estimates
#########################
# Calculate the Pearson chi-square dispersion estimate
dispersion_under <- sum(residuals(model_under, type="pearson")^2) / model_under$df.residual
dispersion_over  <- sum(residuals(model_over, type="pearson")^2) / model_over$df.residual

cat("Estimated dispersion (underdispersed data):", dispersion_under, "\n")
cat("Estimated dispersion (overdispersed data):", dispersion_over, "\n")
