library(glmmTMB);
library(DHARMa)
library(MASS)
library(ggeffects)
library(tidyverse);
library(cowplot)
library(patchwork)

# Set seed for reproducibility
set.seed(123)
n <- 100
x <- runif(n, 0, 1)

# Define true model parameters
beta0 <- 1
beta1 <- 0.7
mu <- exp(beta0 + beta1 * x)

#-
#### UNDERDISPERSION EXAMPLE ####
#-
# For underdispersed data we can use a binomial distribution.
# A binomial count with n_trials and success probability p has mean = n_trials*p
# and variance = n_trials*p*(1-p). By setting p = mu/n_trials (with mu < n_trials),
# the variance becomes mu*(1 - mu/n_trials), which is smaller than mu.
n_trials <- 10  # Choose a fixed number of trials so that mu < n_trials for our range of x
y_under <- rbinom(n, size = n_trials, prob = mu / n_trials)

dataU <- data.frame(y_under = y_under, x = x)


# Fit Poisson GLM (ignoring that the true variance is lower than mu)
model_under <- glm(y_under ~ x, family = poisson, data=dataU)
summary(model_under)

res_under <- simulateResiduals(model_under)
plot(res_under)
testDispersion(res_under, type = "PearsonChisq")
testDispersion(res_under, type = "DHARMa")

pred_under <- as.data.frame(ggpredict(model_under, terms="x[all]"))
# altering the pred to make it look larger
pred_under$conf.high2 <- pred_under$conf.high+0.1
pred_under$conf.low2 <- pred_under$conf.low-0.1

# using Conway-Maxwell-Poisson
munder <- glmmTMB(y_under ~ x, data=dataU, family = compois())
summary(munder)
predUok <- as.data.frame(ggpredict(munder, terms="x[all]"))
r_ok <- simulateResiduals(munder)

# FIGURE UNDERDISP ----
theme_set(theme_cowplot() + 
            theme(panel.background = element_rect(color="black"),
                  axis.text.x = element_blank(),
                  axis.ticks = element_blank(),
                  axis.text.y = element_blank()))


pu <- ggplot(dataU, aes(y=y_under, x=x))+ #geom_point() +
  geom_line(data=pred_under, aes(x=x, y=predicted), col="red") +
  geom_ribbon(data=pred_under, aes(x=x, y=predicted, 
                                   ymin=conf.low2, ymax=conf.high2),
              fill="red",alpha=0.4) +
  geom_line(data=predUok, aes(x=x, y=predicted), col="blue") +
  geom_ribbon(data=predUok, aes(x=x, y=predicted,
                                ymin=conf.low, ymax=conf.high),
              col="blue", fill="blue",alpha=0.2) +
  ylab("Y") + xlab("X") + ylim(2,7.1) +
  
  annotate("rect", xmin=0.08, xmax=0.12, ymin=6.9, ymax=7.1, 
           col="red", alpha=0.4,fill="red")+
  
  annotate("text", x=0.13, y=7, label="Underdispersed model", hjust=0)+
  
  annotate("rect", xmin=0.08, xmax=0.12, ymin=6.4, ymax=6.6, 
           col="blue", alpha=0.4,fill="blue")+
  annotate("text", x=0.13, y=6.5, label="Correct model", hjust=0)
pu
ggsave("figures/00_fig1_underdisp_example.png",device="png", height = 7.5, width = 7.5,units="cm")


# FIGURE RESIDUALS ----

plot(res_under)
plot(r_ok)

pear_under <- data.frame(model = rep(c("underdisp", "correct"), each=100),
                         residual = c(residuals(model_under, type="pearson"),
                                      residuals(munder, type="pearson")),
                         predicted = c(predict(model_under),
                                        predict(munder)))
ggplot(pear_under, aes(x=predicted, y=residual, col=model)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype="dashed")+
  theme_cowplot() +
  theme(legend.position = "top")
ggsave("figures/Pearson_residuals_underdisp.jpeg", device="png", height = 10,
       width = 12,units="cm")

quant_under <- data.frame(model = rep(c("underdisp", "correct"), each=100),
                         residual = c(res_under$scaledResiduals,
                                       r_ok$scaledResiduals),
                         predicted = c(res_under$fittedPredictedResponse,
                                       r_ok$fittedPredictedResponse))

ggplot(quant_under, aes(x=predicted, y=residual, col=model)) +
  geom_point() +
  geom_hline(yintercept = 0.5, linetype="dashed")+
  theme_cowplot() +
  theme(legend.position = "top")
ggsave("figures/Quantile_residuals_underdisp.jpeg", device="png", height = 10,
       width = 12,units="cm")


#########################
# OVERDISPERSION EXAMPLE
#########################
# For overdispersed data we simulate counts from a negative binomial distribution.
# The negative binomial variance is mu + mu^2/size; with a small size parameter,
# the variance exceeds the mean.

size_param <- 1  # small value to induce overdispersion
y_over <- rnbinom(n, mu = mu, size = size_param)


dataO <- data.frame(y_over = y_over, x= x)

# Fit Poisson GLM (ignoring the extra variance)
model_over <- glm(y_over ~ x, family = poisson, data=dataO)
summary(model_over)

res_over <- simulateResiduals(model_over)
plot(res_over)
testDispersion(res_over, type = "PearsonChisq")
testDispersion(res_over, type = "DHARMa")

pred_over <- as.data.frame(ggpredict(model_over, terms="x[all]"))

# using negative binomial
mover <- glm.nb(y_over ~ x, data=dataO)
summary(mover)
predOok <- as.data.frame(ggpredict(mover, terms="x[all]"))
res_ok <- simulateResiduals(mover)

# FIGURE OVERDISPERSION ----

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
  geom_point() +
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
  geom_point() +
  geom_hline(yintercept = 0.5, linetype="dashed")+
  theme_cowplot()+
  theme(legend.position = "top")
ggsave("figures/Quantile_residuals_overdisp.jpeg", device="png", height = 10,
       width = 12,units="cm")

resdata <- data.frame(res_over = ,
                      pred_over = ,
                      res_ok = ,
                      pred_ok = )
ggplot(resdata, aes(x=pred_over, y=res_over))+
  geom_point() +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  theme_cowplot() +
  ggplot(resdata, aes(x=pred_ok, y=res_ok))+
  geom_point(col="red") +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  theme_cowplot()


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
