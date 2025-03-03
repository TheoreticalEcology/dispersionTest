library(glmmTMB);
library(DHARMa)
library(MASS)
library(ggeffects)
library(tidyverse);
library(cowplot)
theme_set(theme_cowplot() + theme(panel.background = element_rect(color="black")))
library(patchwork)

# Set seed for reproducibility
set.seed(123)
n <- 100
x <- runif(n, 0, 1)

# Define true model parameters
beta0 <- 1
beta1 <- 0.7
mu <- exp(beta0 + beta1 * x)

#########################
# UNDERDISPERSION EXAMPLE
#########################
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

# using Conway-Maxwell-Poisson
munder <- glmmTMB(y_under ~ x, data=dataU, family = compois())
summary(munder)
predUok <- as.data.frame(ggpredict(munder, terms="x[all]"))


pu <- ggplot(dataU, aes(y=y_under, x=x))+geom_point() +
  #geom_smooth(method="glm", method.args=list(family=poisson())) +
  geom_line(data=pred_under, aes(x=x, y=predicted), col="red") +
  geom_ribbon(data=pred_under, aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high),
              fill="red",alpha=0.4) +
  geom_line(data=predUok, aes(x=x, y=predicted), col="blue") +
  geom_ribbon(data=predUok, aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high),
              col="blue", fill="blue",alpha=0.2) +
  ylab("Y") + xlab("X") +ylim(0,10)
pu
#ggsave("underdisp_example.jpeg", height = 4, width = 5)




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


po<-ggplot(dataO, aes(y=y_over, x=x))+geom_point() +
  #geom_smooth(method="glm", method.args=list(family=poisson())) +
  geom_line(data=pred_over, aes(x=x, y=predicted), col="red") +
  geom_ribbon(data=pred_over, aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high),
              fill="red",alpha=0.4) +
  geom_line(data=predOok, aes(x=x, y=predicted), col="blue") +
  geom_ribbon(data=predOok, aes(x=x, y=predicted, ymin=conf.low, ymax=conf.high),
              col="blue", fill="blue",alpha=0.2) +
  ylab("Y") + xlab("X") +ylim(0,15)
po
#ggsave("overdisp_example.jpeg", height = 4, width = 5)


### FIgures ####

pu+po


slopes = data.frame(model = c("Underdispersed \n Poisson", "Corrected \n Conway-Maxwell-Poisson", 
                              "Overdispersed \n Poisson", "Corrected \n Negative Binomial"),
                    estimate = c(coef(model_under)[2],
                                 fixef(munder)$cond[2],
                                 coef(model_over)[2],
                                 coef(mover)[2]),
                    conf.low = c(confint(model_under)[2,1],
                                 confint(munder)[2,1],
                                 confint(model_over)[2,1],
                                 confint(mover)[2,1]),
                    conf.up = c(confint(model_under)[2,2],
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


pu + ylim(0,13) + inset_element(insU, left=-0.05,bottom = 0.5,right = 0.5,top=1) 
ggsave("figures/00_fig1_underdisp_example.jpeg", height = 5, width = 5)

insO <- slopes %>% slice(3:4)  %>%
  ggplot(aes(x=estimate, y=model, col=col)) + geom_point() +
  scale_color_manual(values = c( "blue","red"))+
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

po + ylim(0,15) + inset_element(insO, left=-0.05,bottom = 0.5,right = 0.5,top=1) 
ggsave("figures/00_fig1_overdisp_example.jpeg", height = 5, width = 5)



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
