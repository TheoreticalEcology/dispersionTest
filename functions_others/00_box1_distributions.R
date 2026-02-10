
# Figure with distributions for the example Box1
library(tidyverse)
library(COMPoissonReg)
library(mpcmp)

set.seed(2)

distribs <- data.frame(distrib = factor(rep(c("Conway-Maxwell-Poisson",
                                              "Poisson", "Negative Binomial"), 
                                            each=1000), 
                                        levels = c("Conway-Maxwell-Poisson",
                                        "Poisson", "Negative Binomial")), 
                       values = c(rcomp(1000,1,2), rpois(1000,1),rnbinom(1000,mu=1,size=1.5) )) 


notes  <- data.frame(distrib = factor(c("Conway-Maxwell-Poisson",
                                 "Poisson", "Negative Binomial"),
                                 levels = c("Conway-Maxwell-Poisson",
                                            "Poisson", "Negative Binomial")),           
  label = c("mean = 1 \n disp = 0.5", "mean = 1 \n disp = 1", "mean = 1 \n disp = 1.5"),
  x = c(4, 4, 4),                      # Coordenada X no gráfico
  y = c(400, 400, 400)                    # Coordenada Y no gráfico
)

ggplot(distribs, aes(x=values, fill=distrib)) + geom_bar() +
  facet_grid(~distrib) +
  scale_y_continuous(name="")+
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10), 
                     name="Number of individuals per sample")  +
  scale_fill_manual(values=c("#0C639E","#26730D", "#C26A42")) +
  theme(legend.position = "none",
       # axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  geom_text(data = notes, 
            aes(x = x, y = y, label = label))
  
ggsave("figures/box1_distributions.png", heigh=9, width = 28, units="cm")



library(glmmTMB)
library(DHARMa)

m1 <- glmmTMB(values ~ 1, family=compois(), data=distribs[distribs$distrib=="Conway-Maxwell-Poisson",])
summary(m1)
exp(fixef(m1)[[1]])
res1 <- simulateResiduals(m1)
testDispersion(res1, type="PearsonChisq")
testDispersion(res1)


m2 <- glmmTMB(values ~ 1, family=poisson(), data=distribs[distribs$distrib=="Poisson",])
exp(fixef(m2)[[1]])
res2 <- simulateResiduals(m2)
testDispersion(res2, type="PearsonChisq")
testDispersion(res2)

m3 <- glmmTMB(values ~ 1, family=nbinom2(), data=distribs[distribs$distrib=="Negative Binomial",])
summary(m3)
exp(fixef(m3)[[1]])
res3 <- simulateResiduals(m3)
testDispersion(res3, type="PearsonChisq")
testDispersion(res3)
