#####      Study case 2: wild and zoo-housed orangutan behavior      #####
# Melina Leite
# March 2026

# The data was available in the supplementary material of the original study:
# Laumer, I. B., Kansal, S., Van Cauwenberghe, A., Rahmaeti, T., Setia, T. M., Mundry, R., Haun, D., & Schuppli, C. (2025). Wild and zoo-housed orangutans differ in how they explore objects. Scientific Reports, 15(1), 14853. https://doi.org/10.1038/s41598-025-97926-z

# Models' structure followed the supplementary material of the original study (Table S4a, model 2c). 


# Loading Packages ----
library(glmmTMB) # modeling
library(DHARMa) # dispersion tests
library(tidyverse) # coding 
library(readxl) # read excel file


#### Dataset ----

dat2c <- read_excel("case_studies/data2_laumer_etal_2025_datasets.xlsx",
                            sheet = "model 2c")
dat2c$z.log.AgeAtThatTime2 <- dat2c$z.log.AgeAtThatTime^2



#### Original model: truncated Poisson ----

mod2c <- glmmTMB(NumberOfBodyParts ~ 
                   (z.log.AgeAtThatTime + I(z.log.AgeAtThatTime^2)) * Wild_Zoo + 
                   (1 + z.log.AgeAtThatTime + z.log.AgeAtThatTime2||FocalName) + 
                   (1|follow.in.subj), 
                 family=truncated_poisson,
                 data=dat2c)
summary(mod2c)


###### dispersion tests -----

res2c <- simulateResiduals(mod2c)
res2cRefit <- simulateResiduals(mod2c, refit=T, n=100) 


(t2cSim <- testDispersion(res2c)) # simulation-based residual variance
(t2cPear <- testDispersion(res2c, type = "PearsonChisq", 
                           alternative = "greater")) # parametric Pearson
(t2cNPear <- testDispersion(res2cRefit)) # non-parametric Pearson



(tab2c <- data.frame(Model = rep("Zero-Truncated Poisson", 3),
                    Test = c("Sim-based residual variance", 
                             "Parametric Pearson residuals", 
                             "Nonparametric Pearson Residuals"),
                    Dispersion.parameter = round(c(t2cSim$statistic, 
                                                   t2cPear$statistic,
                                                   t2cNPear$statistic),3),
                    P.value = round(c(t2cSim$p.value, t2cPear$p.value, 
                                      t2cNPear$p.value),2)))



##### Alternative: Conway-Maxwel-Poisson -----

mod2cCMP <- glmmTMB(NumberOfBodyParts ~ 
                   (z.log.AgeAtThatTime + I(z.log.AgeAtThatTime^2)) * Wild_Zoo + 
                   (1 + z.log.AgeAtThatTime + z.log.AgeAtThatTime2||FocalName) + 
                   (1|follow.in.subj), 
                 family=truncated_compois,
                 data=dat2c)
summary(mod2cCMP)


###### dispersion tests -----
res2cCMP <- simulateResiduals(mod2cCMP)

(t2cCMPSim <- testDispersion(res2cCMP)) # simulation-based residual variance
#pearson residuals are not implemented in truncated_compois models


tab2cCMP <- data.frame(Model = rep("Zero-Truncated CMP", 3),
                    Test = c("Sim-based residual variance", 
                             "Parametric Pearson residuals", 
                             "Nonparametric Pearson Residuals"),
                    Dispersion.parameter = round(c(t2cCMPSim$statistic, 
                                                   NA,
                                                   NA),3),
                    P.value = round(c(t2cCMPSim$p.value, NA, 
                                      NA),2))


#### RESULTS and tables #####


bind_rows(list(m2c=tab2c, m2cCMP=tab2cCMP),.id="model")


AIC(mod2cCMP) - AIC(mod2c)


