### Dispersion Tests Project
## Melina Leite
# Nov 24

library(DHARMa)
library(tidyverse)
library(here)

#####################################
##### Instructions & Simulating #####
#####################################

## KS GOF test to compare Pearson Stats and Chi-squared distribution


# 1) Simulating 1000 Poisson datasets with different sample sizes and intercepts
    # -sampleSize: c(10,50,100,500)
    # -intercept:  c(-3,-1,0,2,4)
# 2) fitting them to correct GLM model
# 3) calculating pearson statistics of the pearson residuals (sum(resˆ2)) 
# 4) comparing distribution of these statistics with the Chi-squared distribution with the same DF.
# 5) Repeating these steps 100 times to get the proportion of significant results per combination of parameters.


# function to varying N
calculateStatistics <- function(control = 10){
  testData <- DHARMa::createData(sampleSize = control,
                                 numGroups = 1,
                                 family = poisson())
  
  fittedModel <- stats::glm(observedResponse ~ Environment1, 
                            data = testData, family = poisson()) 
  
  pearson <- residuals(fittedModel, "pearson")
  
  out <- list()
  out$Pear.stat <- sum(residuals(fittedModel, "pearson")^2)
  out$rdf <- df.residual(fittedModel)
  return(unlist(out))
}

# sampleSizes
sampleSize = c(10,50,100,500)

# intercept (in a loop)
intercept <- c(-3,-1,0,2,4)


# KS TEST
ks.p <- function(x, df) ks.test(x,"pchisq", df=df)$p.value

final.res <- list()

for (k in 1:100) { # MANY SIMULATIONS TO HAVE A PROP OF SIG RESULTS
  set.seed(k)  
  result <- list()
    
  for (i in intercept){
    out <- runBenchmarks(calculateStatistics, controlValues = sampleSize,
                         intercept = i,
                         nRep=1000, parallel = 15)

    stats <- as.data.frame(apply(as.data.frame(out$simulations), 2, unlist ))
    
    res <- stats %>% group_by(controlValues) %>%
      summarise(ks.p = ks.p(Pear.stat,rdf))
    res$intercept <- i
    result <- rbind(result,res)
  }
  final.res[[k]] <- result
}

# saving results
save(final.res, file = here("data", "1_glm_pearson_KStests.Rdata"))
#load(here("data", "1_glm_pearson_KStests.Rdata"))


# manipulaing results
names(final.res) <- 1:length(final.res)

final <- bind_rows(final.res, .id="sim") %>% group_by(controlValues,intercept) %>%
  summarise(ks.sig = sum(ks.p<0.05))

for (i in 1:nrow(final)) {
  confs <- binom.test(final$ks.sig[i], 100)$conf.int
  final$conf.low[i] <-  round(confs[1],2)
  final$conf.up[i] <-  round(confs[2],2)
}


### FIGURES

ggplot(final, aes(y=ks.sig/100, x=as.factor(controlValues),
                      col=as.factor(intercept))) +
  geom_point() + geom_line(aes(x=as.numeric(as.factor(controlValues)))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up),width = 0.1)+
  scale_color_discrete("intercept")+
  xlab("sampleSize") + ylab("Prop of significant KS test")

ggplot(final, aes(y=ks.sig/100, col=as.factor(controlValues),
                      x=as.factor(intercept))) +
  geom_point() + geom_line(aes(x=as.numeric(as.factor(intercept)))) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up),width = 0.1)+
  scale_color_discrete("sampleSize")+
  xlab("intercept") + ylab("Prop of significant KS test")

# Final figure

ggplot(final, aes(y=rev(as.factor(controlValues)),
                      x=as.factor(intercept), 
                      fill=ks.sig/100)) +
  geom_tile(show.legend = F) +
  scale_y_discrete("sampleSize", labels=c(500,100,50,10)) +
  scale_x_discrete("intercept", position = "top")+
  geom_text(aes(label=round(ks.sig/100,2))) +
  geom_text(aes(y=as.numeric(rev(as.factor(controlValues)))-0.3, 
                label= paste0("(",conf.low," - ",conf.up,")")), size=3) +
  scale_fill_gradient(low="white", high=2) +
  ggtitle("Prop of significant KS tests for Pearson Stats x Chisq distr",
          sub="100 tests each combination")
ggsave(here("figures", "1_glm_pearson_ksP_propsig.jpeg"))


















