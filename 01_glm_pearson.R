### Dispersion Tests Project
## Melina Leite
# Nov 24

library(DHARMa)
library(tidyverse)
library(here)

#####################################
##### Instructions & Simulating #####
#####################################

# 1) Simulating 1000 Poisson datasets with diferent parameters (see below)
# 2) fitting them to correct GLM model -> 
# 3) calculating pearson statistics of the pearson residuals (sum(resˆ2)) -> 
# 4) comparing distribution of these statistics with the Chisquared distribution with the same DF.

#Simulting Poisson:  
#  - varying N 10,20,30,40,50,100,500. 
#  - varying Lambda/intercept 0,1,2,4,8,16,32. 

# function to varying N
calculateStatistics <- function(control = 10){
  testData <- DHARMa::createData(sampleSize = control,
                                 numGroups = 1,
                                 family = poisson())
  
  fittedModel <- stats::glm(observedResponse ~ Environment1, 
                            data = testData, family = poisson()) 
  
  pearson <- residuals(fittedModel, "pearson")
  
  out <- list()
  out$pearson.stat <- sum(pearson^2)
  out$rdf <- df.residual(fittedModel)
  out$disp <- out$pearson.stat/out$rdf
  out$p.value <- DHARMa::testDispersion(fittedModel, 
                                        type="PearsonChisq")$p.value
  return(out)
  }

# varying intercept in a loop (don't know how to have 2 pars in runBenchmarks)
intercept <- c(0,1,2,4,8,16,32)

result <- list()
out.out <- list()
for (i in intercept){
  sampleSize = c(10,20,30,40,50,100,500)
  
  out <- runBenchmarks(calculateStatistics, controlValues = sampleSize,
                       intercept = i,
                       nRep=1000, parallel = 7)
 
  chistats <- as.data.frame(apply(as.data.frame(out$simulations), 2, unlist ))
  chistats$intercept <- i
  result <- rbind(result,chistats)
  out.out[[length(out.out) + 1]] <- out
}


# saving sim results
save(result,out.out,sampleSize,intercept, file=here("data", 
                                                    "glm_pearson_restest.Rdata"))

#######################################
##### Figures and summary results #####
#######################################

#### KS test #### 

ks.p <-list()
ks.d <- list()
for (i in sampleSize){
  for (n in intercept){
    dat <- result[result$controlValues == i & result$intercept == n,]
    ks <- ks.test(dat$pearson.stat,rchisq(100000,df=i-2))
    ks.p <- rbind(ks.p,ks$p.value)
    ks.d <- rbind(ks.d,ks$statistic)
  }
}

res.ks <- data.frame(sampleSize = rep(sampleSize, each=7),
                     intercept = rep(intercept),
                     ks.pvalue = unlist(ks.p),
                     ks.d = unlist(ks.d))
res.ks$col = ifelse(res.ks$ks.pvalue<0.05,"#F7C6C6","#B9DBFA" )


ggplot(res.ks, aes(y=rev(as.factor(sampleSize)),x=as.factor(intercept), fill=col)) +
  geom_tile(show.legend = F) +
  scale_y_discrete("sampleSize", labels=c(500,100,50,40,30,20,10)) +
  scale_x_discrete("intercept", position = "top")+
  geom_text(aes(label=round(ks.pvalue,3))) +
  scale_fill_manual(values= c("#B9DBFA","#F7C6C6")) +
  ggtitle("P-values of KS test for pearson Stats and chisq distr")
ggsave(here("figures", "1_glm_pearson_ksP-values.jpeg"))


#### Distribution plots save as PDF ####


pdf(here("figures", "1_glm_pearson_distributions.pdf"), width=15, height=15)
par(mfrow=c(length(sampleSize),length(intercept)))

for(i in sampleSize){
  for (n in intercept){
   figs2 <- result[result$controlValues == i & result$intercept == n,]
   
   ks <- res.ks[res.ks$sampleSize == i & res.ks$intercept == n,]
   
   plot(density(figs2$pearson.stat), ann=F, axes=F)
  if(ks$ks.pvalue<0.05) rect(par("usr")[1],par("usr")[3],par("usr")[2],
                           par("usr")[4], col = "#F7C6C6")
  par(new=T)
  plot(density(figs2$pearson.stat), 
       main = paste0("sampleSize=",i, "; intercept=",n),
       xlab = paste0("KS: D=", round(ks$ks.d,3), 
                     ", p=", round(ks$ks.pvalue,3)))
  lines(density(rchisq(100000,df=i-2)),col=2)
  }
}
dev.off()

# caption: Black line: Distribution of the Pearson statistics (sum of squares of the Pearson Residuals) from 1000 simulated Poisson data fitted to glm Poisson models. Red line: 100,000 sampled values from a Chisquared distribution with the same residuals degrees of freedm of the models. Red panels are the ones with KS significative results.  

#### Type I error rate for the dispersion tests #####


prop.sig <- list()
for (i in 1:length(out.out)) {
 prop.sig <- rbind(prop.sig, out.out[[i]]$summaries$propSignificant)
}

prop.sig$sampleSize = prop.sig$controlValues
prop.sig$intercept = rep(intercept, each=7)

ggplot(prop.sig, aes(y=rev(as.factor(sampleSize)),x=as.factor(intercept), 
                     fill=p.value)) +
  geom_tile() +
  scale_y_discrete("sampleSize", labels=c(500,100,50,40,30,20,10)) +
  scale_x_discrete("intercept", position = "top")+
  geom_text(aes(label=round(p.value,3))) +
  scale_fill_gradient2(name="type I", low="white", high=2, midpoint=0.05) +
  ggtitle("Proportion of significant Pvalues for Pearson-Chisquare dispersion test")
ggsave(here("figures", "1_glm_pearson_typeI_dispersion.jpeg"))


dat <- left_join(res.ks, prop.sig, by=c("sampleSize","intercept" ))

##### Relationship type 1 error and KS significance result ####

ggplot(dat, aes(x=ks.pvalue, y=p.value, col=as.factor(sampleSize))) + 
  geom_text(aes(label=intercept)) +
  geom_hline(yintercept = 0.05, linetype="dashed") +
  geom_vline(xintercept = 0.05, linetype="dashed") +
  ylab("Type 1 error for dispersion test") +
  ggtitle("Comparing type I dispersion test with KS test",
          subtitle= "Number are intercept values")
ggsave(here("figures", "1_glm_pearson_typeI_ks.jpeg"))

ggplot(dat, aes(x=ks.pvalue, y=p.value, col=as.factor(intercept))) + 
  geom_text(aes(label=sampleSize)) +
  geom_hline(yintercept = 0.05, linetype="dashed") +
  geom_vline(xintercept = 0.05, linetype="dashed") +
  ylab("Type 1 error for dispersion test") 








