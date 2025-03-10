
# function to approximate the residual degrees of freedom of a glmm
# using the package glmmrBase 
# only for lme4 Poisson

approximateDFpearson <- function(model, data, type=c("naive", "KR", "KR2", "sat"),
                              alternative=c("two.sided", "greater", "less")){
  
  if(class(model)[1] != "glmerMod") stop("model is not a glmerMod from lme4")
  if(model@resp$family$family != "poisson") stop("model is not a Poisson")
  
  if(type == "naive"){
    rdf <- df.residual(model) 
    } else{
    f1 <- glmmrBase::lme4_to_glmmr(formula(model), colnames(data))
    mod <- glmmrBase::Model$new(f1, data=data, family=poisson())
    rdf <- dim(data)[1]- mod$small_sample_correction(type)$dof[1]
    }
  
  rp <- residuals(model, "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  
  if(alternative == "greater") pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  else if (alternative == "less") pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=TRUE)
  else if (alternative == "two.sided") pval <- min(min(pchisq(Pearson.chisq, df=rdf, lower.tail=TRUE), pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)) * 2,1)
  
  out = list()
  out$statistic = prat
  names(out$statistic) = "dispersion"
  out$parameter = rdf
  names(out$parameter) = paste0("df-",type) 
  out$method = "Parametric dispersion test via mean Pearson-chisq statistic"
  out$alternative = alternative
  out$p.value = pval
  class(out) = "htest"
  return(out)
}


# testing
# data <- createData(sampleSize = 500, intercept = 5, numGroups=100) 
# 
# model <- glmer(observedResponse ~ Environment1 + (1|group), data=data,
#                family=poisson())
# 
# approximateDFpearson(model, data, type="naive", alternative = "two.sided")
# approximateDFpearson(model, data, type="sat", alternative = "two.sided")
# approximateDFpearson(model, data, type="KR", alternative = "two.sided")
# approximateDFpearson(model, data, type="KR2", alternative = "two.sided")
# 
#  
 
 
 
 
 
 
 
 
 
 
 
 