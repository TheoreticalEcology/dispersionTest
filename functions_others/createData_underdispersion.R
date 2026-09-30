## createData_ud(): DHARMa::createData() (DHARMa master, GitHub) extended with
## (a) underdispersion for Poisson and binomial responses
##     new arguments: underdispersion (u in [0,1)), underdispersionModel
## (b) mean-preserving overdispersion with a constant dispersion factor phi
##     (negative binomial for Poisson, beta-binomial for binomial; Rev. 2, Comment 2)
##     new arguments: phi (>= 1), overdispersionModel
## Everything else is identical to DHARMa::createData(). In particular, the
## argument 'overdispersion' is still the DHARMa default mechanism (Gaussian
## noise on the linear predictor), which also shifts the marginal mean.

## Samplers for underdispersed counts / proportions with a FIXED MEAN
## (mean = family$linkinv(eta), so the mean structure is untouched)

## 1) Conway-Maxwell-Poisson, mean parametrisation (as in glmmTMB::compois)
##    nu > 1 -> underdispersion; Var ~ mu/nu for moderate/large mu
rcmp_mu <- function(n, mu, nu, maxY = NULL) {
  mu <- rep_len(mu, n)
  if (nu == 1) return(rpois(n, mu))
  if (is.null(maxY)) maxY <- ceiling(max(mu) + 20 * sqrt(max(mu)/nu) + 30)
  y <- 0:maxY
  lfy <- lfactorial(y)
  cmp_mean <- function(loglam) {           # mean of CMP(lambda, nu)
    lp <- y * loglam - nu * lfy
    p <- exp(lp - max(lp)); sum(y * p) / sum(p)
  }
  ## lookup table loglambda -> mean, then invert by interpolation (fast)
  grid <- seq(-30, nu * log(maxY + 1) + 5, length.out = 4000)
  mgrid <- vapply(grid, cmp_mean, 0)
  keep <- !duplicated(mgrid)
  loglam <- approx(mgrid[keep], grid[keep], xout = mu, rule = 2)$y
  out <- integer(n)
  for (i in seq_len(n)) {
    lp <- y * loglam[i] - nu * lfy
    p <- exp(lp - max(lp))
    out[i] <- sample.int(maxY + 1, 1, prob = p) - 1L
  }
  out
}

## 2) Conway-Maxwell-binomial (Shmueli et al. 2005; Kadane 2016)
##    P(y) ∝ choose(N,y)^nu p*^y (1-p*)^(N-y); p* solved so that E(y) = N*prob
##    nu > 1 -> underdispersion
rcmb_mu <- function(n, prob, size, nu) {
  prob <- rep_len(prob, n)
  if (nu == 1) return(rbinom(n, size, prob))
  y <- 0:size
  lc <- nu * lchoose(size, y)
  cmb_mean <- function(lo) { lp <- lc + y * lo; p <- exp(lp - max(lp)); sum(y*p)/sum(p) }
  grid <- seq(-40, 40, length.out = 4000)
  mgrid <- vapply(grid, cmb_mean, 0) / size
  keep <- !duplicated(mgrid)
  lo <- approx(mgrid[keep], grid[keep], xout = prob, rule = 2)$y
  out <- integer(n)
  for (i in seq_len(n)) {
    lp <- lc + y * lo[i]; p <- exp(lp - max(lp))
    out[i] <- sample.int(size + 1, 1, prob = p) - 1L
  }
  out
}

## 3) Hypergeometric "sampling without replacement" binomial:
##    N trials drawn from a finite unit of M >= N items, K = round(M*prob) successes
##    Var = N p(1-p) (M-N)/(M-1)  -> constant dispersion factor phi = (M-N)/(M-1)
rhyperbin <- function(n, prob, size, phi) {
  M <- max(size + 1, round((size - phi) / (1 - phi)))   # solves phi = (M-N)/(M-1)
  Mp <- M * rep_len(prob, n)
  K <- floor(Mp) + rbinom(n, 1, Mp - floor(Mp))          # stochastic rounding keeps E(y) = N*prob
  rhyper(n, m = K, n = M - K, k = size)
}


## Samplers for OVERdispersed counts / proportions with a FIXED MEAN and a
## CONSTANT dispersion factor phi: Var(y|eta) = phi * V(mu) exactly, with
## mu = family$linkinv(eta). phi is therefore the same for every observation,
## intercept and slope (the manuscript's Var(Y) = phi * V(mu)).

## 4) Negative binomial with constant dispersion factor ("NB1", quasi-Poisson
##    variance function). Gamma-Poisson mixture:
##    lambda ~ Gamma(shape = k, rate = k/mu), y ~ Poisson(lambda), k = mu/(phi - 1)
##    -> E(y) = mu, Var(y) = mu + mu^2/k = phi * mu
rnb_phi <- function(n, mu, phi) {
  mu <- rep_len(mu, n)
  if (phi == 1) return(rpois(n, mu))
  rnbinom(n, size = mu / (phi - 1), mu = mu)
}

## 5) Beta-binomial: p_i ~ Beta(a, b) with mean prob and intra-class correlation
##    rho = 1/(a + b + 1) = (phi - 1)/(N - 1), y ~ Binomial(N, p_i)
##    -> E(y) = N prob, Var(y) = N prob (1 - prob) [1 + (N - 1) rho] = phi * N prob (1 - prob)
##    Requires N > 1 and 1 <= phi < N (phi = N would be rho = 1, i.e. all-or-nothing).
rbetabin_phi <- function(n, prob, size, phi) {
  prob <- rep_len(prob, n)
  if (phi == 1) return(rbinom(n, size, prob))
  if (size < 2) stop("beta-binomial overdispersion needs binomialTrials > 1")
  if (phi >= size) stop("phi must be smaller than binomialTrials for the beta-binomial")
  rho <- (phi - 1) / (size - 1)
  s <- 1 / rho - 1                                    # a + b
  p <- rbeta(n, shape1 = prob * s, shape2 = (1 - prob) * s)
  rbinom(n, size, p)
}


#' Simulate test data
#' @description This function creates synthetic dataset with various problems such as overdispersion, zero-inflation, etc.
#' @param sampleSize sample size of the dataset.
#' @param intercept intercept (linear scale).
#' @param fixedEffects vector of fixed effects (linear scale).
#' @param quadraticFixedEffects vector of quadratic fixed effects (linear scale).
#' @param numGroups number of groups for the random effect.
#' @param randomEffectVariance variance of the random effect (intercept).
#' @param overdispersion if this is a numeric value, it will be used as the sd of a random normal variate that is added to the linear predictor. Alternatively, a random function can be provided that takes as input the linear predictor. (DHARMa default mechanism; note that it also changes the marginal mean, e.g. E(y) = exp(eta + sd^2/2) for the Poisson.)
#' @param family a family function for the error distribution and link function to be used in the model to simulate data from. (See [stats::family()] for details of family functions for GLMs.)
#' @param scale scale if the distribution has a scale (e.g. sd for the Gaussian).
#' @param cor correlation between predictors.
#' @param roundPoissonVariance if set, this creates a uniform noise on the Poisson response. The aim of this is to create heteroscedasticity.
#' @param pZeroInflation probability to set any data point to zero.
#' @param binomialTrials number of trials for the binomial. Only active if family == binomial.
#' @param temporalAutocorrelation strength of temporal autocorrelation.
#' @param spatialAutocorrelation strength of spatial autocorrelation.
#' @param factorResponse should the response be transformed to a factor (intended to be used for 0/1 data).
#' @param underdispersion strength of underdispersion, u in [0, 1). 0 = no underdispersion. The data are generated so that the conditional mean is unchanged (= family$linkinv(eta)) and Var(y|eta) is approx. (1 - u) times the Poisson / binomial variance (the approximation is exact for "hypergeometric", and good for CMP/CMB except at very small means or probabilities close to 0/1, where strong underdispersion is mathematically impossible).
#' @param underdispersionModel data-generating mechanism for underdispersion: "CMP" (Conway-Maxwell-Poisson, mean parametrisation, nu = 1/(1-u); Poisson only), "CMB" (Conway-Maxwell-binomial, nu = 1/(1-u); binomial only) or "hypergeometric" (binomial only; trials drawn without replacement from a finite unit, phi = 1-u). Defaults to "CMP" for Poisson and "CMB" for binomial.
#' @param phi dispersion factor for mean-preserving overdispersion, phi >= 1. 1 = no overdispersion. The data are generated so that the conditional mean is unchanged (= family$linkinv(eta)) and Var(y|eta) = phi times the Poisson / binomial variance, exactly and for every observation. Cannot be combined with underdispersion > 0.
#' @param overdispersionModel data-generating mechanism used when phi > 1: "NB" (negative binomial with constant dispersion factor, size = mu/(phi-1); Poisson only) or "betabinomial" (beta-binomial with intra-class correlation rho = (phi-1)/(binomialTrials-1); binomial only, requires binomialTrials > 1 and phi < binomialTrials). Defaults to "NB" for Poisson and "betabinomial" for binomial.
#' @param replicates number of datasets to create.
#' @param hasNA should an NA be added to the environmental predictor (for test purposes).
#' @export
#' @example /inst/examples/createDataHelp.R
createData_ud <- function(sampleSize = 100, intercept = 0, fixedEffects = 1,
                       quadraticFixedEffects = NULL, numGroups = 10,
                       randomEffectVariance = 1, overdispersion = 0,
                       family = poisson(), scale = 1, cor = 0,
                       roundPoissonVariance = NULL,  pZeroInflation = 0,
                       binomialTrials = 1, temporalAutocorrelation = 0,
                       spatialAutocorrelation = 0, factorResponse = FALSE,
                       replicates = 1, hasNA = FALSE,
                       underdispersion = 0, underdispersionModel = NULL,
                       phi = 1, overdispersionModel = NULL){

  if (underdispersion < 0 || underdispersion >= 1) stop("underdispersion must be in [0, 1)")
  if (underdispersion > 0 && !(family$family %in% c("poisson", "binomial"))) stop("underdispersion only implemented for poisson and binomial")
  if (is.null(underdispersionModel)) underdispersionModel = if (family$family == "poisson") "CMP" else "CMB"

  if (phi < 1) stop("phi must be >= 1 (use 'underdispersion' for dispersion factors < 1)")
  if (phi > 1 && underdispersion > 0) stop("use either underdispersion > 0 or phi > 1, not both")
  if (phi > 1 && !(family$family %in% c("poisson", "binomial"))) stop("phi > 1 only implemented for poisson and binomial")
  if (is.null(overdispersionModel)) overdispersionModel = if (family$family == "poisson") "NB" else "betabinomial"


  nPredictors = length(fixedEffects)

  out = list()

  time = sample.int(sampleSize) #change to random order because of issue #436
  x = runif(sampleSize)
  y = runif(sampleSize)

  for (i in 1:replicates){

    ########################################################################
    # Create predictors

    predictors = matrix(runif(nPredictors*sampleSize, min = -1), ncol = nPredictors)

    if (cor != 0){
      predTemp <- runif(sampleSize, min = -1)
      predictors  = (1-cor) * predictors + cor * matrix(rep(predTemp, nPredictors), ncol = nPredictors)
    }

    colnames(predictors) = paste("Environment", 1:nPredictors, sep = "")

    ########################################################################
    # Create random effects

    group = rep(1:numGroups, each = sampleSize/numGroups)
    groupRandom = rnorm(numGroups, sd = sqrt(randomEffectVariance))

    ########################################################################
    # Creation of linear prediction

    linearResponse = intercept + predictors %*% fixedEffects + groupRandom[group]

    if(!is.null(quadraticFixedEffects)){
      linearResponse = linearResponse + predictors^2 %*% quadraticFixedEffects
    }

    ########################################################################
    # Overdispersion on linear predictor


    if(is.numeric(overdispersion)) linearResponse = linearResponse + rnorm(sampleSize, sd = overdispersion)
    if(is.function(overdispersion)) linearResponse = linearResponse + overdispersion(linearResponse)


    ########################################################################
    # Autocorrelation

    if(!(temporalAutocorrelation == 0)){
      distMat <- as.matrix(dist(time))

      invDistMat <- 1/distMat * 5000
      diag(invDistMat) <- 0
      invDistMat = sfsmisc::posdefify(invDistMat)

      temporalError <- MASS::mvrnorm(n = 1, mu = rep(0,sampleSize), Sigma = invDistMat)

      linearResponse = linearResponse + temporalAutocorrelation * temporalError
    }


    if(!(spatialAutocorrelation == 0)) {
      distMat <- as.matrix(dist(cbind(x, y)))

      invDistMat <- 1/distMat * 5000
      diag(invDistMat) <- 0
      invDistMat = sfsmisc::posdefify(invDistMat)

      spatialError <- MASS::mvrnorm(n = 1, mu = rep(0,sampleSize), Sigma = invDistMat)

      linearResponse = linearResponse + spatialAutocorrelation * spatialError
    }


    ########################################################################
    # Link and distribution

    linkResponse = family$linkinv(linearResponse)

    if (family$family == "gaussian") observedResponse = rnorm(n = sampleSize, mean = linkResponse, sd = scale)
    # need checking else if (family$family == "gamma") observedResponse = rgamma(n = sampleSize, shape = linkResponse / scale, scale = scale)
    else if (family$family == "binomial"){
      if (phi > 1) {
        if (overdispersionModel != "betabinomial") stop("overdispersionModel for binomial must be 'betabinomial'")
        observedResponse = rbetabin_phi(sampleSize, linkResponse, binomialTrials, phi = phi)
      }
      else if (underdispersion == 0) observedResponse = rbinom(n = sampleSize, binomialTrials, prob = linkResponse)
      else if (underdispersionModel == "CMB") observedResponse = rcmb_mu(sampleSize, linkResponse, binomialTrials, nu = 1/(1 - underdispersion))
      else if (underdispersionModel == "hypergeometric") observedResponse = rhyperbin(sampleSize, linkResponse, binomialTrials, phi = 1 - underdispersion)
      else stop("underdispersionModel for binomial must be 'CMB' or 'hypergeometric'")
      if (binomialTrials > 1) observedResponse = cbind(observedResponse1 = observedResponse, observedResponse0 = binomialTrials - observedResponse)
    }
    else if (family$family == "poisson") {
      if (phi > 1) {
        if (overdispersionModel != "NB") stop("overdispersionModel for poisson must be 'NB'")
        observedResponse = rnb_phi(sampleSize, linkResponse, phi = phi)
      }
      else if (underdispersion > 0) {
        if (underdispersionModel != "CMP") stop("underdispersionModel for poisson must be 'CMP'")
        observedResponse = rcmp_mu(sampleSize, linkResponse, nu = 1/(1 - underdispersion))
      }
      else if(is.null(roundPoissonVariance)) observedResponse = rpois(n = sampleSize, lambda = linkResponse)
      else observedResponse = round(rnorm(n = length(linkResponse), mean = linkResponse, sd = roundPoissonVariance))
    }
    else if (grepl("Negative Binomial",family$family)) {
      theta = as.numeric(gsub("[\\(\\)]", "", regmatches(family$family, gregexpr("\\(.*?\\)", family$family))[[1]]))
      observedResponse = MASS::rnegbin(linkResponse, theta = theta)
    }
    else stop("wrong link argument supplied")

    ########################################################################
    # Zero-inflation

    if(pZeroInflation != 0){
      artificialZeros = rbinom(n = length(observedResponse), size = 1, prob = 1-pZeroInflation)
      observedResponse = observedResponse * artificialZeros
    }


    if(factorResponse) observedResponse = factor(observedResponse)

    # add spatialError?

    out[[i]] <- data.frame(ID = 1:sampleSize, 
                           observedResponse, 
                           predictors, 
                           group = as.factor(group), 
                           time, 
                           x, y, 
                           expectedMean = linkResponse)
  }
  
  if(length(out) == 1) out = out[[1]]

  if(hasNA) out[1,3] = NA

  return(out)
}
#createData()
