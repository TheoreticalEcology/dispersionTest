### Dispersion Tests Project
# Melina Leite
# Sep 26

## Monte Carlo standard errors (MCSE) for summaries over simulation replicates
# Following Morris, White & Crowther (2019), Stat Med 38:2074-2102.
#
# mcse_mean:   MCSE of the mean   = sd / sqrt(nsim)
# mcse_median: MCSE of the median = bootstrap SD of the median over the replicates
#
# For proportions (type I error, power) the MCSE is reported as the exact
# binomial 95% confidence interval given by binom.test(), computed directly in 
# each script.
#
# In the figures, error bars are 95% intervals throughout: the exact binomial
# CI for proportions, and estimate +/- 1.96 * MCSE for means and medians.

mcse_mean <- function(x) {
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

mcse_median <- function(x, nboot = 1000, seed = 42) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  # fixed seed so the figures are reproducible; the global RNG stream is
  # saved and restored so that nothing else in the script is affected
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = .GlobalEnv)) 
      get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv), 
            add = TRUE)
    set.seed(seed)
  }
  sd(replicate(nboot, median(sample(x, replace = TRUE))))
}
