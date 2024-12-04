# changing the runBenchmarks to accept more than ONE parameter to change.


runBenchmarks2 <- function(calculateStatistics, controlValues = NULL, nRep = 10, alpha = 0.05, parallel = FALSE, exportGlobal = FALSE,
                           summaries = FALSE, ...){
  
  
  start_time <- Sys.time()
  
  same.names <- identical(names(formals(calculateStatistics)), names(controlValues))
  
  if(is.list(controlValues) & same.names){
      
      controlValues <- expand.grid(controlValues)
      
      argsims <- list()
      
      for(k in 1:nrow(controlValues)){
        
        argus <- paste0(names(formals(calculateStatistics))[1], "=", 
                        controlValues[k,1])
        
        for (i in 2:length(formals(calculateStatistics))){
          argus <- paste0(argus, ",", names(formals(calculateStatistics))[i], "=",
                          controlValues[k,i])
        }
        argsims[[k]] <-argus
      }
  } else{ stop("controlValues should be a list of named vectors in the same order of the arguments of the function")}
    
  
  # Sequential Simulations
  simulations = list()
  
  if(parallel == FALSE){
    if(is.null(controlValues)) simulations[[1]] = replicate(nRep, calculateStatistics(), simplify = "array")
    else for(j in 1:nrow(controlValues)){
      simulations[[j]] = replicate(nRep, 
                                   eval(parse(text=paste0("calculateStatistics(",
                                                          argsims[[j]],")"))),
                                   simplify = "array")
    }
    
    # Parallel Simulations
    
  }else{
    
    if (parallel == TRUE | parallel == "auto"){
      cores <- parallel::detectCores() - 1
      message("parallel, set cores automatically to ", cores)
    } else if (is.numeric(parallel)){
      cores <- parallel
      message("parallel, set number of cores by hand to ", cores)
    } else stop("wrong argument to parallel")
    
    cl <- parallel::makeCluster(cores)
    
    # for each
    # doParallel::registerDoParallel(cl)
    #
    # `%dopar%` <- foreach::`%dopar%`
    #
    # if(is.null(controlValues)) simulations[[1]] =  t(foreach::foreach(i=1:nRep, .packages=c("lme4", "DHARMa"), .combine = rbind) %dopar% calculateStatistics())
    #
    # else for(j in 1:length(controlValues)){
    #   simulations[[j]] = t(foreach::foreach(i=1:nRep, .packages=c("lme4", "DHARMa"), .combine = rbind) %dopar% calculateStatistics(controlValues[j]))
    # }
    #
    # End for each
    
    # doesn't see to work properly
    loadedPackages = (.packages())
    parExectuer = function(x = NULL, control = NULL) calculateStatistics(control)
    if (exportGlobal == TRUE) parallel::clusterExport(cl = cl, varlist = ls(envir = .GlobalEnv))
    parallel::clusterExport(cl = cl, c("parExectuer", "calculateStatistics", "loadedPackages"), 
                            envir = environment())
    parallel::clusterEvalQ(cl, {for(p in loadedPackages) library(p, character.only=TRUE)})
    
    # parallel::clusterExport(cl = cl, varlist = ls(envir = .GlobalEnv))
    
    # parallel::clusterExport(cl=cl,varlist = c("calculateStatistics"), envir=environment())
    # parallel::clusterExport(cl=cl,varlist = c("controlValues", "alpha", envir=environment())
    
    # parallel::clusterExport(cl=cl,varlist = c("calculateStatistics"), envir=environment())
    # parallel::clusterExport(cl=cl,varlist = c("controlValues", "alpha", envir=environment())
    
    if(is.null(controlValues)) simulations[[1]] = parallel::parSapply(cl, 1:nRep, parExectuer)
    else for(j in 1:nrow(controlValues)){
      simulations[[j]] = parallel::parSapply(cl, 1:nRep, parExectuer, control = controlValues[j])
    }
    
    parallel::stopCluster(cl)
  }
  
  # Calculations of summaries
  
  if(is.null(controlValues)) controlValues = c("N")
  
  nOutputs = nrow(simulations[[1]])
  nControl = nrow(controlValues)
  
  # reducing the list of outputs to a data.frame
  x = Reduce(rbind, lapply(simulations, t))
  x = data.frame(x)
  x$replicate = rep(1:nRep, nrow(controlValues))
  contrs <- cbind(controlValues, simID = rep(row.names(controlValues),nRep))
  x <- cbind(contrs,x)
  
  if(summaries==T){
    summary = list()
    
    # function for aggregation
    aggreg <- function(f,...) {
      ret <- aggregate(x[,- c(ncol(x) - 1, ncol(x))], by = list(x$controlValues), f)
      colnames(ret)[1] = "controlValues"
      return(ret)
    }
    
    if(length(is.na(x)>0)){
      warning("NA values in the output, this might be a problem for the summaries.")
    }
    
    sig <- function(x) mean(x < alpha, na.rm=T)
    isUnif <- function(x) ks.test(x, "punif")$p.value
    
    summary$propSignificant = aggreg(sig)
    summary$meanP = aggreg(mean, na.rm=T)
    summary$isUnifP = aggreg(isUnif)
    
  }
  

  out = list()
  out$controlValues = controlValues
  out$simulations = x
  if(summaries==T) out$summaries = summary
  if(summaries==T) out$nSummaries = ncol(x) - 2
  out$time = Sys.time() - start_time

  out$nReplicates = nRep
  
  class(out) = "DHARMaBenchmark"
  
  return(out)
}
