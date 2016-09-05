## distribution functions of the multi GEV

# hidden heloer function
.setupParaMultiGEV <- function(dur, par, durLevels, nDl, int) {
  wght <- (dur-durLevels[int])/(sapply(int, function(x) diff(durLevels[x+0:1])))
  
  loc <- sapply(1:length(int),
                function(x) sum(par[1:(int[x]-1)]) * (1-wght[x]) + sum(par[1:int[x]]) * wght[x])
  scale <- sapply(1:length(int),
                  function(x) sum(par[(nDl+1):(nDl + int[x] - 1)]) * (1-wght[x]) + sum(par[(nDl+1):(nDl + int[x])]) * wght[x])
  shape <- sapply(1:length(int),
                  function(x) sum(par[(2*nDl+1):(2 * nDl + int[x] - 1)]) * (1-wght[x]) + sum(par[(2*nDl+1):(2 * nDl + int[x])]) * wght[x])
  
  return(list(loc=loc, scale=scale, shape=shape))
}

## density of the multi GEV
dMultiGEV <- function(x, dur, par, log=F) {
  durLevels <- attr(par, "gevDl")
  nDl <- length(durLevels)
  durLevels <- c(0, durLevels)
  int <- findInterval(dur, durLevels)
  
  parList <- .setupParaMultiGEV(dur, par, durLevels, nDl, int)
  
  if (length(dur) == 1) {
    if (int == 1)
      return(dgev(x, par[1], par[nDl+1], par[2*nDl+1], log = log))
    
    if (int > nDl)
      return(dgev(x, sum(par[1:nDl]), sum(par[(nDl+1):(2*nDl)]), 
                  sum(par[(2*nDl+1):(3*nDl)]), log = log))
    
    return(dgev(x, parList$loc, parList$scale, parList$shape, log = log))
  } else {
    res <- matrix(NA, nrow = length(x), ncol = length(int))
    
    for (i in 1:length(int)) {
      if (int[i] == 1)
        res[,i] <- dgev(x, par[1], par[nDl+1], par[2*nDl+1], log = log)
      
      if (int[i] > nDl)
        res[,i] <- dgev(x, sum(par[1:nDl]), sum(par[(nDl+1):(2*nDl)]),
                        sum(par[(2*nDl+1):(3*nDl)]), log = log)
      
      if (int[i] > 1 & int[i] <= nDl)
        res[,i] <- dgev(x, parList$loc[i], parList$scale[i], parList$shape[i], 
                        log = log)
    }
    
    colnames(res) <- paste("d", dur, sep="")
  }
  
  return(res)
}

# CDF of the multi GEV
pMultiGEV <- function(q, dur, par) {
  durLevels <- attr(par, "gevDl")
  nDl <- length(durLevels)
  durLevels <- c(0, durLevels)
  int <- findInterval(dur, durLevels)
  
  parList <- .setupParaMultiGEV(dur, par, durLevels, nDl, int)
  
  if (length(dur) == 1) {
    if (int == 1)
      return(pgev(q, par[1], par[nDl+1], par[2*nDl+1]))
    
    if (int > nDl)
      return(pgev(q, sum(par[1:nDl]), sum(par[(nDl+1):(2*nDl)]), sum(par[(2*nDl+1):(3*nDl)])))
    
    return(pgev(q, parList$loc, parList$scale, parList$shape))
  } else {
    res <- matrix(NA, nrow = length(q), ncol = length(int))
    
    for (i in 1:length(int)) {
      if (int[i] == 1)
        res[,i] <- pgev(q, par[1], par[nDl+1], par[2*nDl+1])
      
      if (int[i] > nDl)
        res[,i] <- pgev(q, sum(par[1:nDl]), sum(par[(nDl+1):(2*nDl)]), sum(par[(2*nDl+1):(3*nDl)]))
      
      if (int[i] > 1 & int[i] <= nDl)
        res[,i] <- pgev(q, parList$loc[i], parList$scale[i], parList$shape[i])
    }
    
    colnames(res) <- paste("d", dur, sep="")
    return(res)
  }
}

# qunatile function of multi GEV
qMultiGEV <- function(p, dur, par, interval=c(0,500)) {
  durLevels <- attr(par, "gevDl")
  nDl <- length(durLevels)
  durLevels <- c(0, durLevels)
  int <- findInterval(dur, durLevels)
  
  parList <- .setupParaMultiGEV(dur, par, durLevels, nDl, int)
  
  if (length(dur) == 1) {
    if (int == 1)
      return(qgev(p, par[1], par[nDl+1], par[2*nDl+1]))
    if (int > nDl)
      return(qgev(p, sum(par[1:nDl]), sum(par[(nDl+1):(2*nDl)]), sum(par[(2*nDl+1):(3*nDl)])))
    
    return(qgev(p, parList$loc, parList$scale, parList$shape))
  } else {
    res <- matrix(NA, nrow = length(p), ncol = length(int))
    
    for (i in 1:length(int)) {
      if (int[i] == 1)
        res[,i] <- qgev(p, par[1], par[nDl+1], par[2*nDl+1])
      
      if (int[i] > nDl)
        res[,i] <- qgev(p, sum(par[1:nDl]), sum(par[(nDl+1):(2*nDl)]), sum(par[(2*nDl+1):(3*nDl)]))
      
      if (int[i] > 1 & int[i] <= nDl)
        res[,i] <- qgev(p, parList$loc[i], parList$scale[i], parList$shape[i])
    }
    
    colnames(res) <- paste("d", dur, sep="")
    return(res)
  }
}

# fitment of the multi GEV
fitMultiGEV <- function(annMax, durLevels, start, lower, gevDl = durLevels){
  stopifnot(length(start)==length(lower))
  
  nDl <- min(ncol(annMax), length(durLevels))
  
  optFun <- function(par) {
    attr(par, "gevDl") <- gevDl
    res <- -sum(sapply(1:nDl, function(cn) dMultiGEV(annMax[,cn], durLevels[cn], par,log=TRUE)))
    if(is.infinite(res))
      return(1e5)
    else
      res
  }
  
  optVal <- optim(start, optFun, lower=lower, method="L-BFGS-B") 
  
  res <- optVal$par
  attr(res, "optim") <- optVal
  attr(res, "gevDl") <- gevDl
  
  return(res)
}

# random numbers of the multi GEV
rMultiGEV <- function(n, dur, par) {
  durLevels <- attr(par, "gevDl")
  nDl <- length(durLevels)
  durLevels <- c(0, durLevels)
  int <- findInterval(dur, durLevels)
  
  parList <- .setupParaMultiGEV(dur, par, durLevels, nDl, int)
  
  if (length(dur) == 1) {
    if (int == 1)
      return(rgev(n, par[1], par[nDl+1], par[2*nDl+1]))
    
    if (int > nDl)
      return(rgev(n, sum(par[1:nDl]), sum(par[(nDl+1):(2*nDl)]),
                  sum(par[(2*nDl+1):(3*nDl)])))
    
    return(rgev(n, parList$loc, parList$scale, parList$shape))
  } else {
    res <- matrix(NA, nrow = n, ncol = length(int))
    
    for (i in 1:length(int)) {
      if (int[i] == 1)
        res[,i] <- rgev(n, par[1], par[nDl+1], par[2*nDl+1])
      
      if (int[i] > nDl)
        res[,i] <- rgev(n, sum(par[1:nDl]), sum(par[(nDl+1):(2*nDl)]), sum(par[(2*nDl+1):(3*nDl)]))
      
      if (int[i] > 1 & int[i] <= nDl)
        res[,i] <- rgev(n, parList$loc[i], parList$scale[i], parList$shape[i])
    }
    
    colnames(res) <- paste("d", dur, sep="")
    return(res)
  }
}