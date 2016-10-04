## distribution functions of the multi GEV

dMultiGEV <- function(x, dur, par=c(0, 1, 0), 
                      parFun=function(dur, par) par,
                      log=FALSE, oneColMat = FALSE) {
  stopifnot(is.function(parFun))
  
  if (is.matrix(x)) {
    stopifnot(ncol(x) == length(dur) | length(dur) == 1)
    res <- x
    dur <- rep(dur, length.out = ncol(res))
    for (i in 1:length(dur)) {
      pars <- parFun(dur[i], par)
      res[,i] <- dgev(x[,i], pars[1], pars[2], pars[3], log)
    }
    colnames(res) <- paste("d", dur, sep="")
  } else {
    if (length(dur) == 1 & !oneColMat) {
      pars <- parFun(dur, par)
      return(dgev(x, pars[1], pars[2], pars[3], log))
    } else {
      res  <- matrix(NA, length(x), length(dur))
      for (i in 1:length(dur)) {
        pars <- parFun(dur[i], par)
        res[,i] <- dgev(x, pars[1], pars[2], pars[3], log)
    }
      colnames(res) <- paste("d", dur, sep="")
    }
  }
  
  return(res)
}

pMultiGEV <- function(q, dur, par=c(0, 1, 0), 
                      parFun=function(dur, par) par,
                      oneColMat = FALSE) {
  stopifnot(is.function(parFun))
  
  if (is.matrix(q)) {
    stopifnot(ncol(q) == length(dur) | length(dur) == 1)
    res <- q
    dur <- rep(dur, length.out = ncol(res))
    for (i in 1:length(dur)) {
      pars <- parFun(dur[i], par)
      res[,i] <- pgev(q[,i], pars[1], pars[2], pars[3], log)
    }
    colnames(res) <- paste("d", dur, sep="")
  } else {
    if (length(dur) == 1 & !oneColMat) {
      pars <- parFun(dur, par)
      return(pgev(q, pars[1], pars[2], pars[3], log))
    } else {
      res  <- matrix(NA, length(q), length(dur))
      for (i in 1:length(dur)) {
        pars <- parFun(dur[i], par)
        res[,i] <- pgev(q, pars[1], pars[2], pars[3], log)
      }
      colnames(res) <- paste("d", dur, sep="")
    }
  }
  
  return(res)
}

qMultiGEV <- function(p, dur, par=c(0, 1, 0), 
                      parFun=function(dur, par) par,
                      oneColMat = FALSE) {
  stopifnot(is.function(parFun))
  
  if (is.matrix(p)) {
    stopifnot(ncol(p) == length(dur) | length(dur) == 1)
    res <- p
    dur <- rep(dur, length.out = ncol(res))
    for (i in 1:length(dur)) {
      pars <- parFun(dur[i], par)
      res[,i] <- qgev(p[,i], pars[1], pars[2], pars[3])
    }
    colnames(res) <- paste("d", dur, sep="")
  } else {
    if (length(dur) == 1 & !oneColMat) {
      pars <- parFun(dur, par)
      return(qgev(p, pars[1], pars[2], pars[3]))
    } else {
      res  <- matrix(NA, length(p), length(dur))
      for (i in 1:length(dur)) {
        pars <- parFun(dur[i], par)
        res[,i] <- qgev(p, pars[1], pars[2], pars[3])
      }
      colnames(res) <- paste("d", dur, sep="")
    }
  }
  
  return(res)
}

rMultiGEV <- function(n, dur,  par=c(0, 1, 0), 
                      parFun=function(dur, par) par,
                      oneColMat = FALSE) {
  stopifnot(is.function(parFun))

  if (length(dur) == 1 & !oneColMat) {
    pars <- parFun(dur, par)
    return(rgev(n, pars[1], pars[2], pars[3]))
  } else {
    res  <- matrix(NA, n, length(dur))
    for (i in 1:length(dur)) {
      pars <- parFun(dur[i], par)
      res[,i] <- rgev(n, pars[1], pars[2], pars[3])
    }
    colnames(res) <- paste("d", dur, sep="")
  }
  
  return(res)
}

# fitment of the multi GEV
fitMultiGEV <- function(annMax, durLevels, start, parFun, ..., inf.res = 1e6) {
  stopifnot(is.function(parFun))
  nDl <- length(durLevels)
  stopifnot(ncol(annMax) == nDl)
  
  optFun <- function(par) {
    cat(par, "\n")
    res <- -sum(sapply(1:nDl, function(cn) dMultiGEV(annMax[,cn], durLevels[cn], 
                                                     par, parFun, log=TRUE)))
    if(is.infinite(res))
      return(inf.res)
    
    res
  }
  
  optVal <- optim(start, optFun, ...) 
  
  res <- optVal$par
  attr(res, "optim") <- optVal
  attr(res, "parFun") <- parFun
  
  return(res)
}


# trimmed moments
TLMOM2parGEV <- function(data, t1, t2){
  TLM <- TLmoms(data, nmom=3, leftrim=t1, rightrim=t2)
  t3 <- TLM$lambdas[3]/TLM$lambdas[2]
  
  if(t1==0 && t2==1){
    z <- 1/(2+t3) * 10/9 - (2*log(2)-log(3))/(3*log(3)-2*log(4))
    xi <- 8.567394*z - 0.675969*z^2
    
    gamxi <- gamma(xi)
    C3 <- (1/3)^xi
    C2 <- (1/2)^xi
    
    sigma <- 2/3 *TLM$lambdas[2]*1/gamxi * 1/(C3 - 2*C2 + 1)
    mu <- TLM$lambdas[1] - sigma/xi - sigma*gamxi*(C2 - 2)
    
    return(c(mu=mu, sigma=sigma, xi=-xi))
  }
  if(t1==1 && t2==0){
    z <- (t3-40/3) * 9/20 - (log(3)-log(4))/(log(2)-log(3))
    xi <- -9.066941*z - 3.374925*z^2 - 0.303208*z^3
    
    gamxi <- gamma(xi)
    C3 <- (1/3)^xi
    C2 <- (1/2)^xi
    
    sigma <- TLM$lambdas[2] * 2 * 1/gamxi * 1/(3*C2 - 3*C3)
    mu <- TLM$lambdas[1] - sigma/xi + sigma*gamxi*C2
    
    return(c(mu=mu, sigma=sigma, xi=-xi))
  }
  if(t1==0 && t2==2) {
    z <- (t3+5/3) * 6/5 - (3*log(5)-8*log(4)+6*log(3))/(log(4)-3*log(3)+3*log(2))
    xi <- -2.468959*z + 1.130074*z^2 - 0.635912*z^3
    
    gamxi <- gamma(xi)
    C4 <- (1/4)^xi
    C3 <- (1/3)^xi
    C2 <- (1/2)^xi
    
    sigma <- TLM$lambdas[2] * 2 * 1/gamxi * 1/(-4*C4 + 12*C3 - 12*C2 + 4)
    mu <- TLM$lambdas[1] - sigma/xi - sigma*gamxi*(-C3 + 3*C2 - 3)
    
    return(c(mu=mu, sigma=sigma, xi=-xi))
  }
  # and if(t1==2 && t2==0) {
  if(t1==0 && t2==0){
    z <- 2*TLM$lambdas[2]/(TLM$lambdas[3]+3*TLM$lambdas[2])-log(2)/log(3)
    xi <- 7.8590 * z + 2.9554 * z^2
    sigma <- (xi * TLM$lambdas[2])/(gamma(1+xi)*(1-2^(-xi)))
    mu <- TLM$lambdas[1]+sigma*(gamma(1+xi)-1)/xi
    
    return(c(mu=mu, sigma=sigma, xi=-xi))
  }
  if(t1==1 && t2==1){
    z <- 9/20*(t3)+(log(3)-2*log(4)+log(5))/(log(2)-2*log(3)+log(4))
    xi <- 25.31711*z - 91.5507*z^2 + 110.0626*z^3 - 46.5518*z^4
    
    gamxi <- gamma(xi)
    C3 <- (1/3)^xi
    C2 <- (1/2)^xi
    C4 <- (1/4)^xi
    
    sigma <- TLM$lambdas[2]*1/gamxi*1/(3*C2-6*C3+3*C4)
    mu <- TLM$lambdas[1]-sigma/xi-sigma*gamxi*(-3*C2+2*C3)
    
    return(c(mu=mu, sigma=sigma, xi=-xi))
  }
  stop("trimming not implemented")
}