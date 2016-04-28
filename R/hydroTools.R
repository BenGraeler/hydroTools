library(caTools)
library(evd)
library(xts)

durationLevels <- function(series, timeRes=5/60, durLevels=c(1,3,6,12,24), NaN.rm=FALSE) {
  if(any(durLevels/timeRes < ceiling(durLevels/timeRes)))
    warning("The durLevels are not a multiple of the temporal resolution timeRes. Multiplicies are rounded to the next larger integer.")
  
  res <- NULL
  for(i in 1:length(durLevels)) {
    dl <- durLevels[i]
    res[[i]] <- runmean(series, dl/timeRes, align="left", endrule = "NA")*dl/timeRes
  }
  
  if(any(unlist(lapply(res,is.nan)))) {
    if(!NaN.rm)
      warning("NaNs have been produced.")
    else
      for (elem in 1:length(res))
       res[[elem]][is.nan(res[[elem]])] <- NA
  }
  
  res <- as.data.frame(res)
  colnames(res) <- paste("dl",durLevels, sep="")
  return(res)
}

hydroYears <- function(start, end) {
  stopifnot(end - start > 0)
  
  yrs <- c(start, rep((start+1):(end-1),each=2), end)
  apply(matrix(paste(yrs,c("-11-01","-10-31"),sep=""),ncol = 2, byrow = T),1,function(x) paste(x[1], x[2],sep="/"))
}

aggInt <- function(ts, int, FUN, ...) {
  t(sapply(int, function(x) apply(ts[x], 2, FUN, ...)))
}

## multi GEV
dMultiGEV <- function(x, dur, par, durLevels=c(1,3,6,12,24), log=F) {
  nDl <- length(durLevels)
  durLevels <- c(0, durLevels)
  int <- findInterval(dur, durLevels)
  wght <- (dur-durLevels[int])/(diff(durLevels[int+0:1]))
  
  loc <- sum(par[1:(int-1)]) * (1-wght) + sum(par[1:int]) * wght
  scale <- sum(par[(nDl+1):(nDl + int - 1)]) * (1-wght) + sum(par[(nDl+1):(nDl + int)]) * wght
  shape <- sum(par[(2*nDl+1):(2 * nDl + int - 1)]) * (1-wght) + sum(par[(2*nDl+1):(2 * nDl + int)]) * wght
  
  if (int == 1)
    return(dgev(x, par[1], par[nDl+1], par[2*nDl+1], log=log))
  if (int > nDl)
    return(dgev(x, sum(par[1:nDl]), sum(par[(nDl+1):(2*nDl)]), sum(par[(2*nDl+1):(3*nDl)]), log=log))
  
  return(dgev(x, loc, scale, shape, log))
  
#   loc1 <- loc # sum(par[1:(int-1)])
#   scale1 <- sum(par[(nDl+1):(nDl + int - 1)])
#   shape1 <- sum(par[(2*nDl+1):(2 * nDl + int - 1)])
#   
#   loc2 <- loc # loc1 + par[int]
#   scale2 <- scale1 + par[nDl + int]
#   shape2 <- shape1 + par[2 * nDl + int]
#   
#   if(log)
#     return(log((1-wght) * dgev(x, loc1, scale1, shape1) + wght * dgev(x, loc2, scale2, shape2)))
#   else
#     return((1-wght) * dgev(x, loc1, scale1, shape1) + wght * dgev(x, loc2, scale2, shape2))
}

pMultiGEV <- function(x, dur, par, durLevels=c(1,3,6,12,24)) {
  nDl <- length(durLevels)
  durLevels <- c(0, durLevels)
  int <- findInterval(dur, durLevels)
  wght <- (dur-durLevels[int])/(diff(durLevels[int+0:1]))

  loc <- sum(par[1:(int-1)]) * (1-wght) + sum(par[1:int]) * wght
  scale <- sum(par[(nDl+1):(nDl + int - 1)]) * (1-wght) + sum(par[(nDl+1):(nDl + int)]) * wght
  shape <- sum(par[(2*nDl+1):(2 * nDl + int - 1)]) * (1-wght) + sum(par[(2*nDl+1):(2 * nDl + int)]) * wght
  
  if (int == 1)
    return(pgev(x, par[1], par[nDl+1], par[2*nDl+1]))
  if (int > nDl)
    return(pgev(x, sum(par[1:nDl]), sum(par[(nDl+1):(2*nDl)]), sum(par[(2*nDl+1):(3*nDl)])))

  return(pgev(x, loc, scale, shape))
#   loc1 <- loc # sum(par[1:(int-1)])
#   scale1 <- sum(par[(nDl+1):(nDl + int - 1)])
#   shape1 <- sum(par[(2*nDl+1):(2 * nDl + int - 1)])
#   
#   loc2 <- loc # loc1 + par[int]
#   scale2 <- scale1 + par[nDl + int]
#   shape2 <- shape1 + par[2 * nDl + int]
#   
#   return((1-wght) * pgev(x, loc1, scale1, shape1) + wght * pgev(x, loc2, scale2, shape2))
}

qMultiGEV <- function(p, dur, par, durLevels=c(1,3,6,12,24), interval=c(0,500)) {
  nDl <- length(durLevels)
  durLevels <- c(0, durLevels)
  int <- findInterval(dur, durLevels)
  wght <- (dur-durLevels[int])/(diff(durLevels[int+0:1]))
  
  loc <- sum(par[1:(int-1)]) * (1-wght) + sum(par[1:int]) * wght
  scale <- sum(par[(nDl+1):(nDl + int - 1)]) * (1-wght) + sum(par[(nDl+1):(nDl + int)]) * wght
  shape <- sum(par[(2*nDl+1):(2 * nDl + int - 1)]) * (1-wght) + sum(par[(2*nDl+1):(2 * nDl + int)]) * wght
  
  if (int == 1)
    return(qgev(p, par[1], par[nDl+1], par[2*nDl+1]))
  if (int > nDl)
    return(qgev(p, sum(par[1:nDl]), sum(par[(nDl+1):(2*nDl)]), sum(par[(2*nDl+1):(3*nDl)])))
  
  return(qgev(p, loc, scale, shape))
  
  # sapply(p, function(pi) uniroot(function(x) pMultiGEV(x, dur, par, durLevels) - pi, interval)$root)
}

fitMultiGEV <- function(annMax, durLevels, start, lower, gevDl = durLevels){
  stopifnot(length(start)==length(lower))
  
  nDl <- min(ncol(annMax), length(durLevels))
  
  optFun <- function(par) {
    res <- -sum(sapply(1:nDl, function(cn) dMultiGEV(annMax[,cn], durLevels[cn], par, gevDl, log=TRUE)))
    if(is.infinite(res))
      return(1e5)
    else
      res
  }
  
  return(optim(start, optFun, lower=lower, method="L-BFGS-B")$par)
}

rMultiGEV <- function(n, dur, par, durLevels=c(1,3,6,12,24)) {
  nDl <- length(durLevels)
  durLevels <- c(0, durLevels)
  int <- findInterval(dur, durLevels)
  wght <- (dur-durLevels[int])/(diff(durLevels[int+0:1]))
  
  loc <- sum(par[1:(int-1)]) * (1-wght) + sum(par[1:int]) * wght
  scale <- sum(par[(nDl+1):(nDl + int - 1)]) * (1-wght) + sum(par[(nDl+1):(nDl + int)]) * wght
  shape <- sum(par[(2*nDl+1):(2 * nDl + int - 1)]) * (1-wght) + sum(par[(2*nDl+1):(2 * nDl + int)]) * wght
  
  if (int == 1)
    return(rgev(n, par[1], par[nDl+1], par[2*nDl+1]))
  if (int > nDl)
    return(rgev(n, sum(par[1:nDl]), sum(par[(nDl+1):(2*nDl)]), sum(par[(2*nDl+1):(3*nDl)])))
  
  return(rgev(n, loc, scale, shape))
} 

## extract events from time series

# input: ts, lower threshold, time gap, tc/area
# output: data.frame with: time, endTime, Volume, peak intensity, duration

ts2events <- function(ts, var, AE=0, tc=0.5*AE^0.6, gap=max(4*60,2*tc*60), thd=0.01, tunit="mins"){
  if(ncol(ts) > 1 & missing(var)) {
    warning("First column of xts time series is used.")
    ts <- ts[,1]
  }
  if (!missing(var))
    ts <- ts[,var]
  
  aboveThd <- which(ts > thd)
  eStart <- c(aboveThd[1], aboveThd[which(diff(aboveThd) > 1)+1])
  eEnd <- c(aboveThd[which(diff(aboveThd) > 1)],tail(aboveThd,1))
  
  timeTs <- index(ts)
  boolInd <- which(difftime(timeTs[eStart[-1]], timeTs[eEnd[-length(eEnd)]], units = tunit) > gap)
  eStart <- eStart[c(1, boolInd+1)]
  eEnd <- eEnd[c(boolInd,length(eEnd))]
  
  aboveZero <- which(ts > 0)
  eStartZero <- c(aboveZero[1], aboveZero[which(diff(aboveZero) > 1)+1])
  eEndZero <- c(aboveZero[which(diff(aboveZero) > 1)],tail(aboveZero,1))
  
  for (i in 1:length(eStart)) { # i <- 22
    eStart[i] <- min(eStart[i], tail(eStartZero[eStartZero <= eStart[i] & eStartZero > c(0,eEnd)[i]],1))
    eEnd[i] <- max(eEnd[i], eEndZero[eEndZero >= eEnd[i] & eEndZero < c(eStart,length(ts))[i+1]][1], na.rm = T)
  }
  
  res <- data.frame(time=timeTs[eStart], endTime=timeTs[eEnd])
  stats <- do.call(rbind, aggInt(ts, paste(res$time,res$endTime,sep="/"), 
                  function(x) data.frame(volume=sum(x), peakTime=attr(x,"names")[which.max(x)], peak=max(x))))
  res$volume <- stats$volume
  res$peakTime <- as.POSIXct(as.character(stats$peakTime), format="%Y-%m-%d %H:%M:%S", tz=attributes(index(ts))$tzone)
  res$peak <- stats$peak
  res$duration <- difftime(res$endTime,res$time, units = tunit)
  
  return(res)
}
# 
# startTime <- proc.time()
# ebachDf <- ts2events(ebach, "V3")
# dTime <- proc.time() - startTime
# 
# cat(dTime)
# 
# startTime <- proc.time()
# essenDf <- ts2events(essen, "N")
# dTime2 <- proc.time() - startTime
# 
# 
# cat(dTime2)

###################
# # calculate mowing window quantiles
# ts <- cbind(ts, runquantile(ts, w.width, probs, endrule = "NA", align = w.align))


## count quantile exceedances for a given window
exQt <- function(ts, probs = 0.5, qt,
                 w.width = 10, w.align = "right", 
                 exceed = TRUE, relative = TRUE) {
  if(missing(qt))
    qt <- quantile(ts, probs)
  
  # check whether either qt or w.width has length 1
  stopifnot(length(qt) == 1 | length(w.width) == 1)
  
  exQtSingle <- function(.qt, .w.width) {
    # count exceedances
    if(exceed)
      cExc <- apply(embed(ts, .w.width), 1, function(x) sum(x > .qt))
    else
      cExc <- apply(embed(ts, .w.width), 1, function(x) sum(x < .qt))
    
    if(relative)
      cExc <- cExc/.w.width
    
    switch(w.align,
           left = c(cExc, rep(NA, .w.width-1)),
           center = c(rep(NA, floor((.w.width-1)/2)), cExc, rep(NA, ceiling((.w.width-1)/2))),
           right = c(rep(NA, .w.width-1), cExc),
           stop("Alignment method \"", w.align, "\" is not implemented!"))
  }
  
  if(length(w.width) == 1) {
    res <- do.call(cbind, lapply(qt, function(x) exQtSingle(x, w.width)))
    colnames(res) <- names(qt)
  } else {
    res <- do.call(cbind, lapply(w.width, function(x) exQtSingle(qt, x)))
    colnames(res) <- as.character(w.width)
  }
  
  return(xts(res, index(ts)))
}
